# This script runs PHITS simulations for all generated input files in parallel.
#
# Features:
#   - Use phantom selection from params["phantom"]
#   - Select source organ(s) from SOURCE_CSV
#   - Display source-organ IDs and names
#   - Select all source organs
#   - Run new PHITS simulations
#   - Skip completed simulations
#   - Re-run completed simulations
#   - Automatically determine PHITS parallelization from input files
#   - Run multiple PHITS jobs in parallel

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import pandas as pd

import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom_group

def run_phits(params):

    # ==========================================================
    # Configuration
    # ==========================================================
    
    phits_executable = config.PHITS_EXECUTABLE
    input_root = config.GENERATED_INPUTS_DIR
    source_csv = config.SOURCE_CSV

    # Pattern for OMP/MPI threads
    parallelization_pattern = re.compile(
        r"\$(OMP|MPI)\s*=\s*(\d+)",
        re.IGNORECASE
    )

    # ==========================================================
    # User Interface
    # ==========================================================

    print("\n")
    print("=" * 90)
    print("PHITS SAF Simulation Launcher")
    print("=" * 90)

    phantom_selection = params["phantom"]

    selected_phantoms = get_phantom_group(phantom_selection)
        
    # ==========================================================
    # Source Organ Selection
    # ==========================================================

    if not source_csv.is_file():

        raise FileNotFoundError(
            f"Source organ file not found:\n"
            f"{source_csv}"
        )

    source_df = pd.read_csv(source_csv)

    print("\n")
    print("=" * 90)
    print("Source Organ Selection")
    print("=" * 90)

    print(
        f"Source organ file:\n"
        f"{source_csv}"
    )

    # ------------------------------------------------------
    # Validate the updated source-organ CSV schema
    # ------------------------------------------------------
    # Current source_organs.csv columns:
    #   Source region | Acronym | ID number(s)
    required_columns = {"Source region", "Acronym", "ID number(s)"}
    missing_columns = required_columns.difference(source_df.columns)

    if missing_columns:
        raise ValueError(
            "The source-organ CSV is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}. "
            f"Found columns: {source_df.columns.tolist()}\n"
            f"File: {source_csv}"
        )

    # Clean up labels and IDs while preserving multi-ID groups
    for column in required_columns:
        source_df[column] = source_df[column].astype("string").str.strip()

    source_df = source_df.dropna(subset=["Acronym", "ID number(s)"])
    source_df = source_df[
        source_df["Acronym"].ne("") & source_df["ID number(s)"].ne("")
    ].drop_duplicates(subset=["Acronym"], keep="first").reset_index(drop=True)

    if source_df.empty:
        raise ValueError(f"No valid source-organ rows were found in {source_csv}")

    # ------------------------------------------------------
    # CSV acronyms are the same labels used in generated input filenames.
    # The IDs column may contain multiple IDs, e.g. "9700_9900".
    available_sources = []
    for _, row in source_df.iterrows():
        acronym = str(row["Acronym"]).strip()
        ids_text = str(row["ID number(s)"]).strip()
        source_region = str(row["Source region"]).strip()
        available_sources.append({
            "source_region": source_region,
            "acronym": acronym,
            "ids_text": ids_text,
        })

    # ------------------------------------------------------
    # Display source-organ groups
    # ------------------------------------------------------
    for i, source in enumerate(available_sources, start=1):
        print(
            f"[{i}] {source['source_region']} "
            f"(Acronym: {source['acronym']}; IDs: {source['ids_text']})"
        )

    print("[A] All source organs")

    # ------------------------------------------------------
    # Select source-organ groups by CSV row
    # ------------------------------------------------------
    while True:
        source_choice = input(
            "\\nSelect source organ group(s): "
        ).strip().upper()

        if source_choice == "A":
            selected_sources = available_sources.copy()
            break

        try:
            indices = [int(x.strip()) - 1 for x in source_choice.split(",")]
            if not indices or any(
                index < 0 or index >= len(available_sources)
                for index in indices
            ):
                raise ValueError

            selected_sources = [
                available_sources[index]
                for index in dict.fromkeys(indices)
            ]
            break

        except ValueError:
            print(
                "\\nInvalid selection. Enter valid row numbers "
                "separated by commas, or A."
            )

    print("\\nSelected source organ group(s):")
    for source in selected_sources:
        print(
            f"  - {source['source_region']} "
            f"(Acronym: {source['acronym']}; IDs: {source['ids_text']})"
        )

    selected_source_acronyms = {
        source["acronym"] for source in selected_sources
    }

    # ==========================================================
    # Find input files
    # ==========================================================

    input_files = []

    # The launcher matches the selected CSV acronym directly against
    # the source label in each generated filename. This supports both
    # single-ID sources (e.g. Brain) and multi-ID groups (e.g. Lungs).

    # ------------------------------------------------------
    # Search generated PHITS input files
    # ------------------------------------------------------

    for phantom in selected_phantoms:

        phantom_dir = (
            input_root /
            phantom
        )

        if phantom_dir.exists():
            # Standard layout: one subdirectory per phantom.
            phantom_input_files = sorted(
                phantom_dir.rglob("*.inp")
            )
        else:
            # Also support a flat/nested input directory where files are
            # named with the phantom in the filename, e.g.
            # phits_MRCP_AF_H150W55_source_Lungs_...inp.
            print(
                f"\\nDirectory not found for {phantom}: {phantom_dir}"
            )
            print(
                "Searching the configured input root for matching "
                "phantom filenames instead."
            )
            phantom_token = f"_{phantom}_source_"
            phantom_input_files = sorted(
                infile
                for infile in input_root.rglob("*.inp")
                if phantom_token in infile.name
                or infile.name.startswith(f"{phantom}_source_")
            )

        for infile in phantom_input_files:

            filename = infile.name

            # --------------------------------------------------
            # Match source acronym in filename
            # --------------------------------------------------
            matched = any(
                f"_source_{acronym.replace(',', '').replace(' ', '_')}_"
                in filename
                for acronym in selected_source_acronyms
            )

            if matched:

                input_files.append(
                    infile
                )

    input_files = sorted(
        set(input_files)
    )

    # ------------------------------------------------------
    # Check whether input files were found
    # ------------------------------------------------------

    if not input_files:

        print(
            "\nNo PHITS input files were found "
            "for the selected phantom/source organs."
        )

        return

    # ======================================================
    # Existing simulation check
    # ======================================================

    print("=" * 90)
    print("Existing Simulation Check")
    print("=" * 90)

    print(
        "[1] Skip jobs that have already been simulated"
    )

    print(
        "[2] Re-run completed jobs "
        "(overwrite outputs)"
    )

    while True:

        overwrite_choice = input(
            "Enter your choice (1-2): "
        ).strip()

        if overwrite_choice in {"1", "2"}:
            break

        print(
            "Invalid choice. "
            "Please enter 1 or 2."
        )

    skip_completed = (
        overwrite_choice == "1"
    )

    # ==========================================================
    # Check PHITS job completion
    # ==========================================================

    def job_completed(inp_file):

        job_dir = inp_file.parent

        batch_exists = (job_dir / "batch.out").exists()

        stem = inp_file.stem

        required_outputs = [
            job_dir / f"{stem}.out",
            job_dir / f"phits_deposit_{stem}.out",
            job_dir / f"phits_fluence_{stem}.out",
        ]

        return all(output.is_file() and output.stat().st_size > 0
                   for output in required_outputs)

    # ==========================================================
    # Filter completed jobs
    # ==========================================================

    if skip_completed:

        original = len(input_files)

        input_files = [
            inp
            for inp in input_files
            if not job_completed(inp)
        ]

        skipped = (
            original -
            len(input_files)
        )

        print(
            f"\nSkipping {skipped} "
            f"completed job(s)."
        )

        print(
            f"Remaining jobs: "
            f"{len(input_files)}"
        )

    # ==========================================================
    # Check whether anything remains
    # ==========================================================

    if not input_files:

        print(
            "\nNo simulations to run. "
            "All jobs have already been completed."
        )

        return

    # ==========================================================
    # Determine PHITS threads
    # ==========================================================

    def get_threads(input_file):

        """
        Read the number of OMP/MPI threads
        from a PHITS input file.
        """

        with open(
            input_file,
            "r"
        ) as f:

            for line in f:

                match = (
                    parallelization_pattern
                    .search(line)
                )

                if match:

                    return int(
                        match.group(2)
                    )

        raise ValueError(
            f"No $OMP or $MPI definition "
            f"found in {input_file}"
        )

    # Use the first input file to determine
    # the number of threads per PHITS job.
    threads = get_threads(
        input_files[0]
    )

    # ==========================================================
    # Parallelization
    # ==========================================================

    total_threads = os.cpu_count() or 1

    max_parallel = max(
        1,
        total_threads // threads
    )

    print(
        f"\nFound {len(input_files)} "
        f"input files."
    )

    print(
        f"Total CPU threads : "
        f"{total_threads}"
    )

    print(
        f"Threads per PHITS : "
        f"{threads}"
    )

    print(
        f"Parallel PHITS jobs: "
        f"{max_parallel}"
    )

    print(
        f"\nRunning up to "
        f"{max_parallel} jobs simultaneously.\n"
    )

    # ==========================================================
    # Run one PHITS simulation
    # ==========================================================

    def run_phits_job(input_file):

        print(
            f"Starting {input_file.name}"
        )

        process = subprocess.Popen(
            [
                phits_executable,
                input_file.name
            ],
            cwd=input_file.parent
        )

        returncode = process.wait()

        if returncode == 0:

            print(
                f"Finished {input_file.name}"
            )

        else:

            print(
                f"Failed {input_file.name} "
                f"(return code {returncode})"
            )

        return (input_file, returncode)

    # ==========================================================
    # Run jobs in parallel
    # ==========================================================

    with ThreadPoolExecutor(
        max_workers=max_parallel
    ) as executor:

        futures = [
            executor.submit(
                run_phits_job,
                inp
            )
            for inp in input_files
        ]

        for future in as_completed(
            futures
        ):

            input_file, returncode = (
                future.result()
            )

            if returncode == 0:

                print(
                    f"✓ {input_file.name} "
                    f"completed successfully"
                )

            else:

                print(
                    f"✗ {input_file.name} "
                    f"failed "
                    f"(return code "
                    f"{returncode})"
                )

    print("\nAll PHITS jobs finished.")