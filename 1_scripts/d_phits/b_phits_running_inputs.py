# This script runs PHITS simulations for all generated input files in parallel.
#
# Features:
#   - Use phantom selection from params["phantom"]
#   - Select source organ(s) from SOURCE_CSV
#   - Display source-organ IDs and names
#   - Select all source organs
#   - Run new PHITS simulations
#   - Re-run failed simulations from RERUN_CSV
#   - Skip completed simulations
#   - Re-run completed simulations
#   - Automatically determine PHITS parallelization from input files
#   - Run multiple PHITS jobs in parallel

from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import pandas as pd

import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS

def run_phits(params):

    # ==========================================================
    # Configuration
    # ==========================================================
    
    phits_executable = config.PHITS_EXECUTABLE
    input_root = config.GENERATED_INPUTS_DIR
    rerun_csv = config.PHITS_RERUN_CSV_FILE
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
    print("=" * 50)
    print("PHITS SAF Simulation Launcher")
    print("=" * 50)

    # ----------------------------------------------------------
    # Select operation
    # ----------------------------------------------------------

    print("[1] Run new simulations")
    print("[2] Re-run failed simulations")

    while True:

        choice = input("Enter your choice (1-2): ").strip()

        if choice in {"1", "2"}:
            break

        print("Invalid choice. Please enter 1 or 2.")

    phantom_selection = params["phantom"]

    if phantom_selection == "MRCP_AM":
        selected_phantoms = ["MRCP_AM"]

    elif phantom_selection == "MRCP_AF":
        selected_phantoms = ["MRCP_AF"]

    elif phantom_selection == "MRCP_AF_AM":
        selected_phantoms = ["MRCP_AM", "MRCP_AF"]

    elif phantom_selection == "MFCP_AM":
        selected_phantoms = ["MFCP_AM"]

    elif phantom_selection == "MFCP_AF":
        selected_phantoms = ["MFCP_AF"]

    elif phantom_selection == "MFCP_AF_AM":
        selected_phantoms = ["MFCP_AM", "MFCP_AF"]

    else:
        raise ValueError(f"Unknown phantom selection: {phantom_selection}")
    
    # ==========================================================
    # Source Organ Selection
    # ==========================================================

    if choice == "1":

        if not source_csv.is_file():

            raise FileNotFoundError(
                f"Source organ file not found:\n"
                f"{source_csv}"
            )

        source_df = pd.read_csv(source_csv)

        print("\n")
        print("=" * 50)
        print("Source Organ Selection")
        print("=" * 50)

        print(
            f"Source organ file:\n"
            f"{source_csv}"
        )

        # ------------------------------------------------------
        # Source-organ column
        # ------------------------------------------------------

        source_organ_column = "source_organ_ID"

        if source_organ_column not in source_df.columns:

            raise ValueError(
                f"Column '{source_organ_column}' was not found "
                f"in {source_csv}"
            )

        # ------------------------------------------------------
        # Get source-organ IDs
        # ------------------------------------------------------

        available_source_organs = (
            source_df[source_organ_column]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
            .tolist()
        )

        if not available_source_organs:

            raise ValueError(
                "No source organs were found in "
                f"{source_csv}"
            )

        # ------------------------------------------------------
        # Build source-organ ID → name lookup
        # ------------------------------------------------------

        source_organ_names = {}

        for phantom in selected_phantoms:

            phantom_sources = SOURCE_ORGANS.get(
                phantom,
                {}
            )

            for organ_id, organ_name in phantom_sources.items():

                source_organ_names[str(organ_id)] = organ_name

        # ------------------------------------------------------
        # Display source organs
        # ------------------------------------------------------

        for i, organ_id in enumerate(
            available_source_organs,
            start=1
        ):

            organ_name = source_organ_names.get(
                organ_id,
                "Unknown"
            )

            print(
                f"[{i}] {organ_id} ({organ_name})"
            )

        print("[A] All source organs")

        # ------------------------------------------------------
        # Select source organs
        # ------------------------------------------------------

        while True:

            source_choice = input(
                "\nSelect source organ(s): "
            ).strip().upper()

            if source_choice == "A":

                selected_source_organs = (
                    available_source_organs.copy()
                )

                break

            try:

                indices = [
                    int(x.strip()) - 1
                    for x in source_choice.split(",")
                ]

                if not indices:
                    raise ValueError

                if any(
                    index < 0
                    or index >= len(
                        available_source_organs
                    )
                    for index in indices
                ):
                    raise ValueError

                selected_source_organs = list(
                    dict.fromkeys(
                        available_source_organs[index]
                        for index in indices
                    )
                )

                break

            except ValueError:

                print(
                    "\nInvalid selection. "
                    "Please enter valid numbers "
                    "separated by commas, or A."
                )

        # ------------------------------------------------------
        # Display selected source organs
        # ------------------------------------------------------

        print("\nSelected source organ(s):")

        for organ_id in selected_source_organs:

            organ_name = source_organ_names.get(
                organ_id,
                "Unknown"
            )

            print(
                f"  - {organ_id} ({organ_name})"
            )

    # ==========================================================
    # Find input files
    # ==========================================================

    if choice == "2":

        # ======================================================
        # Re-run failed simulations
        # ======================================================

        if not rerun_csv.exists():

            raise FileNotFoundError(
                f"\nRerun file not found:\n"
                f"{rerun_csv}\n\n"
                "Run the uncertainty-checking script first."
            )

        rerun = pd.read_csv(rerun_csv)

        if rerun.empty:

            print(
                "\nThere are no failed simulations "
                "to rerun."
            )

            return

        if "Input File" not in rerun.columns:

            raise ValueError(
                "The rerun file does not contain "
                "an 'Input File' column.\n"
                "Run the latest uncertainty-checking "
                "script."
            )

        input_files = []

        for relative_path in rerun["Input File"]:

            input_file = (
                input_root /
                Path(relative_path)
            )

            if input_file.exists():

                input_files.append(
                    input_file
                )

            else:

                print(
                    f"Missing input file: "
                    f"{relative_path}"
                )

        input_files = sorted(
            set(input_files)
        )

        if not input_files:

            print(
                "\nNo valid failed simulations "
                "were found to rerun."
            )

            return

        print(
            f"\nRunning {len(input_files)} failed "
            f"simulation(s) listed in "
            f"{rerun_csv.name}."
        )

        # Rerun failed jobs regardless of
        # whether output files currently exist.
        skip_completed = False

    else:

        # ======================================================
        # Normal simulation mode
        # ======================================================

        input_files = []

        # Convert selected IDs → organ names
        selected_source_organ_names = []

        for organ_id in selected_source_organs:

            organ_name = source_organ_names.get(
                organ_id
            )

            if organ_name is None:

                print(
                    f"[WARNING] No organ name found "
                    f"for source organ ID {organ_id}."
                )

                continue

            selected_source_organ_names.append(
                organ_name
            )

        # ------------------------------------------------------
        # Search generated PHITS input files
        # ------------------------------------------------------

        for phantom in selected_phantoms:

            phantom_dir = (
                input_root /
                phantom
            )

            if not phantom_dir.exists():

                print(
                    f"\nSkipping {phantom}: "
                    f"input directory does not exist."
                )

                continue

            phantom_input_files = sorted(
                phantom_dir.rglob("*.inp")
            )

            for infile in phantom_input_files:

                filename = infile.name

                # --------------------------------------------------
                # Match source organ name in filename
                # --------------------------------------------------

                matched = False

                for organ_name in selected_source_organ_names:

                    safe_name = (
                        organ_name
                        .replace(",", "")
                        .replace(" ", "_")
                    )

                    if f"_source_{safe_name}_" in filename:

                        matched = True
                        break

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

        print("=" * 50)
        print("Existing Simulation Check")
        print("=" * 50)

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

    total_threads = (
        os.cpu_count()
    )

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