# This script automates the running of Geant4 simulations using the generated input files.

# Import necessary libraries
#from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import pandas as pd

import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS
from e_geant4.a_geant4_setup_and_build_executable import find_geant4make

def run_geant4(params):

    geant4make_path = find_geant4make()

    executable = config.GEANT4_EXECUTABLE_FILE
    input_root = config.GEANT4_GENERATED_INPUTS_DIR
    current_working_dir = config.GEANT4_BUILD_DIR
    
    rerun_csv = config.RERUN_CSV
    
    #########################################
    # User Interface
    #########################################

    print("\n")
    print("=" * 50)
    print("Geant4 SAF Simulation Launcher")
    print("=" * 50)

    print("[1] Adult Male (MRCP-AM)")
    print("[2] Adult Female (MRCP-AF)")
    print("[3] Both")
    print("[4] Re-run failed simulations")

    while True:

        choice = input("Enter your choice (1-4): ").strip()
        if choice in {"1", "2", "3", "4"}:
            break
        print("Invalid choice. Please enter 1, 2, 3, or 4.")

    # ==========================================================
    # Determine phantom(s)
    # ==========================================================

    if choice == "1":

        selected_phantoms = ["AM"]

    elif choice == "2":

        selected_phantoms = ["AF"]

    elif choice == "3":

        selected_phantoms = ["AM", "AF"]

    else:

        # Rerun mode will determine the phantom from
        # the input-file paths.
        selected_phantoms = ["AM", "AF"]


    #########################################
    # Source Organ Selection
    #########################################

    if choice != "4":

        source_csv = config.SOURCE_CSV

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
        # Source organ column
        # ------------------------------------------------------

        source_organ_column = "source_organ_ID"

        if source_organ_column not in source_df.columns:

            raise ValueError(
                f"Column '{source_organ_column}' was not found "
                f"in {source_csv}"
            )

        # ------------------------------------------------------
        # Get source-organ IDs from SOURCE_CSV
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
        # Build source-organ name lookup
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
        # User selects source organs
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
                    or index >= len(available_source_organs)
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

    if choice == "4":

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
                relative_path
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

        # Rerun jobs regardless of existing output
        skip_completed = False

    else:

        # ======================================================
        # Normal simulation mode
        # ======================================================

        input_files = []

        for phantom in selected_phantoms:

            phantom_dir = input_root / phantom

            # ------------------------------------------------------
            # Skip phantom if its generated-input directory
            # does not exist
            # ------------------------------------------------------

            if not phantom_dir.exists():

                print(
                    f"\nSkipping {phantom}: "
                    f"input directory does not exist."
                )

                continue

            phantom_input_files = sorted(
                phantom_dir.rglob("*.in")
            )

            for infile in phantom_input_files:

                # --------------------------------------------------
                # Read source ID from source_id.txt
                # --------------------------------------------------

                source_id_file = (
                    infile.parent / "source_id.txt"
                )

                if not source_id_file.exists():

                    print(
                        f"[WARNING] Missing source_id.txt for "
                        f"{infile.name}"
                    )

                    continue

                source_id = (
                    source_id_file
                    .read_text()
                    .strip()
                )

                # --------------------------------------------------
                # Match against selected source-organ IDs
                # --------------------------------------------------

                if source_id in selected_source_organs:

                    input_files.append(infile)

        # Remove duplicates and sort
        input_files = sorted(set(input_files))

        # ----------------------------------------------------------
        # Check whether anything was found
        # ----------------------------------------------------------

        if not input_files:

            print(
                "\nNo Geant4 input files were found "
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
            "[2] Re-run completed jobs (overwrite outputs)"
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
    # Check completion
    # ==========================================================

    def job_completed(infile):

        outfile = infile.with_name(
            infile.stem.replace(
                "Geant4_",
                "Geant4_deposit_"
            ) + ".out"
        )

        return outfile.exists()

    # ==========================================================
    # Filter completed jobs
    # ==========================================================

    if skip_completed:

        original = len(input_files)

        input_files = [
            infile
            for infile in input_files
            if not job_completed(infile)
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
    # Parallelization
    # ==========================================================

    total_threads = os.cpu_count()

    threads_per_job = params["threads"]

    max_parallel = max(
        1,
        total_threads // threads_per_job
    )

    print(
        f"\nFound {len(input_files)} "
        f"input files."
    )

    print(
        f"Total CPU threads  : "
        f"{total_threads}"
    )

    print(
        f"Threads per Geant4 : "
        f"{threads_per_job}"
    )

    print(
        f"Parallel Geant4 jobs: "
        f"{max_parallel}"
    )

    print(
        f"\nRunning up to "
        f"{max_parallel} jobs simultaneously.\n"
    )

    # ==========================================================
    # Run one Geant4 job
    # ==========================================================

    def run_job(infile):

        relative = infile.relative_to(
            input_root
        )

        phantom = relative.parts[0]

        outfile = infile.with_name(
            infile.stem.replace(
                "Geant4_",
                "Geant4_deposit_"
            ) + ".out"
        )

        source_id_file = (
            infile.parent /
            "source_id.txt"
        )

        if not source_id_file.exists():

            print(
                f"Missing source_id.txt for "
                f"{infile.name}"
            )

            return infile, 1

        source_id = (
            source_id_file
            .read_text()
            .strip()
        )

        cmd = (
            f'source "{geant4make_path}" && '
            f'"{executable}" '
            f'-i {source_id} '
            f'-m "{infile.resolve()}" '
            f'-o "{outfile.resolve()}"'
        )

        # ------------------------------------------------------
        # AF phantom flag
        # ------------------------------------------------------

        if phantom == "AF":

            cmd += " -f"

        print(
            f"Starting {infile.name}"
        )

        process = subprocess.Popen(
            ["bash", "-c", cmd],
            cwd=current_working_dir
        )

        returncode = process.wait()

        if returncode == 0:

            print(
                f"Finished {infile.name}"
            )

        else:

            print(
                f"Failed {infile.name} "
                f"(return code {returncode})"
            )

        return infile, returncode

    # ==========================================================
    # Run jobs in parallel
    # ==========================================================

    with ThreadPoolExecutor(
        max_workers=max_parallel
    ) as executor:

        futures = [
            executor.submit(
                run_job,
                infile
            )
            for infile in input_files
        ]

        for future in as_completed(
            futures
        ):

            infile, returncode = (
                future.result()
            )

            if returncode == 0:

                print(
                    f"✓ {infile.name} "
                    f"completed successfully"
                )

            else:

                print(
                    f"✗ {infile.name} "
                    f"failed "
                    f"(return code "
                    f"{returncode})"
                )

    print(
        "\nAll Geant4 jobs finished."
    )