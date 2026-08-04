# This script runs PHITS simulations for all generated input files in parallel, using a specified number of threads.
# It runs programs simultaneously depending on the number of available CPU cores and the threads specified in c_generating_inputs.py.
# After one simulation is complete, it runs the next one immediately until all input files have been processed.

from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import pandas as pd
import b_config.a_config as config

def run_phits():
    # Batch file to run PHITS
    phits = config.PHITS_BAT

    # Root directory containing all generated input files and the csv for rerunning simulations with uncertainty greater than 5%
    input_root = config.GENERATED_INPUTS_DIR 
    rerun_csv = config.RERUN_CSV

    # Pattern in the input files to determine the number of threads
    patterns = {
        "parallelization_pattern": re.compile(r"\$(OMP|MPI)\s*=\s*(\d+)", re.IGNORECASE)
    }

    # User interface for selecting which phantom(s) to run simulations for
    print("\n")
    print("=" * 50)
    print("PHITS Specific Absorbed Fraction (SAF) Simulation Launcher")
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

    # Determine which input files to run based on user choice
    if choice == "1":

        input_files = sorted((input_root / "AM").rglob("*.inp"))

    elif choice == "2":

        input_files = sorted((input_root / "AF").rglob("*.inp"))

    elif choice == "3":

        input_files = sorted(input_root.rglob("*.inp"))

    else:

        if not rerun_csv.exists():
            raise FileNotFoundError(
                "5_rerun_required.csv not found.\n"
                "Run h_check_uncertainty.py first."
            )

        rerun = pd.read_csv(rerun_csv)

        if rerun.empty:
            print("\nThere are no failed simulations to rerun.")
            raise SystemExit

        if "Input File" not in rerun.columns:
            raise ValueError(
                "The rerun file does not contain an 'Input File' column.\n"
                "Run the latest version of h_check_uncertainty.py."
            )

        input_files = []

        for relative_path in rerun["Input File"]:

            input_file = input_root / Path(relative_path)

            if input_file.exists():
                input_files.append(input_file)
            else:
                print(f"Missing input file: {relative_path}")

        input_files = sorted(set(input_files))

    # Check that there are input files to simulate2

    if not input_files:
        if choice == "4":
            print("\nNo failed simulations found in 5_rerun_required.csv.")
        else:
            print("\nNo PHITS input files were found.")

        exit()
        
    if choice != "4":

        print("=" * 50)
        print("Existing Simulation Check")
        print("=" * 50)
        print("[1] Skip jobs that have already been simulated")
        print("[2] Re-run completed jobs (overwrite outputs)")

        while True:
            overwrite_choice = input("Enter your choice (1-2): ").strip()
            if overwrite_choice in {"1", "2"}:
                break
            print("Invalid choice. Please enter 1 or 2.")

        SKIP_COMPLETED = overwrite_choice == "1"

    else:
        print(
            f"\nRunning {len(input_files)} failed simulation(s) "
            f"listed in {rerun_csv.name}."
            )
        SKIP_COMPLETED = False

    # Filter input files based on whether they have already been completed
    def job_completed(inp_file):
        job_dir = inp_file.parent

        batch_exists = (job_dir / "batch.out").exists()

        deposit_exists = any(
            f.name.startswith("deposit_") and f.suffix == ".out"
            for f in job_dir.iterdir()
        )

        phits_exists = any(
            f.name.startswith("phits_") and f.suffix == ".out"
            for f in job_dir.iterdir()
        )

        return batch_exists and deposit_exists and phits_exists

    # Find the number of threads for an input file
    def get_threads(input_file):
        """Read the number of OMP/MPI threads from a PHITS input file."""

        with open(input_file, "r") as f:
            for line in f:
                match = patterns["parallelization_pattern"].search(line)
                if match:
                    return int(match.group(2))

        raise ValueError(
            f"No $OMP or $MPI definition found in {input_file}"
        )

    # Run one PHITS simulation
    def run_phits(input_file):

        print(f"Starting {input_file.name}")

        process = subprocess.Popen(
            [phits, input_file.name],
            cwd=input_file.parent
        )

        returncode = process.wait()

        print(f"Finished {input_file.name}")

        return input_file, returncode

    # Filter input files based on completion status
    if SKIP_COMPLETED:
        original = len(input_files)

        input_files = [
            inp for inp in input_files
            if not job_completed(inp)
        ]

        skipped = original - len(input_files)

        print(f"\nSkipping {skipped} completed job(s).")
        print(f"Remaining jobs: {len(input_files)}")

    if len(input_files) == 0:
        print("\nNo simulations to run. All jobs have already been completed.")
        exit()

    # Maximum number of PHITS jobs to run simultaneously
    TOTAL_THREADS = os.cpu_count()

    threads = get_threads(input_files[0])

    MAX_PARALLEL = max(1, TOTAL_THREADS // threads)

    print(f"Found {len(input_files)} input files.")
    print(f"Total CPU threads : {TOTAL_THREADS}")
    print(f"Threads per PHITS : {threads}")
    print(f"Parallel PHITS jobs: {MAX_PARALLEL}\n")
    print(f"Running up to {MAX_PARALLEL} jobs simultaneously.\n")

    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as executor:

        futures = [executor.submit(run_phits, inp) for inp in input_files]

        for future in as_completed(futures):

            input_file, returncode = future.result()

            if returncode == 0:
                print(f"✓ {input_file.name} completed successfully")
            else:
                print(f"✗ {input_file.name} failed (return code {returncode})")

    print("\nAll jobs finished.")