# This script automates the running of Geant4 simulations using the generated input files.

# Import necessary libraries
#from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import time

import b_config.a_config as config
from e_geant4.a_geant4_setup_and_build_executable import find_geant4make

def run_geant4():

    overall_start = time.perf_counter()

    geant4make_path = find_geant4make()

    executable = config.GEANT4_BUILD_DIR / "Internal"
    input_root = config.GEANT4_GENERATED_INPUTS_DIR

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

    while True:
        choice = input("Enter your choice (1-3): ").strip()

        if choice in {"1", "2", "3"}:
            break

        print("Invalid choice.")

    if choice == "1":
        input_files = sorted((input_root / "AM").rglob("*.in"))

    elif choice == "2":
        input_files = sorted((input_root / "AF").rglob("*.in"))

    else:
        input_files = sorted(input_root.rglob("*.in"))

    if not input_files:
        print("\nNo Geant4 input files found.")
        return

    #########################################
    # Existing simulation check
    #########################################

    print("=" * 50)
    print("Existing Simulation Check")
    print("=" * 50)
    print("[1] Skip jobs that have already been simulated")
    print("[2] Re-run completed jobs")

    while True:

        overwrite = input("Enter your choice (1-2): ").strip()

        if overwrite in {"1", "2"}:
            break

    skip_completed = overwrite == "1"

    #########################################
    # Check completion
    #########################################

    def job_completed(infile):

        outfile = infile.with_name(
            infile.stem.replace("Geant4_", "Geant4_deposit_") + ".out"
        )

        return outfile.exists()

    if skip_completed:

        original = len(input_files)

        input_files = [
            f for f in input_files
            if not job_completed(f)
        ]

        print(f"Skipping {original-len(input_files)} completed job(s).")

    if not input_files:
        print("\nNo simulations to run.")
        return

    #########################################
    # Parallelization
    #########################################

    total_threads = os.cpu_count()

    threads_per_job = config.THREADS

    max_parallel = max(1, total_threads // threads_per_job)

    print(f"\nFound {len(input_files)} jobs.")
    print(f"CPU Threads      : {total_threads}")
    print(f"Threads / Geant4 : {threads_per_job}")
    print(f"Parallel jobs    : {max_parallel}\n")

    #########################################
    # Run one Geant4 job
    #########################################

    def run_job(infile):

        relative = infile.relative_to(input_root)
        phantom = relative.parts[0]

        outfile = infile.with_name(
            infile.stem.replace("Geant4_", "Geant4_deposit_") + ".out"
        )

        source_id = (infile.parent / "source_id.txt").read_text().strip()

        cmd = (
            f'source "{geant4make_path}" && '
            f'"{executable}" '
            f'-i {source_id} '
            f'-m "{infile.resolve()}" '
            f'-o "{outfile.resolve()}"'
        )

        if phantom == "AF":
            cmd += " -f"

        print(f"Starting {infile.name}")

        process = subprocess.Popen(
            ["bash", "-c", cmd],
            cwd=config.GEANT4_BUILD_DIR
        )

        returncode = process.wait()

        if returncode == 0:
            print(f"Finished {infile.name}")
        else:
            print(f"Failed {infile.name} (return code {returncode})")

        print(f"Finished {infile.name}")

        return infile, returncode

    #########################################
    # Run jobs
    #########################################

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:

        futures = [
            executor.submit(run_job, infile)
            for infile in input_files
        ]

        for future in as_completed(futures):

            infile, returncode = future.result()

            if returncode == 0:
                print(f"✓ {infile.name}")
            else:
                print(f"✗ {infile.name} failed")

    overall_end = time.perf_counter()

    elapsed = overall_end - overall_start

    print("\n" + "="*60)
    print(f"Total elapsed wall-clock time: {elapsed:.2f} s")
    print("="*60)

    print("\nAll Geant4 jobs completed.")