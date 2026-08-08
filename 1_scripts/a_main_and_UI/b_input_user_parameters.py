# This script asks the user for various input files and PHITS parameters

# Import necessary libraries
import b_config.a_config as config
from pathlib import Path

def get_user_parameters():

    # User Interface
    print("\n")
    print("=" * 50)
    print("Internal Dosimetry Pipeline Configuration")
    print("=" * 50)

    print(f"\nOperating System : {config.SYSTEM}")
    print(f"Running in WSL   : {config.IS_WSL}")

    print("\nMonte Carlo particle transport code options for the working OS:")

    if config.IS_WINDOWS:
        print("[1] PHITS")
        print("[2] Geant4 (Not available on Windows)")

        while True:
            choice = input("Select transport code (1): ").strip()

            if choice in ("", "1"):
                simulation_code = "PHITS"
                break

            elif choice == "2":
                print("\nGeant4 is only supported on Linux or WSL2.")
                print("Please use PHITS instead.\n")

            else:
                print("Invalid choice. Please enter 1.")

    elif config.IS_LINUX or config.IS_WSL:
        print("[1] PHITS")
        print("[2] Geant4")

        while True:
            choice = input("Select transport code (1-2): ").strip()

            if choice == "1":
                simulation_code = "PHITS"
                break

            elif choice == "2":
                simulation_code = "GEANT4"
                break

            print("Invalid choice. Please enter 1 or 2.")

    else:
        raise RuntimeError("Unsupported operating system.")

    if simulation_code == "PHITS":

        phits_root = input(f"\nPHITS root [Current = {config.PHITS_ROOT}]: ").strip()
        phits_root = Path(phits_root) if phits_root else config.PHITS_ROOT

        parallelization = (
            input(f"Parallelization (OMP/MPI) [Current = {config.PARALLELIZATION}]: ")
            .strip()
            .upper()
        )
        if parallelization == "":
            parallelization = config.PARALLELIZATION

        threads = input(f"Parallelization Threads [Current = {config.THREADS}]: ").strip()
        threads = int(threads) if threads else config.THREADS

        maxcas = input(f"PHITS maxcas (no. of particle histories per batch) [Current = {config.MAXCAS}]: ").strip()
        maxcas = int(maxcas) if maxcas else config.MAXCAS

        maxbch = input(f"PHITS maxbch (no. of batches) [Current = {config.MAXBCH}]: ").strip()
        maxbch = int(maxbch) if maxbch else config.MAXBCH

    elif simulation_code == "GEANT4":

        threads = input(f"\nParallelization Threads [Current = {config.THREADS}]: ").strip()
        threads = int(threads) if threads else config.THREADS

        nps = input(f"GEANT4 nps (no. of particle histories) [Current = {config.NPS}]: ").strip()
        nps = int(nps) if nps else config.NPS


    while True:
        source_dir = input(
            f"Directory containing source_organs.csv "
            f"[Current = {config.SOURCE_CSV.parent}]: "
        ).strip()

        source_csv = (
            Path(source_dir) / "source_organs.csv"
            if source_dir
            else config.SOURCE_CSV
        )

        if source_csv.is_file():
            break

        print(f"\nError: '{source_csv}' was not found. Please try again.\n")

    while True:
        phantom_input_generation = input(
            f"Phantom Input Generation (AM/AF/Both) [Current = {config.PHANTOM_INPUT_GENERATION}]: "
        ).strip().upper()

        if phantom_input_generation == "":
            phantom_input_generation = config.PHANTOM_INPUT_GENERATION

        if phantom_input_generation in ("AM", "AF", "Both"):
            break

        print("\nError: Please enter AM, AF, or Both.\n")

    uncertainty_limit = input(
        f"Maximum allowed statistical uncertainty (%) "
        f"[Current = {config.UNCERTAINTY_LIMIT}]: "
    ).strip()

    uncertainty_limit = (
        float(uncertainty_limit)
        if uncertainty_limit
        else config.UNCERTAINTY_LIMIT
    )

    config.update_config("SIMULATION_CODE", simulation_code)
    config.update_config("UNCERTAINTY_LIMIT", uncertainty_limit)
    config.update_config("THREADS", threads)
    config.update_config("SOURCE_CSV", source_csv)
    config.update_config("PHANTOM_INPUT_GENERATION", phantom_input_generation)
    
    if simulation_code == "PHITS":
        config.update_config("PHITS_ROOT", phits_root)
        config.update_config("PARALLELIZATION", parallelization)
        config.update_config("MAXCAS", maxcas)
        config.update_config("MAXBCH", maxbch)

    elif simulation_code == "GEANT4":
        config.update_config("NPS", nps)

    ############################
    # Create required directories
    ############################

    DIRECTORIES = [
        config.GENERATED_INPUTS_DIR,
        config.RESULTS_DIR,
    ]

    for directory in DIRECTORIES:
        directory.mkdir(parents=True, exist_ok=True)

    params = {
    "simulation_code": simulation_code,
    "uncertainty_limit": uncertainty_limit,
    "threads": threads,
    "source_csv": source_csv,
    "phantom": phantom_input_generation,
    }

    if simulation_code == "PHITS":
        params.update({
        "phits_root": phits_root,
        "parallelization": parallelization,
        "maxcas": maxcas,
        "maxbch": maxbch,
    })

    elif simulation_code == "GEANT4":
        params.update({
        "nps": nps,
    })

    return params
    