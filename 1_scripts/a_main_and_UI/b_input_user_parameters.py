# This script asks the user for various input files and PHITS parameters

# Import necessary libraries
import b_config.a_config as config
from pathlib import Path
import shutil

def get_user_parameters():

    # User Interface
    print("\n")
    print("=" * 50)
    print("Internal Dosimetry Pipeline Configuration")
    print("=" * 50)

    print(f"\nOperating System : {config.SYSTEM}")

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
                print("\nGeant4 is only supported on Linux")
                print("Please use PHITS instead.\n")

            else:
                print("Invalid choice. Please enter 1.")

    elif config.IS_LINUX:
        print("[1] PHITS")
        print("[2] GEANT4")

        while True:
            choice = input(f"Select transport code (1-2) [Current = {config.SIMULATION_CODE}]: ").strip()

            if choice == "":
                simulation_code = config.SIMULATION_CODE
                break

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

        phits_root = input(f"\nPHITS installation directory [Current = {config.PHITS_INSTALLATION_DIR}]: ").strip()
        phits_root = Path(phits_root) if phits_root else config.PHITS_INSTALLATION_DIR

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

        # ========================================================
        # SOURCE TYPE
        # ========================================================

        print("\nSource type:")
        print("[1] photon / gamma")
        print("[2] electron / e-")

        current_source_type = (
            "photon / gamma"
            if config.SELECTED_SOURCE_TYPE in {"gamma", "photon"}
            else "electron / e-"
        )

        while True:

            choice = input(
                f"Select source type (1-2) "
                f"[Current = {current_source_type}]: "
            ).strip()

            if choice == "":
                if config.SELECTED_SOURCE_TYPE in {"gamma", "photon"}:
                    source_type = "photon"
                else:
                    source_type = "electron"
                break

            if choice == "1":
                source_type = "photon"
                break

            elif choice == "2":
                source_type = "electron"
                break

            print("Invalid choice. Please enter 1 or 2.")

    elif simulation_code == "GEANT4":

        threads = input(f"\nParallelization Threads [Current = {config.THREADS}]: ").strip()
        threads = int(threads) if threads else config.THREADS

        nps = input(f"GEANT4 nps (no. of particle histories) [Current = {config.NPS}]: ").strip()
        nps = int(nps) if nps else config.NPS

        # ========================================================
        # SOURCE TYPE
        # ========================================================

        print("\nSource type:")
        print("[1] photon / gamma")
        print("[2] electron / e-")

        current_source_type = (
            "photon / gamma"
            if config.SELECTED_SOURCE_TYPE in {"gamma", "photon"}
            else "electron / e-"
        )

        while True:

            choice = input(
                f"Select source type (1-2) "
                f"[Current = {current_source_type}]: "
            ).strip()

            if choice == "":
                if config.SELECTED_SOURCE_TYPE in {"gamma", "photon"}:
                    source_type = "gamma"
                else:
                    source_type = "e-"
                break

            if choice == "1":
                source_type = "gamma"
                break

            elif choice == "2":
                source_type = "e-"
                break

            print("Invalid choice. Please enter 1 or 2.")

    # ============================================================
    # SOURCE ENERGIES
    # ============================================================

    while True:

        energy_input = input(
            "Source energies in MeV "
            f"[Current = {config.SOURCE_ENERGIES}]: "
        ).strip()

        if energy_input == "":
            source_energies = config.SOURCE_ENERGIES.copy()
            break

        try:

            source_energies = [
                float(x.strip())
                for x in energy_input.split(",")
            ]

            if not source_energies:
                raise ValueError

            if any(energy <= 0 for energy in source_energies):
                raise ValueError

            # Remove duplicate energies while preserving order
            source_energies = list(
                dict.fromkeys(source_energies)
            )

            break

        except ValueError:

            print(
                "\nError: Please enter positive "
                "comma-delimited energies."
            )

    while True:
        source_dir = input(
            f"\nDirectory containing source_organs.csv "
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
        print("\nPhantom Input Generation:")
        print("[1] AM")
        print("[2] AF")
        print("[3] Both")

        choice = input(
            f"Select phantom input generation "
            f"[Current = {config.PHANTOM_INPUT_GENERATION}]: "
        ).strip()


        if choice == "":
            phantom_input_generation = (config.PHANTOM_INPUT_GENERATION)
            break

        if choice == "1":
            phantom_input_generation = "AM"
            break

        elif choice == "2":
            phantom_input_generation = "AF"
            break

        elif choice == "3":
            phantom_input_generation = "Both"
            break

        print("\nError: Please enter 1, 2, or 3.\n")

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
    config.update_config("SOURCE_ENERGIES", source_energies)
    config.update_config("SELECTED_SOURCE_TYPE", source_type)

    if simulation_code == "PHITS":
        config.update_config("PHITS_INSTALLATION_DIR", phits_root)
        config.update_config("PARALLELIZATION", parallelization)
        config.update_config("MAXCAS", maxcas)
        config.update_config("MAXBCH", maxbch)

    elif simulation_code == "GEANT4":
        config.update_config("NPS", nps)

    ############################
    # Fresh start
    ############################

    if simulation_code == "PHITS":
        generated_inputs_dir = config.GENERATED_INPUTS_DIR

    elif simulation_code == "GEANT4":
        generated_inputs_dir = config.GEANT4_GENERATED_INPUTS_DIR

    if generated_inputs_dir.exists():

        while True:

            choice = input(
                f"\nWARNING: Existing generated inputs were found in:\n"
                f"    {generated_inputs_dir}\n\n"
                f"Delete this entire folder and create a fresh start? "
                f"(Y/N): "
            ).strip().upper()

            if choice in ("Y", "YES"):
                shutil.rmtree(generated_inputs_dir)
                print("\nGenerated inputs folder deleted.")

                break

            elif choice in ("N", "NO"):
                print("\nKeeping existing generated inputs.")

                break

            else:
                print("\nInvalid choice. Please enter Y or N.")

    ############################
    # Create required directories
    ############################

    DIRECTORIES = [
        config.GENERATED_INPUTS_DIR,
        config.RESULTS_DIR,
        config.RESULTS_PHITS_DIR,
        config.RESULTS_GEANT4_DIR,
        config.RESULTS_SAF_DATABASE_DIR,
        config.RESULTS_S_VALUES_DIR,
    ]

    for directory in DIRECTORIES:
        directory.mkdir(parents=True, exist_ok=True)

    params = {
    "simulation_code": simulation_code,
    "uncertainty_limit": uncertainty_limit,
    "threads": threads,
    "source_csv": source_csv,
    "phantom": phantom_input_generation,
    "source_type": source_type,
    "source_energies": source_energies,
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
    