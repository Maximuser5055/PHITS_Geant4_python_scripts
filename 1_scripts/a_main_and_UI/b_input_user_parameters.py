# This script asks the user for various input files and PHITS parameters

# Import necessary libraries
from pathlib import Path
import shutil

import b_config.a_config as config
from f_simulation_and_SAFs_further_analysis.c_check_uncertainty import check_existing_saf_database

def display_existing_saf_database_status(status, uncertainty_limit, publishable_dir,):
    """Display existing SAF database status and uncertainty."""

    print()
    print("=" * 90)
    print("EXISTING SAF DATABASE CHECK")
    print("=" * 90)

    print(
        f"\nSAF database directory:\n"
        f"    {publishable_dir}"
    )

    if not status["exists"]:
        print(
            "\nNo existing SAF database was found for "
            "the selected phantom."
        )

        return

    if not status["complete"]:
        print("\nAn incomplete SAF database was found.")

        print("\nExisting files:")
        for file in status["existing_files"]:
            print(f"    {file.name}")

        print("\nMissing files:")
        for file in status["missing_files"]:
            print(f"    {file.name}")

        return

    print("\nComplete SAF database found.")

    print("\nFiles:")
    for file in status["existing_files"]:
        print(f"    {file.name}")

    print("\nMaximum statistical uncertainty:")

    for file_name, maximum in status["uncertainty_by_file"].items():

        file_status = (
            "PASS"
            if maximum < uncertainty_limit
            else "FAIL"
        )

        print(
            f"    {file_name:<35}"
            f"{maximum:8.2f} %   [{file_status}]"
        )

    overall_status = (
        "PASS"
        if status["uncertainty_pass"]
        else "FAIL"
    )

    print()
    print(f"Uncertainty limit : {uncertainty_limit:.2f} %")
    print(
        f"Overall maximum   : "
        f"{status['max_uncertainty']:.2f} %"
    )
    print(f"Overall status    : {overall_status}")

def get_user_parameters():

    # ============================================================
    # Start of Internal Dosimetry UI
    # ============================================================
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
    
    while True:
            print("\nWhich phantom will you do internal dosimetry on?:")
            print("ICRP 145 Mesh-type Reference Computational Phantoms")
            print("[1] MRCP AF (Adult Female)")
            print("[2] MRCP AM (Adult Male)")
            print("\nFilipino-based Mesh-type Computational Phantoms")
            print("[3] MFCP AF")
            print("[4] MFCP AM")
    
            print("\nPhantom sets")
            print("[5] Both MRCP phantoms (AF + AM)")
            print("[6] Both MFCP phantoms (AF + AM)")
    
            current_phantom_display = config.PHANTOM_DISPLAY_NAMES.get(
                config.PHANTOM_INPUT_GENERATION,
                config.PHANTOM_INPUT_GENERATION)
            
            choice = input(
                f"Select phantom input generation "
                f"[Current = {current_phantom_display}]: "
            ).strip()
    
    
            if choice == "":
                phantom_input_generation = (config.PHANTOM_INPUT_GENERATION)
                break
    
            if choice == "1":
                phantom_input_generation = "MRCP_AF"
                break
    
            elif choice == "2":
                phantom_input_generation = "MRCP_AM"
                break
    
            elif choice == "3":
                phantom_input_generation = "MFCP_AF"
                break
    
            elif choice == "4":
                phantom_input_generation = "MFCP_AM"
                break
    
            elif choice == "5":
                phantom_input_generation = "MRCP_AF_AM"
                break
    
            elif choice == "6":
                phantom_input_generation = "MFCP_AF_AM"
                break
    
            print("\nError: Please enter 1, 2, 3, 4, 5, or 6.\n")

    # ============================================================
    # EXISTING SAF DATABASE
    # ============================================================

    if simulation_code == "PHITS":

        publishable_dir = config.RESULTS_PHITS_PUBLISHABLE_SAF_DATABASE_DIR

    elif simulation_code == "GEANT4":

        publishable_dir = config.RESULTS_GEANT4_PUBLISHABLE_SAF_DATABASE_DIR

    else:

        raise ValueError(f"Unsupported simulation code: {simulation_code}")

    saf_database_status = check_existing_saf_database(
        phantom_input_generation,
        simulation_code,
        config.UNCERTAINTY_LIMIT
    )

    display_existing_saf_database_status(
        saf_database_status,
        config.UNCERTAINTY_LIMIT,
        publishable_dir,)

    use_existing_saf_database = False
    redo_saf_calculations = True

    if saf_database_status["complete"]:

        print()
        print("=" * 90)
        print("EXISTING SAF DATABASE AVAILABLE")
        print("=" * 90)

        if saf_database_status["uncertainty_pass"]:
            print(
                "\nThe existing SAF database satisfies "
                "the selected uncertainty limit."
            )
        else:
            print(
                "\nWARNING: The existing SAF database does NOT "
                f"satisfy the selected uncertainty limit of {config.UNCERTAINTY_LIMIT} %."
            )

        print("\nWhat would you like to do?")
        print(
            "[1] Use the existing SAF database and "
            "proceed to S-value calculation"
        )
        print("[2] Redo the SAF calculations")

        while True:

            choice = input("\nSelect option (1-2): ").strip()

            if choice == "1":
                use_existing_saf_database = True
                redo_saf_calculations = False
                print("\nUsing the existing SAF database.")
                
                if use_existing_saf_database:

                    print()
                    print("=" * 90)
                    print("SKIPPING SAF CALCULATION")
                    print("=" * 90)

                    print(
                        "\nThe existing publishable SAF database will be used."
                    )

                    print(
                        f"Simulation code : {simulation_code}"
                    )

                    print(
                        f"SAF database    : {publishable_dir}"
                    )

                    print(
                        "\nProceeding directly to S-value calculation..."
                    )

                    return {
                        "use_existing_saf_database": True,
                        "redo_saf_calculations": False,
                        "saf_database_status": saf_database_status,
                        "simulation_code": simulation_code,
                        "uncertainty_limit": config.UNCERTAINTY_LIMIT,
                        "phantom": phantom_input_generation,
                        "saf_database_dir": publishable_dir,
                    }

            elif choice == "2":
                use_existing_saf_database = False
                redo_saf_calculations = True
                print("\nThe SAF calculations will be redone.")
                break

            print("Invalid choice. Please enter 1 or 2.")

    else:

        database_phantom_display = config.SAF_DATABASE_PHANTOM_NAMES.get(
            phantom_input_generation,
            phantom_input_generation,
        )
        
        print()
        print("=" * 90)
        print(f"NO COMPLETE SAF DATABASE AVAILABLE FOR {simulation_code} {database_phantom_display}")
        print("=" * 90)

        print("\nS-value calculation cannot proceed yet.")
        print("The SAF calculations must be performed first.")

        use_existing_saf_database = False
        redo_saf_calculations = True

    # ============================================================
    # Start of SAF calculation pipeline
    # ============================================================

    print("\n")
    print("=" * 50)
    print("Specific Absorbed Fraction (SAFs) Pipeline Configuration")
    print("=" * 50)

    print(f"\n{simulation_code} is currently selected as the Monte Carlo particle transport code option.")

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

    print(f"\nMaximum allowed statistical uncertainty: {config.UNCERTAINTY_LIMIT} %")

    config.update_config("SIMULATION_CODE", simulation_code)
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

    if redo_saf_calculations and generated_inputs_dir.exists():

        while True:

            
            choice = input(
                f"\nWARNING: Existing generated inputs were found in:\n"
                f"    {generated_inputs_dir}\n\n"
                "Note: It is best to delete the entire inputs folder."
                f"\nDelete this entire folder and create a fresh start? "
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
        config.RESULTS_GEANT4_PUBLISHABLE_SAF_DATABASE_DIR,
        config.RESULTS_PHITS_PUBLISHABLE_SAF_DATABASE_DIR,
        config.RESULTS_S_VALUES_DIR,
    ]

    for directory in DIRECTORIES:
        directory.mkdir(parents=True, exist_ok=True)

    params = {
    # Check SAF database workflow
    "use_existing_saf_database": use_existing_saf_database,
    "redo_saf_calculations": redo_saf_calculations,
    "saf_database_status": saf_database_status,

    # SAF pipeline configuration
    "simulation_code": simulation_code,
    "uncertainty_limit": config.UNCERTAINTY_LIMIT,
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
    