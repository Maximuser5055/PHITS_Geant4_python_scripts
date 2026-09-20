# This script asks the user for various input files and PHITS parameters

# Import necessary libraries
from pathlib import Path
import shutil
import os

import b_config.a_config as config
from b_config.b_phantom_registry import PHANTOMS, PHANTOM_GROUPS
from f_simulation_and_SAFs_further_analysis.c_check_existing_saf_database import check_existing_saf_database

def display_existing_saf_database_status(status, publishable_dir):
    """Display existing SAF database status."""

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

def get_saf_database_display_name(phantom_code):

    if phantom_code in PHANTOM_GROUPS:

        return PHANTOM_GROUPS[phantom_code].display_name

    if phantom_code in PHANTOMS:

        return PHANTOMS[phantom_code].display_name

    raise ValueError(f"Unknown phantom or phantom group: {phantom_code}")

def get_thread_count(current_threads):

    max_threads = os.cpu_count() or 1

    if current_threads > max_threads:
        print(f"\nWarning: Configured threads ({current_threads}) exceed the available logical CPU threads ({max_threads}).")
        print(f"Using {max_threads} threads instead.")
        current_threads = max_threads

    while True:

        threads_input = input(f"\nParallelization Threads [Current = {current_threads}, Max = {max_threads}]: ").strip()

        if threads_input == "":
            return current_threads

        try:
            threads = int(threads_input)

        except ValueError:
            print(
                f"\nError: Please enter an integer "
                f"between 1 and {max_threads}."
            )
            continue

        if 1 <= threads <= max_threads:
            return threads

        print(f"\nError: Number of threads must be between 1 and {max_threads}.")

def get_user_parameters():

    # ============================================================
    # Start of Internal Dosimetry UI
    # ============================================================
    print("\n")
    print("=" * 90)
    print("Internal Dosimetry Pipeline Configuration")
    print("=" * 90)

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
    
    # ============================================================
    # PHANTOM SELECTION
    # ============================================================

    phantom_options = list(PHANTOMS.keys())
    phantom_group_options = list(PHANTOM_GROUPS.keys())

    phantom_menu = []

    # Add individual phantoms
    for phantom_code in phantom_options:

        phantom_menu.append({
            "code": phantom_code,
            "type": "phantom",
        })

    # Add phantom groups
    for group_code in phantom_group_options:

        phantom_menu.append({
            "code": group_code,
            "type": "group",
        })


    print("\nWhich phantom will you do internal dosimetry on?")

    print("\nAvailable phantoms:")

    for index, item in enumerate(phantom_menu, start=1):

        code = item["code"]

        if item["type"] == "phantom":

            display_name = PHANTOMS[code].display_name

        else:

            display_name = get_saf_database_display_name(code)

        print(f"[{index}] {display_name}")


    current_phantom_display = get_saf_database_display_name(
        config.PHANTOM_INPUT_GENERATION
    )


    while True:

        choice = input(
            f"\nSelect phantom input generation "
            f"[Current = {current_phantom_display}]: "
        ).strip()

        # Press Enter to keep current selection
        if choice == "":

            phantom_input_generation = (
                config.PHANTOM_INPUT_GENERATION
            )

            break

        try:

            index = int(choice) - 1

            if 0 <= index < len(phantom_menu):

                phantom_input_generation = (
                    phantom_menu[index]["code"]
                )

                break

        except ValueError:

            pass

        print(
            f"\nError: Please enter a number from "
            f"1 to {len(phantom_menu)}.\n"
        )

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
        simulation_code
    )

    display_existing_saf_database_status(
        saf_database_status,
        publishable_dir)

    use_existing_saf_database = False
    redo_saf_calculations = True

    if saf_database_status["complete"]:

        print()
        print("=" * 90)
        print("EXISTING SAF DATABASE AVAILABLE")
        print("=" * 90)

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

        database_phantom_display = get_saf_database_display_name(phantom_input_generation)
        
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
    print("=" * 90)
    print("Specific Absorbed Fraction (SAFs) Pipeline Configuration")
    print("=" * 90)

    print(f"\n{simulation_code} is currently selected as the Monte Carlo particle transport code option.")

    # Initialize variables
    threads = config.THREADS
    phits_root = config.PHITS_INSTALLATION_DIR
    parallelization = config.PARALLELIZATION
    maxcas = config.MAXCAS
    maxbch = config.MAXBCH
    nps = config.NPS
    source_type = config.SELECTED_SOURCE_TYPE
    
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

        threads = get_thread_count(config.THREADS)

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

        threads = get_thread_count(config.THREADS)

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

    else:
        raise ValueError(f"Unsupported simulation code: {simulation_code}")

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
    