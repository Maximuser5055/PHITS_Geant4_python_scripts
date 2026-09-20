# This script extracts the dose and SAFs from the Geant4 CSV output files and merging them while deleting the individual CSV files.

# Import necessary libraries
import re
import pandas as pd

import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom, get_phantom_group
from c_database.b_organ_database import ORGANS

def geant4_calculate_dose_and_SAFs(params):

    # -------------------------------------------------------------
    # Constants
    # -------------------------------------------------------------

    MeV_to_J = config.MEV_TO_J

    # -------------------------------------------------------------
    # Directories and Configs
    # -------------------------------------------------------------

    input_root = config.GEANT4_GENERATED_INPUTS_DIR
    output_root = config.RESULTS_GEANT4_DIR
    source_type_map = config.GEANT4_SOURCE_TYPE_MAP
    
    # -------------------------------------------------------------
    # File names
    # -------------------------------------------------------------

    geant4_output_file = "geant4_dose_and_SAFs.csv"

    filename_pattern = re.compile(
        r"geant4_deposit_(.+?)_source_"
        r"(.+?)_"
        r"(.+?)_energy_"
        r"([0-9Ee.+-]+)\.csv",
        re.IGNORECASE
    )

    # -------------------------------------------------------------------------
    # Determine selected phantom family
    # -------------------------------------------------------------------------

    phantom_selection = params["phantom"]
    selected_phantoms = get_phantom_group(phantom_selection)

    output_file = output_root / geant4_output_file

    # -------------------------------------------------------------------------
    # Find all deposit tally files
    # -------------------------------------------------------------------------

    deposit_files = sorted(file
        for file in input_root.rglob("geant4_deposit_*.csv")
        if (
            not file.stem.lower().endswith("_photon_fluence")
            and any(
                file.name.lower().startswith(
                    f"geant4_deposit_{phantom.lower()}_"
                )
                for phantom in selected_phantoms
            )
        )
    )

    if not deposit_files:

        print("\n[WARNING] No Geant4 deposit CSV files found.")

        return

    print(f"\nFound {len(deposit_files)} Geant4 deposit CSV file(s).\n")

    # -------------------------------------------------------------
    # Organ name -> organ ID lookup
    # -------------------------------------------------------------

    ORGAN_NAME_TO_ID = {
        phantom: {
            organ["name"]: organ["organ_id"]
            for organ in ORGANS[phantom].values()
        }
        for phantom in ORGANS
    }

    # -------------------------------------------------------------
    # Read one CSV
    # -------------------------------------------------------------

    processed_files = []

    all_tallies = []

    for deposit_file in deposit_files:

        print(f"Processing {deposit_file.name}")

        match = filename_pattern.match(deposit_file.name)

        if not match:

            raise ValueError(
                f"Cannot parse filename:\n{deposit_file.name}"
            )

        phantom_code = match.group(1).upper()
        phantom_spec = get_phantom(phantom_code)

        phantom = phantom_spec.display_name

        source_organ = match.group(2)
        source_type = source_type_map.get(match.group(3).lower(), match.group(3).lower())
        source_energy = float(match.group(4))

        number_of_particles = params["nps"]

        source_organ_id = ORGAN_NAME_TO_ID[
            phantom_code
        ][source_organ]

        organ_database = ORGANS[phantom_code]

        df = pd.read_csv(deposit_file)

        df.rename(
            columns={
                "Organ ID": "Target Organ ID",
                "Organ Mass (g)": "Target Organ Mass (g)"
            },
            inplace=True
        )

        # Keep only organs in the database

        df = df[
            df["Target Organ ID"].isin(
                organ_database.keys()
            )
        ].copy()

        # Organ names

        df["Target Organ Name"] = df[
            "Target Organ ID"
        ].map(
            lambda x: organ_database[int(x)]["name"]
            if pd.notna(x) else pd.NA
        )

        # Use masses from database

        df["Target Organ Mass (g)"] = df[
            "Target Organ ID"
        ].map(
            lambda x: organ_database[int(x)]["mass"]
            if pd.notna(x) else pd.NA
        )

        # SAF

        df["SAF (kg^-1)"] = (

            df["Dose (Gy/source)"] /

            (source_energy * MeV_to_J)

        )

        # Statistical uncertainty

        df["Statistical Uncertainty (%)"] = (

            df["Relative Error"] * 100

        )

        # Insert metadata columns

        df.insert(
            0,
            "Phantom",
            phantom
        )

        df.insert(
            1,
            "Source Organ ID",
            source_organ_id
        )

        df.insert(
            2,
            "Source Organ Name",
            source_organ
        )

        df.insert(
            3,
            "Source Type",
            source_type
        )

        df.insert(
            4,
            "Source Energy (MeV)",
            source_energy
        )

        df.insert(
            5,
            "Number of Particles",
            number_of_particles
        )
        
        # Reorder

        df = df[
            [

                "Phantom",

                "Source Organ ID",
                "Source Organ Name",

                "Source Type",
                "Source Energy (MeV)",

                "Number of Particles",

                "Target Organ ID",
                "Target Organ Name",
                "Target Organ Mass (g)",

                "Dose (Gy/source)",

                "SAF (kg^-1)",

                "Relative Error",

                "Statistical Uncertainty (%)"

            ]
        ]

        all_tallies.append(df)

        processed_files.append(deposit_file)

    # -------------------------------------------------------------
    # Combine
    # -------------------------------------------------------------

    combined_df = pd.concat(
        all_tallies,
        ignore_index=True
    )

    sort_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Type",
        "Source Energy (MeV)",
        "Target Organ ID"
    ]

    combined_df.sort_values(by=sort_columns, inplace=True)
    combined_df.reset_index(drop=True, inplace=True)

    # -------------------------------------------------------------
    # Save
    # -------------------------------------------------------------

    combined_df.to_csv(output_file, index=False)

    # -------------------------------------------------------------
    # Delete individual CSVs
    # -------------------------------------------------------------

    deleted = 0

    for file in processed_files:

        try:

            file.unlink()
            deleted += 1

        except Exception as e:

            print(
                f"[WARNING] Could not delete "
                f"{file.name}: {e}"
            )

    print("\nFinished extracting all Geant4 deposit CSV files.")

    print(f"Wrote:\n  {output_file}\n")

    print(f"Deleted {deleted} individual deposit CSV file(s).")