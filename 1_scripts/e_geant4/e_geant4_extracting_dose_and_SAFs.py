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

    geant4_output_file = "c_geant4_dose_and_SAFs.csv"

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
    # Load source-region metadata from the same CSV used by the
    # input generators. Keep IDs as strings because one source region
    # may contain multiple IDs (e.g. 9700_9900).
    # -------------------------------------------------------------

    source_csv = config.SOURCE_CSV
    if not source_csv.is_file():
        raise FileNotFoundError(f"Source-organ CSV not found: {source_csv}")

    source_df = pd.read_csv(source_csv, dtype=str, encoding="utf-8-sig")
    source_df.columns = [str(column).strip() for column in source_df.columns]

    required_columns = {"Acronym", "ID number(s)"}
    missing_columns = required_columns - set(source_df.columns)
    if missing_columns:
        raise ValueError(
            f"Source-organ CSV is missing column(s): {sorted(missing_columns)}. "
            f"Found columns: {list(source_df.columns)}. File: {source_csv}"
        )

    source_df["Acronym"] = source_df["Acronym"].fillna("").str.strip()
    source_df["ID number(s)"] = source_df["ID number(s)"].fillna("").str.strip()
    source_df = source_df[
        (source_df["Acronym"] != "") & (source_df["ID number(s)"] != "")
    ]

    source_name_to_id = dict(
        zip(source_df["Acronym"], source_df["ID number(s)"])
    )
    if not source_name_to_id:
        raise ValueError(f"No valid source regions found in {source_csv}")

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

        
        parsed_phantom_code = match.group(1)
        phantom_spec = get_phantom(parsed_phantom_code)

        phantom_code = phantom_spec.code
        phantom = phantom_spec.display_name

        source_organ = match.group(2)
        source_type = source_type_map.get(
            match.group(3).lower(),
            match.group(3).lower()
        )
        source_energy = float(match.group(4))

        number_of_particles = params["nps"]

        # The source acronym comes from the output filename and its ID(s)
        # come from source_organs.csv. Do not infer source IDs from the
        # target-organ database: a source may comprise multiple region IDs.
        try:
            source_organ_id = source_name_to_id[source_organ]
        except KeyError as exc:
            available = ", ".join(sorted(source_name_to_id))
            raise KeyError(
                f"Source acronym '{source_organ}' from '{deposit_file.name}' "
                f"was not found in {source_csv}. Available acronyms: {available}"
            ) from exc

        if phantom_code not in ORGANS:
            raise KeyError(
                f"Phantom '{phantom_code}' was not found in the organ database."
            )

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

    if not all_tallies:
        raise RuntimeError(
            "No Geant4 deposit CSV files were successfully processed."
        )

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