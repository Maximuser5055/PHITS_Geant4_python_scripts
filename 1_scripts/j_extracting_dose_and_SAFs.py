# This script extracts necessary information from the tally deposit PHITS output files to compute for the specific absorbed fraction (SAF) for each target organ.
# It extracts: region ID, dose, and relative error and it correlates the extracted region ID to b_organ_database.py to get the organ id, name and mass.
# It then calculates the SAFs for each target organ.

# Import necessary libraries
import re
import pandas as pd
import c_config as config
from e_organ_database import ORGANS

def calculate_dose_and_safs():
    # -------------------------------------------------------------------------
    # Root directories
    # -------------------------------------------------------------------------

    input_root = config.GENERATED_INPUTS_DIR

    output_root = config.RESULTS_DIR

    # -------------------------------------------------------------------------
    # Regex patterns
    # -------------------------------------------------------------------------

    patterns = {

        "file_name":
            re.compile(
                    r"deposit_MRCP_(AM|AF)_source_(.+?)_([A-Za-z0-9+-]+)_energy_([0-9Ee.+-]+)\.out",
                re.IGNORECASE
            ),

        "tally_row":
            re.compile(
                r"^\s*"
                r"(\d+)\s+"                 # row number
                r"(\d+)\s+"                 # region number
                r"([0-9Ee.+-]+)\s+"         # volume
                r"([0-9Ee.+-]+)\s+"         # y(all)
                r"([0-9Ee.+-]+)\s*$"        # r.err
            )
    }

    # -------------------------------------------------------------------------
    # Find all deposit tally files
    # -------------------------------------------------------------------------

    deposit_files = sorted(input_root.rglob("deposit_*.out"))

    print(f"Found {len(deposit_files)} deposit tally file(s).\n")

    # Find source organ ID
    ORGAN_NAME_TO_ID = {
        phantom: {
            organ["name"]: organ["organ_id"]
            for organ in ORGANS[phantom].values()
        }
        for phantom in ORGANS
    }

    # -------------------------------------------------------------------------
    # Extract one deposit tally
    # -------------------------------------------------------------------------

    def extract_deposit(deposit_file):

        filename_match = patterns["file_name"].search(deposit_file.name)

        if not filename_match:
            raise ValueError(f"Could not parse filename:\n{deposit_file.name}")

        phantom = filename_match.group(1)
        source_organ = filename_match.group(2)
        source_type = filename_match.group(3)
        source_energy = float(filename_match.group(4))

        organ_database = ORGANS[phantom]

        results = []

        with open(deposit_file, encoding="utf-8", errors="ignore") as f:

            lines = f.readlines()

        for line in lines:

            match = patterns["tally_row"].match(line)

            if not match:
                continue

            region = int(match.group(2))

            dose = float(match.group(4))

            relative_error = float(match.group(5))

            source_organ_id = ORGAN_NAME_TO_ID[phantom][source_organ]

            target_organ = organ_database.get(region)

            MEV_TO_J = 1.6021766339999e-13

            source_energy_joule = source_energy * MEV_TO_J

            saf = dose / source_energy_joule

            # Skip regions not in organ database
            if target_organ is None:
                continue

            results.append({

                "Phantom":
                    "Adult Male" if phantom == "AM" else "Adult Female",

                "Source Organ ID":
                    source_organ_id,

                "Source Organ Name":
                    source_organ,

                "Source Type":
                    source_type,

                "Source Energy (MeV)":
                    source_energy,

                "Target Organ ID":
                    target_organ["organ_id"],

                "Target Organ Name":
                    target_organ["name"],

                "Target Organ Mass (g)":
                    target_organ["mass"],

                "Dose (Gy/source)":
                    dose,

                "SAF (kg^-1)":
                    saf,

                "Relative Error":
                    relative_error,

                "Statistical Uncertainty (%)":
                    relative_error * 100
            })

        return pd.DataFrame(results)

    # -------------------------------------------------------------------------
    # Extract every tally
    # -------------------------------------------------------------------------

    all_tallies = []

    for deposit_file in deposit_files:

        print(f"Processing {deposit_file.name}")

        df = extract_deposit(deposit_file)

        all_tallies.append(df)

    # -------------------------------------------------------------------------
    # Combine all tallies
    # -------------------------------------------------------------------------

    combined_df = pd.concat(all_tallies, ignore_index=True)

    # Separate Adult Male and Adult Female
    adult_male = combined_df[combined_df["Phantom"] == "Adult Male"].copy()
    adult_female = combined_df[combined_df["Phantom"] == "Adult Female"].copy()

    # Sort each dataframe
    sort_columns = [
        "Source Type",
        "Source Energy (MeV)",
        "Source Organ ID",
        "Target Organ ID"
    ]

    adult_male = adult_male.sort_values(by=sort_columns)
    adult_female = adult_female.sort_values(by=sort_columns)

    adult_male.reset_index(drop=True, inplace=True)
    adult_female.reset_index(drop=True, inplace=True)

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    adult_male.to_csv(
        output_root / "3_all_dose_and_SAFs_AM.csv",
        index=False
    )

    adult_female.to_csv(
        output_root / "4_all_dose_and_SAFs_AF.csv",
        index=False
    )

    print("\nFinished extracting all deposit tallies.")