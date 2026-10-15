# This script extracts necessary information from the tally deposit PHITS output files to compute for the specific absorbed fraction (SAF) for each target organ.
# It extracts: region ID, dose, and relative error and it correlates the extracted region ID to b_organ_database.py to get the organ id, name and mass.
# It then calculates the SAFs for each target organ.

# Import necessary libraries
import re
import pandas as pd
import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom
from c_database.b_organ_database import ORGANS

def phits_calculate_dose_and_safs(params):

    # -------------------------------------------------------------------------
    # Constants
    # -------------------------------------------------------------------------

    MeV_to_J = config.MEV_TO_J

    # -------------------------------------------------------------------------
    # Root directories
    # -------------------------------------------------------------------------

    input_root = config.GENERATED_INPUTS_DIR
    output_root = config.RESULTS_PHITS_DIR

    # -------------------------------------------------------------------------
    # File names
    # -------------------------------------------------------------------------

    phits_saf_file = "c_phits_dose_and_SAFs.csv"
    
    # -------------------------------------------------------------------------
    # Regex patterns
    # -------------------------------------------------------------------------

    patterns = {

        "file_name":
            re.compile(r"phits_deposit_(.+?)_source_(.+?)_([A-Za-z0-9+-]+)_energy_([0-9Ee.+-]+)\.out",
                re.IGNORECASE),

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
    # Determine selected phantom family
    # -------------------------------------------------------------------------

    phantom_selection = params["phantom"]
    output_file = output_root / phits_saf_file

    # -------------------------------------------------------------------------
    # Find all deposit tally files
    # -------------------------------------------------------------------------

    deposit_files = sorted(input_root.rglob("phits_deposit_*_source_*.out"))

    print(f"Found {len(deposit_files)} deposit tally file(s).\n")

    # Load source-region metadata from the same CSV used by input generation.
    # "ID number(s)" is kept as text so multi-ID sources such as 9700_9900
    # remain intact instead of being forced into a single integer.
    source_csv = config.SOURCE_CSV
    if not source_csv.exists():
        raise FileNotFoundError(f"Source-organ CSV not found: {source_csv}")

    source_df = pd.read_csv(source_csv, dtype=str, encoding="utf-8-sig")
    required_columns = {"Acronym", "ID number(s)"}
    missing_columns = required_columns - set(source_df.columns)
    if missing_columns:
        raise ValueError(
            f"Source-organ CSV is missing column(s): {sorted(missing_columns)}. "
            f"Found: {list(source_df.columns)}"
        )

    source_df["Acronym"] = source_df["Acronym"].fillna("").str.strip()
    source_df["ID number(s)"] = source_df["ID number(s)"].fillna("").str.strip()
    source_df = source_df[
        (source_df["Acronym"] != "") & (source_df["ID number(s)"] != "")
    ]
    source_name_to_id = dict(zip(source_df["Acronym"], source_df["ID number(s)"]))

    if not source_name_to_id:
        raise ValueError(f"No valid source regions found in {source_csv}")

    # -------------------------------------------------------------------------
    # Extract one deposit tally
    # -------------------------------------------------------------------------

    def extract_deposit(deposit_file):

        filename_match = patterns["file_name"].search(deposit_file.name)

        if not filename_match:
            raise ValueError(f"Could not parse filename:\n{deposit_file.name}")


        parsed_phantom_code = filename_match.group(1)
        phantom_spec = get_phantom(parsed_phantom_code)

        # Use the canonical registry code for database lookups.
        phantom = phantom_spec.code

        source_organ = filename_match.group(2)
        source_type = filename_match.group(3)
        source_energy = float(filename_match.group(4))

        number_of_particles = params["maxcas"] * params["maxbch"]

        if phantom not in ORGANS:
            raise KeyError(
                f"Phantom '{phantom}' was not found in ORGANS."
            )

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

            target_organ = organ_database.get(region)

            # Skip regions not in organ database
            if target_organ is None:
                continue

            # Source name comes from the input filename and is the Acronym
            # in source_organs.csv. Its ID may contain multiple region IDs.
            try:
                source_organ_id = source_name_to_id[source_organ]
            except KeyError as exc:
                available = ", ".join(sorted(source_name_to_id))
                raise KeyError(
                    f"Source acronym '{source_organ}' from '{deposit_file.name}' "
                    f"was not found in {source_csv}. Available acronyms: {available}"
                ) from exc

            source_energy_joule = source_energy * MeV_to_J

            saf = dose / source_energy_joule
            
            results.append({

                "Phantom":
                    phantom_spec.display_name,

                "Source Organ ID":
                    source_organ_id,

                "Source Organ Name":
                    source_organ,

                "Source Type":
                    source_type,

                "Source Energy (MeV)":
                    source_energy,

                "Number of Particles":
                    number_of_particles,

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

    if not all_tallies:
        raise RuntimeError("No PHITS deposit tally files were successfully processed.")

    combined_df = pd.concat(all_tallies, ignore_index=True)

    sort_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Type",
        "Source Energy (MeV)",
        "Target Organ ID"
    ]

    combined_df.sort_values(by=sort_columns, inplace=True)

    combined_df.reset_index(drop=True, inplace=True)

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    combined_df.to_csv(output_file, index=False)

    print("\nFinished extracting all deposit tallies.")