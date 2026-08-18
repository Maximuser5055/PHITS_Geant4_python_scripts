# This script updates the persistent SAF database by concatenating
# the current PHITS or Geant4 SAF results with the existing database.
#
# The user is asked to inspect the updated database before the
# current result files are deleted.

import pandas as pd

import b_config.a_config as config

# ============================================================
# DATABASE FILES AND CONFIGS
# ============================================================

PHITS_DATABASE_FILE = (config.RESULTS_SAF_DATABASE_DIR / "phits_all_safs_uncertainties.csv")
GEANT4_DATABASE_FILE = (config.RESULTS_SAF_DATABASE_DIR / "geant4_all_safs_uncertainties.csv")

phits_results_dir = config.RESULTS_PHITS_DIR
geant4_results_dir = config.RESULTS_GEANT4_DIR
target_region_csv = config.TARGET_REGION_CSV
target_region_mapping = pd.read_csv(target_region_csv)

# ============================================================
# EXPECTED COLUMNS
# ============================================================

REQUIRED_COLUMNS = [
    "Phantom",

    "Source Organ ID",
    "Source Organ Name",
    "Source Type",
    "Source Energy (MeV)",

    "Target Organ IDs",
    "Target Region Name",
    "Target Region Mass (g)",

    "Dose (Gy/source)",
    "SAF (kg^-1)",

    "Relative Error",
    "Statistical Uncertainty (%)",

    "Calculation Method",
]

# ============================================================
# CURRENT RESULT FILES
# ============================================================

def find_current_result_files(simulation_code):

    simulation_code = simulation_code.upper()

    if simulation_code == "PHITS":

        root = config.RESULTS_PHITS_DIR
        pattern = "*_target_regions_dose_SAFs_*.csv"

    elif simulation_code == "GEANT4":

        root = config.RESULTS_GEANT4_DIR
        pattern = "*_target_regions_dose_SAFs_*.csv"

    else:

        raise ValueError(f"Unsupported simulation code: {simulation_code}")

    files = sorted(root.glob(pattern))

    return files

# ============================================================
# CHECK COLUMNS
# ============================================================

def validate_columns(df, filename):

    missing = [
        column
        for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\nMissing columns in {filename}:\n"
            +
            "\n".join(
                f"  {column}"
                for column in missing
            )
        )

# ============================================================
# REMOVE RESULTS ALREADY PRESENT IN DATABASE
# ============================================================

def filter_new_results(current_results,existing_database):

    # --------------------------------------------------------
    # Load target-region definitions
    # --------------------------------------------------------

    expected_target_regions = set(
        target_region_mapping[
            "Target region"
        ]
        .dropna()
        .astype(str)
    )

    # --------------------------------------------------------
    # Columns that uniquely identify a source simulation
    # --------------------------------------------------------

    simulation_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Type",
        "Source Energy (MeV)",
    ]

    # --------------------------------------------------------
    # If there is no existing database, everything is new
    # --------------------------------------------------------

    if existing_database.empty:

        return current_results.copy()

    new_rows = []

    skipped_groups = 0
    partial_groups = 0

    # --------------------------------------------------------
    # Group current results by source simulation
    # --------------------------------------------------------

    for group_key, current_group in current_results.groupby(simulation_columns, dropna=False):

        # ----------------------------------------------------
        # Find corresponding rows already in database
        # ----------------------------------------------------

        database_group = existing_database.copy()

        for column, value in zip(
            simulation_columns,
            group_key
        ):

            if pd.isna(value):

                database_group = database_group[
                    database_group[column].isna()
                ]

            else:

                database_group = database_group[
                    database_group[column] == value
                ]

        # ----------------------------------------------------
        # No existing results for this source simulation
        # ----------------------------------------------------

        if database_group.empty:

            new_rows.append(current_group)

            continue

        # ====================================================
        # Determine which target regions already exist
        # ====================================================

        existing_regions = set(
            database_group[
                "Target Region Name"
            ]
            .dropna()
            .astype(str)
        )

        # ----------------------------------------------------
        # Check whether all normal target regions exist
        # ----------------------------------------------------

        normal_target_regions = (
            expected_target_regions
            - {
                "Red (active) marrow",
                "50-um endosteal region",
            }
        )

        missing_normal_regions = (
            normal_target_regions
            - existing_regions
        )

        # ----------------------------------------------------
        # RBM/endosteum need to be checked by calculation
        # method because each can legitimately have two rows.
        # ----------------------------------------------------

        skeletal_regions = {
            "Red (active) marrow",
            "50-um endosteal region",
        }

        missing_skeletal_rows = []

        for region in skeletal_regions:

            region_current = current_group[
                current_group[
                    "Target Region Name"
                ].astype(str) == region
            ]

            # -----------------------------------------------
            # If current batch has no row for this region,
            # there is nothing to add.
            # -----------------------------------------------

            if region_current.empty:
                continue

            for method in region_current[
                "Calculation Method"
            ].dropna().unique():

                existing_method = database_group[
                    (
                        database_group[
                            "Target Region Name"
                        ].astype(str) == region
                    )
                    &
                    (
                        database_group[
                            "Calculation Method"
                        ].astype(str) == str(method)
                    )
                ]

                if existing_method.empty:

                    missing_skeletal_rows.append(
                        region_current[
                            region_current[
                                "Calculation Method"
                            ].astype(str) == str(method)
                        ]
                    )

        # ====================================================
        # Determine whether the entire source simulation
        # is already complete
        # ====================================================

        if (
            not missing_normal_regions
            and not missing_skeletal_rows
        ):

            print(
                f"\n[SKIP] Existing database already contains "
                f"complete results:"
            )

            print(
                f"       Phantom       : {group_key[0]}"
            )

            print(
                f"       Source Organ  : {group_key[1]}"
            )

            print(
                f"       Source Type   : {group_key[2]}"
            )

            print(
                f"       Source Energy : {group_key[3]} MeV"
            )

            skipped_groups += 1

            continue

        # ====================================================
        # Existing group is incomplete
        # ====================================================

        partial_groups += 1

        # ----------------------------------------------------
        # Add missing normal target regions
        # ----------------------------------------------------

        if missing_normal_regions:

            missing_normal_rows = current_group[
                current_group[
                    "Target Region Name"
                ].astype(str).isin(
                    missing_normal_regions
                )
            ]

            if not missing_normal_rows.empty:

                new_rows.append(
                    missing_normal_rows
                )

        # ----------------------------------------------------
        # Add missing RBM/endosteum calculation-method rows
        # ----------------------------------------------------

        if missing_skeletal_rows:

            new_rows.extend(
                missing_skeletal_rows
            )

    # --------------------------------------------------------
    # Combine new rows
    # --------------------------------------------------------

    if not new_rows:

        return pd.DataFrame(
            columns=current_results.columns
        )

    filtered_results = pd.concat(
        new_rows,
        ignore_index=True
    )

    print()
    print(f"Existing complete simulations skipped : {skipped_groups}")

    print(f"Partially existing simulations updated : {partial_groups}")

    print(f"New rows to append : {len(filtered_results)}")

    return filtered_results

# ============================================================
# UPDATE SAF DATABASE
# ============================================================

def update_master_saf_database(params):

    simulation_code = (params["simulation_code"].upper())

    # --------------------------------------------------------
    # Select current-results directory
    # and master database
    # --------------------------------------------------------

    if simulation_code == "PHITS":

        current_results_dir = phits_results_dir
        database_file = PHITS_DATABASE_FILE

    elif simulation_code == "GEANT4":

        current_results_dir = geant4_results_dir
        database_file = GEANT4_DATABASE_FILE

    else:

        raise ValueError(f"Unsupported simulation code: {simulation_code}")

    # --------------------------------------------------------
    # Find current result files
    # --------------------------------------------------------

    current_files = find_current_result_files(simulation_code)

    if not current_files:

        print()
        print(f"[WARNING] No current {simulation_code} SAF result files were found.")
        print(f"Search directory:\n{current_results_dir}")

        return None
    
    # --------------------------------------------------------
    # Print files being added
    # --------------------------------------------------------

    print()
    print("=" * 90)
    print(f"UPDATING {simulation_code} SAF DATABASE")
    print("=" * 90)
    print()
    print("Current result files:")

    for file in current_files:

        print(f"  {file.name}")

    # --------------------------------------------------------
    # Read current results
    # --------------------------------------------------------

    current_dataframes = []

    for file in current_files:

        df = pd.read_csv(file)

        validate_columns(df, file.name)

        current_dataframes.append(df)

    current_results = pd.concat(current_dataframes, ignore_index=True)

    # --------------------------------------------------------
    # Read existing database if it exists
    # --------------------------------------------------------

    if database_file.exists():

        print()
        print("Existing SAF database found:")

        print(f"  {database_file}")

        existing_database = pd.read_csv(database_file)

        validate_columns(existing_database, database_file.name)

    else:

        print()
        print("No existing SAF database found.")

        print("A new database will be created.")

        existing_database = pd.DataFrame(columns=REQUIRED_COLUMNS)

    # --------------------------------------------------------
    # Remove results that are already represented in the
    # existing database.
    #
    # For RBM and endosteum, Calculation Method is also checked
    # because each can have:
    #   Direct dose calculation
    #   Fluence-to-dose response functions
    # --------------------------------------------------------

    filtered_current_results = filter_new_results(current_results, existing_database)

    # --------------------------------------------------------
    # Concatenate only genuinely new rows
    # --------------------------------------------------------

    combined_database = pd.concat(
        [existing_database,
         filtered_current_results
        ], ignore_index=True)

    # --------------------------------------------------------
    # Exact duplicate protection
    # --------------------------------------------------------

    rows_before = len(combined_database)

    combined_database = (
        combined_database
        .drop_duplicates()
        .copy()
    )

    duplicate_rows_removed = (
        rows_before
        - len(combined_database)
    )

    # --------------------------------------------------------
    # Sort database using the same target-region order as the
    # TARGET_REGION_CSV, matching the target-region
    # generation script.
    # --------------------------------------------------------

    target_region_order = (
        target_region_mapping[
            "Target region"
        ]
        .tolist()
    )

    combined_database["Target Region Name"] = pd.Categorical(
        combined_database["Target Region Name"],
        categories=target_region_order,
        ordered=True,
    )

    combined_database["Calculation Method"] = pd.Categorical(
        combined_database["Calculation Method"],
        categories=[
            "Direct dose calculation",
            "Fluence-to-dose response functions",
        ],
        ordered=True,
    )

    combined_database.sort_values(
        by=[
            "Phantom",
            "Source Organ ID",
            "Source Type",
            "Source Energy (MeV)",
            "Target Region Name",
            "Calculation Method",
        ],
        inplace=True,
        ignore_index=True,
    )

    # --------------------------------------------------------
    # Write database
    # --------------------------------------------------------

    combined_database.to_csv(database_file, index=False)

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("-" * 90)

    print(f"Current result rows : {len(current_results)}")

    print(f"Existing database rows : {len(existing_database)}"
    )

    print(f"Duplicate rows removed : {duplicate_rows_removed}"
    )

    print(f"Final database rows : {len(combined_database)}"
    )

    print()
    print("SAF database updated:")

    print(f"  {database_file}")

    # --------------------------------------------------------
    # Basic database summary
    # --------------------------------------------------------

    print()
    print("Database summary:")

    print(f"  Phantoms      : "
          f"{combined_database['Phantom'].nunique()}"
    )

    print(f"  Source organs : "
          f"{combined_database['Source Organ ID'].nunique()}"
    )

    print(f"  Source types  : "
          f"{combined_database['Source Type'].nunique()}"
    )

    print(f"  Energies      : "
          f"{combined_database['Source Energy (MeV)'].nunique()}"
    )

    print(f"  Target regions: "
          f"{combined_database['Target Region Name'].nunique()}"
    )

    # --------------------------------------------------------
    # User inspection
    # --------------------------------------------------------

    print()
    print("=" * 90)
    print("DATABASE INSPECTION")
    print("=" * 90)

    print()
    print("Please inspect the updated SAF database before "
          "deleting the current result files."
    )

    print()
    print(f"Database file:\n{database_file}")

    print()
    print("Current result files will NOT be deleted yet.")

    # --------------------------------------------------------
    # Ask whether to delete current results
    # --------------------------------------------------------

    while True:

        answer = input(
            "\nDelete the current result files? [y/n]: "
        ).strip().lower()

        if answer in ("y", "yes"):

            print()
            print("Deleting current result files...")

            deleted = 0

            for file in current_files:

                try:

                    file.unlink()
                    print(f"  Deleted: {file.name}")

                    deleted += 1

                except OSError as error:

                    print(f"  [WARNING] Could not delete {file.name}: {error}")

            print()
            print(f"Deleted {deleted} current result file(s).")

            break

        elif answer in ("", "n", "no"):

            print()
            print("Current result files were kept.")

            break

        else:

            print("Please enter y or n.")

    print()
    print("=" * 90)
    print(f"{simulation_code} SAF DATABASE UPDATE COMPLETE")
    print("=" * 90)

    return combined_database