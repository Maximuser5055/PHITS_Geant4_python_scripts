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

# Output files
PHITS_DATABASE_FILE = config.RESULTS_SAF_DATABASE_DIR / "a_phits_all_safs_and_uncertainties.csv"
GEANT4_DATABASE_FILE = config.RESULTS_SAF_DATABASE_DIR / "b_geant4_all_safs_and_uncertainties.csv"

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
    "Number of Particles",

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
        pattern = "*_phits_*_target_regions_dose_SAFs.csv"

    elif simulation_code == "GEANT4":

        root = config.RESULTS_GEANT4_DIR
        pattern = "*_geant4_*_target_regions_dose_SAFs.csv"

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

def filter_new_results(current_results, existing_database, override_duplicates=False):

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

        current_source_sets = (
            current_results[
                [
                    "Source Organ ID",
                    "Source Type",
                    "Source Energy (MeV)",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        )

        current_phantom_types = (
            current_results[
                "Phantom"
            ]
            .nunique()
        )

        statistics = {

            "current_rows":
                len(current_results),

            "current_source_sets":
                current_source_sets,

            "current_phantom_types":
                current_phantom_types,

            "new_source_sets":
                current_source_sets,

            "partial_source_sets":
                0,

            "skipped_source_sets":
                0,

            "rows_to_append":
                len(current_results),
        }

        return current_results.copy(), statistics

    new_rows = []

    new_source_sets = set()
    partial_source_sets = set()
    skipped_source_sets = set()

    # --------------------------------------------------------
    # Group current results by source simulation
    # --------------------------------------------------------

    for group_key, current_group in current_results.groupby(simulation_columns, dropna=False):

        # ----------------------------------------------------
        # Find corresponding rows already in database
        # ----------------------------------------------------

        database_group = existing_database.copy()

        for column, value in zip(simulation_columns, group_key):

            if pd.isna(value):

                database_group = database_group[database_group[column].isna()]

            else:

                database_group = database_group[database_group[column] == value]

        # ----------------------------------------------------
        # No existing results for this source simulation
        # ----------------------------------------------------

        if database_group.empty:

            new_rows.append(current_group)
            source_set_key = (group_key[1], group_key[2], group_key[3],)
            new_source_sets.add(source_set_key)

            continue

        # ====================================================
        # Determine whether the source + target-organ results
        # are already complete
        # ====================================================

        # ----------------------------------------------------
        # Target regions actually present in the CURRENT
        # result files for this source simulation
        # ----------------------------------------------------

        current_regions = set(
            current_group[
                "Target Region Name"
            ]
            .dropna()
            .astype(str)
        )

        # ----------------------------------------------------
        # Target regions already present in the DATABASE
        # for this same source simulation
        # ----------------------------------------------------

        existing_regions = set(
            database_group[
                "Target Region Name"
            ]
            .dropna()
            .astype(str)
        )

        # ----------------------------------------------------
        # Normal target regions that are present in the
        # current results but missing from the database
        # ----------------------------------------------------

        missing_normal_regions = (
            current_regions
            - {
                "Red (active) marrow",
                "50-um endosteal region",
            }
        ) - existing_regions

        # ----------------------------------------------------
        # RBM and endosteum:
        # check Target Region + Calculation Method
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

            # No current result for this skeletal region
            if region_current.empty:
                continue

            for method in (
                region_current[
                    "Calculation Method"
                ]
                .dropna()
                .astype(str)
                .unique()
            ):

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
                        ].astype(str) == method
                    )
                ]

                if existing_method.empty:

                    missing_skeletal_rows.append(
                        region_current[
                            region_current[
                                "Calculation Method"
                            ].astype(str) == method
                        ]
                    )

        # ====================================================
        # Determine whether the entire source simulation
        # is already complete
        # ====================================================

        if (not missing_normal_regions and not missing_skeletal_rows):

            source_set_key = (
                group_key[1],
                group_key[2],
                group_key[3],
            )

            # --------------------------------------------------------
            # Duplicate exists
            #
            # Replace existing results if requested.
            # Otherwise keep the existing database results.
            # --------------------------------------------------------

            if override_duplicates:

                print(
                    f"\n[REPLACE] Existing results will be "
                    f"replaced:"
                )

                print(
                    f"          Phantom       : {group_key[0]}"
                )

                print(
                    f"          Source Organ  : {group_key[1]}"
                )

                print(
                    f"          Source Type   : {group_key[2]}"
                )

                print(
                    f"          Source Energy : {group_key[3]} MeV"
                )

                new_rows.append(
                    current_group
                )

                partial_source_sets.add(
                    source_set_key
                )

            else:

                print(
                    f"\n[SKIP] Keeping existing results:"
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

                skipped_source_sets.add(
                    source_set_key
                )

            continue

        # ====================================================
        # Existing group is incomplete
        # ====================================================

        source_set_key = (group_key[1], group_key[2], group_key[3],)

        partial_source_sets.add(source_set_key)

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

        filtered_results = pd.DataFrame(
            columns=current_results.columns
        )

    else:

        filtered_results = pd.concat(
            new_rows,
            ignore_index=True
        )

    all_source_sets = (
        new_source_sets
        |
        partial_source_sets
        |
        skipped_source_sets
    )

    skipped_source_sets = skipped_source_sets - partial_source_sets

    new_source_sets = new_source_sets - partial_source_sets - skipped_source_sets

    statistics = {

        "current_rows": len(current_results),

        "current_source_sets": len(all_source_sets),

        "current_phantom_types": current_results["Phantom"].nunique(),

        "new_source_sets": len(new_source_sets),

        "partial_source_sets": len(partial_source_sets),

        "skipped_source_sets": len(skipped_source_sets),

        "rows_to_append": len(filtered_results),
    }

    return filtered_results, statistics

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
    # Check whether any complete source simulations already
    # exist in the database.
    # --------------------------------------------------------

    simulation_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Type",
        "Source Energy (MeV)",
    ]

    if not existing_database.empty:

        current_simulations = (
            current_results[
                simulation_columns
            ]
            .drop_duplicates()
        )

        existing_simulations = (
            existing_database[
                simulation_columns
            ]
            .drop_duplicates()
        )

        duplicate_simulations = (
            current_simulations.merge(
                existing_simulations,
                on=simulation_columns,
                how="inner"
            )
        )

    else:

        duplicate_simulations = pd.DataFrame(
            columns=simulation_columns
        )
        
    override_duplicates = False

    if not duplicate_simulations.empty:

        print()
        print("=" * 90)
        print("DUPLICATE SAF SIMULATIONS DETECTED")
        print("=" * 90)

        print()

        print(
            f"Found {len(duplicate_simulations)} "
            "source simulation(s) already in the database."
        )

        print()
        print(
            "These may be reruns intended to improve "
            "statistical uncertainty."
        )

        print()

        for _, row in duplicate_simulations.iterrows():

            print(
                f"  {row['Phantom']} | "
                f"{row['Source Organ ID']} | "
                f"{row['Source Type']} | "
                f"{row['Source Energy (MeV)']:g} MeV"
            )

        print()
        print(
            "[1] Keep existing database results"
        )

        print(
            "[2] Replace with new results"
        )

        while True:

            choice = input(
                "\nEnter 1 or 2: "
            ).strip()

            if choice == "1":

                override_duplicates = False

                print(
                    "\nKeeping existing database results."
                )

                break

            elif choice == "2":

                override_duplicates = True

                print(
                    "\nNew results will replace "
                    "the existing duplicates."
                )

                break

            else:

                print(
                    "Invalid choice. "
                    "Please enter 1 or 2."
                )
                
    # --------------------------------------------------------
    # Remove results that are already represented in the
    # existing database.
    #
    # For RBM and endosteum, Calculation Method is also checked
    # because each can have:
    #   Direct dose calculation
    #   Fluence-to-dose response functions
    # --------------------------------------------------------

    filtered_current_results, filter_statistics = (
        filter_new_results(current_results, 
                           existing_database, 
                           override_duplicates)
    )

    # --------------------------------------------------------
    # Remove existing simulations that are being replaced
    # --------------------------------------------------------

    database_to_keep = existing_database.copy()

    if override_duplicates:

        for _, key in duplicate_simulations.iterrows():

            mask = (
                (database_to_keep["Phantom"] == key["Phantom"])
                &
                (
                    database_to_keep["Source Organ ID"]
                    == key["Source Organ ID"]
                )
                &
                (
                    database_to_keep["Source Type"]
                    == key["Source Type"]
                )
                &
                (
                    database_to_keep["Source Energy (MeV)"]
                    == key["Source Energy (MeV)"]
                )
            )

            database_to_keep = (
                database_to_keep[
                    ~mask
                ]
            )

    # --------------------------------------------------------
    # Concatenate existing + new results
    # --------------------------------------------------------

    combined_database = pd.concat(
        [
            database_to_keep,
            filtered_current_results
        ],
        ignore_index=True
    )

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
    print("=" * 90)
    print("SAF DATABASE UPDATE SUMMARY")
    print("=" * 90)

    print()

    # --------------------------------------------------------
    # Current batch breakdown
    # --------------------------------------------------------

    print("Current batch:")

    batch_summary = (current_results[
            ["Phantom",
             "Source Type",
             "Source Energy (MeV)",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            ["Phantom",
             "Source Type",
             "Source Energy (MeV)",
            ]
        )
    )

    for _, row in batch_summary.iterrows():

        mask = (
            (current_results["Phantom"] == row["Phantom"])
            &
            (current_results["Source Type"] == row["Source Type"])
            &
            (current_results["Source Energy (MeV)"] == row["Source Energy (MeV)"])
        )

        source_count = (current_results.loc[mask,"Source Organ ID"].nunique())
        row_count = (current_results.loc[mask].shape[0])

        print(f"  {row['Phantom']} | {row['Source Type']} | {row['Source Energy (MeV)']:g} MeV")

        print(f"      Source organs : {source_count}")

        print(f"      Result rows   : {row_count}")

    print()
    
    print(f"Current result rows       : {filter_statistics['current_rows']}")

    print(f"Current source sets       : {filter_statistics['current_source_sets']}")

    print(f"Phantom types             : {filter_statistics['current_phantom_types']}")

    print(f"Existing database rows    : {len(existing_database)}")

    print()

    print(f"New source sets           : {filter_statistics['new_source_sets']}")

    print(f"Partial source sets       : {filter_statistics['partial_source_sets']}")

    print(f"Skipped source sets       : {filter_statistics['skipped_source_sets']}")

    print()

    print(f"Rows to append            : {filter_statistics['rows_to_append']}")

    print(f"Exact duplicates removed  : {duplicate_rows_removed}")

    print(f"Final database rows       : {len(combined_database)}")
    print()

    print("Database:")

    print(f"  {database_file}")
      
    # --------------------------------------------------------
    # Basic database summary
    # --------------------------------------------------------

    print()
    print("Database summary:")

    print(f"  Phantoms           : "
          f"{combined_database['Phantom'].nunique()}"
    )

    print(f"  Source organs      : "
          f"{combined_database['Source Organ ID'].nunique()}"
    )

    print(f"  Source types       : "
          f"{combined_database['Source Type'].nunique()}"
    )

    print(f"  Energies           : "
          f"{combined_database['Source Energy (MeV)'].nunique()}"
    )

    print(f"  Particle history counts    : "
      f"{combined_database['Number of Particles'].nunique()}"
    )
    
    print(f"  Target regions     : "
          f"{combined_database['Target Region Name'].nunique()}"
    )

    print(f"  Total database rows: {len(combined_database)}")

    print()
    print("=" * 90)
    print(f"{simulation_code} SAF DATABASE UPDATE COMPLETE")
    print("=" * 90)

    return combined_database