# ============================================================
# Concatenate SAF databases from different computers
#
# This script combines two master SAF databases produced by
# separate computers.
#
# The script:
#   1. Selects PHITS or Geant4
#   2. Loads the current computer's master SAF database
#   3. Loads a SAF database from another computer
#   4. Checks column compatibility
#   5. Identifies:
#        - new rows
#        - exact duplicates
#        - conflicting results
#   6. Allows the user to resolve conflicts
#   7. Creates a merged database
#   8. Allows inspection before replacing the master
#   9. Optionally creates the publishable SAF database
#
# IMPORTANT:
# Statistical uncertainty is NOT used to automatically choose
# between conflicting results.
#
# ============================================================


# ============================================================
# IMPORTS
# ============================================================
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

import b_config.a_config as config


# ============================================================
# CONSTANTS
# ============================================================

PHITS_DATABASE_NAME = "a_phits_all_safs_and_uncertainties.csv"
GEANT4_DATABASE_NAME = "b_geant4_all_safs_and_uncertainties.csv"

# These columns define one individual SAF calculation.
#
# Number of particles is intentionally NOT included.
# Therefore:
# 10^8 and 10^9 particles
# are considered the same source/target calculation and
# can be identified as a conflict if their results differ.

CALCULATION_KEY = [
    "Phantom",
    "Source Organ ID",
    "Source Type",
    "Source Energy (MeV)",
    "Target Region Name",
    "Calculation Method",
]

# Tolerance for comparing two SAF databases
NUMERICAL_RTOL = 1e-10
NUMERICAL_ATOL = 1e-30

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def print_separator():
    print("=" * 90)

def choose_simulation_code():

    print_separator()
    print("SELECT SAF DATABASE")
    print_separator()

    print()
    print("[1] PHITS")
    print("[2] Geant4")

    while True:

        choice = input(
            "\nEnter 1 or 2: "
        ).strip()

        if choice == "1":
            return "PHITS"

        if choice == "2":
            return "GEANT4"

        print(
            "\nInvalid choice. "
            "Please enter 1 or 2."
        )

def get_database_path(simulation_code):

    if simulation_code == "PHITS":

        database_path = (
            config.RESULTS_SAF_DATABASE_DIR
            / PHITS_DATABASE_NAME
        )

    elif simulation_code == "GEANT4":

        database_path = (
            config.RESULTS_SAF_DATABASE_DIR
            / GEANT4_DATABASE_NAME
        )

    else:

        raise ValueError(
            f"Unsupported simulation code: "
            f"{simulation_code}"
        )

    return database_path


def load_database(path, label):

    if not path.is_file():

        raise FileNotFoundError(
            f"\n{label} database was not found:\n"
            f"    {path}\n"
        )

    try:

        df = pd.read_csv(path)

    except Exception as error:

        raise RuntimeError(
            f"\nCould not read {label} database:\n"
            f"    {path}\n"
            f"\nError:\n"
            f"    {error}\n"
        )

    if df.empty:

        print(
            f"\nWarning: {label} database is empty:"
            f"\n    {path}"
        )

    return df


def check_columns(current_database, other_database):

    current_columns = list(
        current_database.columns
    )

    other_columns = list(
        other_database.columns
    )

    current_set = set(
        current_columns
    )

    other_set = set(
        other_columns
    )

    missing_from_other = (
        current_set - other_set
    )

    missing_from_current = (
        other_set - current_set
    )

    if (
        missing_from_other
        or missing_from_current
    ):

        print_separator()
        print("COLUMN MISMATCH")
        print_separator()

        if missing_from_other:

            print(
                "\nColumns present in the current "
                "database but missing from the "
                "other database:"
            )

            for column in sorted(
                missing_from_other
            ):

                print(
                    f"    - {column}"
                )

        if missing_from_current:

            print(
                "\nColumns present in the other "
                "database but missing from the "
                "current database:"
            )

            for column in sorted(
                missing_from_current
            ):

                print(
                    f"    - {column}"
                )

        print(
            "\nThe databases cannot be safely "
            "concatenated."
        )

        return False

    return True

def normalize_key_columns(df):

    df = df.copy()

    # --------------------------------------------------------
    # Normalize strings
    # --------------------------------------------------------

    string_columns = [
        "Phantom",
        "Source Type",
        "Target Region Name",
        "Calculation Method",
    ]

    for column in string_columns:

        if column in df.columns:

            df[column] = (
                df[column]
                .fillna("")
                .astype(str)
                .str.strip()
            )

    # --------------------------------------------------------
    # Normalize numeric columns
    # --------------------------------------------------------

    if "Source Organ ID" in df.columns:

        df["Source Organ ID"] = pd.to_numeric(
            df["Source Organ ID"],
            errors="coerce"
        )

    if "Source Energy (MeV)" in df.columns:

        df["Source Energy (MeV)"] = pd.to_numeric(
            df["Source Energy (MeV)"],
            errors="coerce"
        )

    return df

def validate_key_columns(df, label):

    missing = [
        column
        for column in CALCULATION_KEY
        if column not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\n{label} database is missing "
            f"required identification columns:\n"
            +
            "\n".join(
                f"    - {column}"
                for column in missing
            )
        )

def records_are_equivalent(
    current_rows,
    other_rows,
    comparison_columns
):
    """
    Determine whether two sets of records are equivalent.

    Numeric values are compared using a relative and absolute
    tolerance to avoid false conflicts caused by floating-point
    representation differences.

    Non-numeric values must match exactly.
    """

    current = (
        current_rows[
            comparison_columns
        ]
        .drop_duplicates()
        .reset_index(drop=True)
        .copy()
    )

    other = (
        other_rows[
            comparison_columns
        ]
        .drop_duplicates()
        .reset_index(drop=True)
        .copy()
    )

    # Different number of unique records means they are different.
    if len(current) != len(other):
        return False

    if current.empty:
        return True

    # --------------------------------------------------------
    # Sort rows so that row ordering does not matter.
    # --------------------------------------------------------

    current = current.sort_values(
        by=comparison_columns,
        kind="stable"
    ).reset_index(drop=True)

    other = other.sort_values(
        by=comparison_columns,
        kind="stable"
    ).reset_index(drop=True)

    # --------------------------------------------------------
    # Compare each column.
    # --------------------------------------------------------

    for column in comparison_columns:

        current_values = current[column]
        other_values = other[column]

        # ----------------------------------------------------
        # Numeric columns
        # ----------------------------------------------------

        if (
            pd.api.types.is_numeric_dtype(
                current_values
            )
            and
            pd.api.types.is_numeric_dtype(
                other_values
            )
        ):

            current_numeric = (
                pd.to_numeric(
                    current_values,
                    errors="coerce"
                )
                .to_numpy(
                    dtype=float
                )
            )

            other_numeric = (
                pd.to_numeric(
                    other_values,
                    errors="coerce"
                )
                .to_numpy(
                    dtype=float
                )
            )

            if not (
                pd.Series(
                    current_numeric
                ).combine(
                    pd.Series(
                        other_numeric
                    ),
                    lambda current_value,
                           other_value:
                        (
                            pd.isna(
                                current_value
                            )
                            and
                            pd.isna(
                                other_value
                            )
                        )
                        or
                        (
                            pd.notna(
                                current_value
                            )
                            and
                            pd.notna(
                                other_value
                            )
                            and
                            abs(
                                current_value
                                - other_value
                            )
                            <=
                            (
                                NUMERICAL_ATOL
                                +
                                NUMERICAL_RTOL
                                *
                                max(
                                    abs(
                                        current_value
                                    ),
                                    abs(
                                        other_value
                                    )
                                )
                            )
                        )
                )
                .all()
            ):

                return False

        # ----------------------------------------------------
        # Non-numeric columns
        # ----------------------------------------------------

        else:

            current_values = (
                current_values
                .fillna("")
                .astype(str)
                .to_numpy()
            )

            other_values = (
                other_values
                .fillna("")
                .astype(str)
                .to_numpy()
            )

            if not (
                current_values
                ==
                other_values
            ).all():

                return False

    return True

def find_conflicting_calculations(current_database, other_database):

    # --------------------------------------------------------
    # Identify calculations appearing in both databases.
    # --------------------------------------------------------

    current_keys = (
        current_database[
            CALCULATION_KEY
        ]
        .drop_duplicates()
    )

    other_keys = (
        other_database[
            CALCULATION_KEY
        ]
        .drop_duplicates()
    )

    common_keys = current_keys.merge(
        other_keys,
        on=CALCULATION_KEY,
        how="inner"
    )

    if common_keys.empty:

        return pd.DataFrame(
            columns=CALCULATION_KEY
        )

    # --------------------------------------------------------
    # Columns to compare.
    #
    # The calculation key identifies the calculation.
    # Everything else is checked for consistency.
    # --------------------------------------------------------

    comparison_columns = [
        column
        for column in current_database.columns
        if column not in CALCULATION_KEY
    ]

    conflicts = []

    # --------------------------------------------------------
    # Compare each common calculation.
    # --------------------------------------------------------

    for _, key in common_keys.iterrows():

        current_mask = pd.Series(
            True,
            index=current_database.index
        )

        other_mask = pd.Series(
            True,
            index=other_database.index
        )

        for column in CALCULATION_KEY:

            current_mask &= (
                current_database[column]
                == key[column]
            )

            other_mask &= (
                other_database[column]
                == key[column]
            )

        current_rows = (
            current_database[
                current_mask
            ]
        )

        other_rows = (
            other_database[
                other_mask
            ]
        )

        if not records_are_equivalent(
            current_rows,
            other_rows,
            comparison_columns
        ):

            conflicts.append(
                {
                    **key.to_dict()
                }
            )

    return pd.DataFrame(
        conflicts,
        columns=CALCULATION_KEY
    )

def find_new_rows(current_database, other_database):

    current_keys = set(
        map(
            tuple,
            current_database[
                CALCULATION_KEY
            ].drop_duplicates().to_numpy()
        )
    )

    other_keys = (
        other_database[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
    )

    new_mask = ~other_keys.isin(
        current_keys
    )

    return other_database[
        new_mask
    ].copy()

def find_exact_duplicate_rows(current_database, other_database):

    # --------------------------------------------------------
    # Find rows from the other database that already exist
    # exactly in the current database.
    # --------------------------------------------------------

    current_records = set(
        tuple(row)
        for row in current_database.to_numpy()
    )

    duplicate_mask = (
        other_database.apply(
            lambda row:
                tuple(row)
                in current_records,
            axis=1
        )
    )

    return other_database[
        duplicate_mask
    ].copy()


def show_conflicts(conflicts):

    if conflicts.empty:

        print(
            "\nNo conflicting calculations "
            "were found."
        )

        return

    print_separator()
    print("CONFLICTING CALCULATIONS")
    print_separator()

    print(
        f"\nNumber of conflicting calculations: "
        f"{len(conflicts)}"
    )

    print(
        "\nThese calculations exist in both "
        "databases but have different result "
        "values or metadata."
    )

    print(
        "\nStatistical uncertainty and particle "
        "count are NOT used to automatically "
        "choose a result."
    )

    print()

    for _, row in conflicts.iterrows():

        energy = row[
            "Source Energy (MeV)"
        ]

        try:
            energy_text = (
                f"{float(energy):g}"
            )

        except Exception:
            energy_text = str(
                energy
            )

        print(
            f"  {row['Phantom']} | "
            f"Source {row['Source Organ ID']} | "
            f"{row['Source Type']} | "
            f"{energy_text} MeV | "
            f"{row['Target Region Name']} | "
            f"{row['Calculation Method']}"
        )


def choose_conflict_resolution():

    print()

    print(
        "[1] Keep the current computer's results"
    )

    print(
        "[2] Use the other computer's results"
    )

    print(
        "[3] Cancel the merge"
    )

    while True:

        choice = input(
            "\nChoose 1, 2, or 3: "
        ).strip()

        if choice == "1":
            return "CURRENT"

        if choice == "2":
            return "OTHER"

        if choice == "3":
            return "CANCEL"

        print(
            "Invalid choice. "
            "Please enter 1, 2, or 3."
        )


def replace_conflicting_rows(
    current_database,
    other_database,
    conflicts,
    resolution
):

    if conflicts.empty:

        return (
            current_database,
            other_database
        )

    if resolution == "CURRENT":

        # ----------------------------------------------------
        # Remove conflicting rows from the other database.
        # ----------------------------------------------------

        other_database = (
            remove_calculation_keys(
                other_database,
                conflicts
            )
        )

    elif resolution == "OTHER":

        # ----------------------------------------------------
        # Remove conflicting rows from the current database.
        # ----------------------------------------------------

        current_database = (
            remove_calculation_keys(
                current_database,
                conflicts
            )
        )

    return (
        current_database,
        other_database
    )


def remove_calculation_keys(
    database,
    keys
):

    if keys.empty:

        return database

    merged = database.merge(
        keys[
            CALCULATION_KEY
        ].drop_duplicates(),
        on=CALCULATION_KEY,
        how="left",
        indicator=True
    )

    return (
        merged[
            merged["_merge"] == "left_only"
        ]
        .drop(
            columns=["_merge"]
        )
        .copy()
    )


def sort_database(database):

    sort_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Type",
        "Source Energy (MeV)",
        "Target Region Name",
    ]

    available_columns = [
        column
        for column in sort_columns
        if column in database.columns
    ]

    return database.sort_values(
        by=available_columns,
        kind="stable"
    ).reset_index(
        drop=True
    )


def save_database(
    database,
    output_path
):

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    database.to_csv(
        output_path,
        index=False
    )


def ask_replace_master():

    print()

    while True:

        choice = input(
            "\nReplace the current master database "
            "with the merged database? [y/N]: "
        ).strip().lower()

        if choice in {
            "",
            "n",
            "no"
        }:

            return False

        if choice in {
            "y",
            "yes"
        }:

            return True

        print(
            "Please enter y or n."
        )


# ============================================================
# MAIN FUNCTION
# ============================================================

def concatenate_master_saf_databases():

    # ========================================================
    # SELECT SIMULATION CODE
    # ========================================================

    simulation_code = (
        choose_simulation_code()
    )

    current_database_path = (
        get_database_path(
            simulation_code
        )
    )

    print()

    print(
        f"Current {simulation_code} "
        "master database:"
    )

    print(
        f"    {current_database_path}"
    )

    # ========================================================
    # LOAD CURRENT DATABASE
    # ========================================================

    current_database = load_database(
        current_database_path,
        "Current"
    )

    # ========================================================
    # ASK FOR OTHER COMPUTER DATABASE
    # ========================================================

    print()

    print(
        "Enter the full path to the "
        f"{simulation_code} SAF database "
        "from the other computer."
    )

    print()

    while True:

        other_path_input = (
            input("Other database path: ")
            .strip()
            .strip('"')
            .strip("'")
            .replace("\\", "/")
        )

        if not other_path_input:

            print("\nA database path is required.")

            continue

        other_database_path = Path(
            other_path_input
        ).expanduser()

        if other_database_path.is_file():

            break

        print(
            f"\nFile not found:"
            f"\n    {other_database_path}"
        )

        print(
            "Please try again."
        )

    # ========================================================
    # LOAD OTHER DATABASE
    # ========================================================

    other_database = load_database(
        other_database_path,
        "Other"
    )

    # ========================================================
    # VALIDATE COLUMNS
    # ========================================================

    if not check_columns(
        current_database,
        other_database
    ):

        return

    # ========================================================
    # VALIDATE CALCULATION KEY
    # ========================================================

    validate_key_columns(
        current_database,
        "Current"
    )

    validate_key_columns(
        other_database,
        "Other"
    )

    # ========================================================
    # NORMALIZE DATA
    # ========================================================

    current_database = (
        normalize_key_columns(
            current_database
        )
    )

    other_database = (
        normalize_key_columns(
            other_database
        )
    )

    # ========================================================
    # BASIC STATISTICS
    # ========================================================

    current_rows = len(
        current_database
    )

    other_rows = len(
        other_database
    )

    # ========================================================
    # CLASSIFY OTHER DATABASE ROWS
    # ========================================================

    # --------------------------------------------------------
    # New calculations:
    # The calculation key does not exist in the current
    # database.
    # --------------------------------------------------------

    new_rows = find_new_rows(
        current_database,
        other_database
    )

    # --------------------------------------------------------
    # Conflicting calculations:
    # The calculation key exists in both databases, but the
    # associated values/metadata are meaningfully different.
    # --------------------------------------------------------

    conflicts = (
        find_conflicting_calculations(
            current_database,
            other_database
        )
    )

    # --------------------------------------------------------
    # Equivalent calculations:
    # The calculation key exists in both databases and the
    # records are numerically/textually equivalent within the
    # defined tolerance.
    # --------------------------------------------------------

    new_calculation_keys = set(
        new_rows[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
    )

    conflict_keys = set(
        conflicts[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
    )

    other_keys = set(
        other_database[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
    )

    duplicate_calculation_keys = (
        other_keys
        - new_calculation_keys
        - conflict_keys
    )

    equivalent_rows = other_database[
        other_database[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
        .isin(
            duplicate_calculation_keys
        )
    ].copy()

    # ========================================================
    # DISPLAY SUMMARY
    # ========================================================

    print()

    print_separator()
    print("MERGE SUMMARY")
    print_separator()

    print(
        f"\nCurrent database rows : "
        f"{current_rows}"
    )

    print(
        f"Other database rows   : "
        f"{other_rows}"
    )

    print(
        f"\nNew calculation rows  : "
        f"{len(new_rows)}"
    )

    print(
        f"Equivalent calculation rows: "
        f"{len(equivalent_rows)}"
    )

    print(
        f"Conflicting calculations: "
        f"{len(conflicts)}"
    )

    # ========================================================
    # SHOW CONFLICTS
    # ========================================================

    show_conflicts(
        conflicts
    )

    # ========================================================
    # RESOLVE CONFLICTS
    # ========================================================

    conflict_resolution = (
        "NONE"
    )

    if not conflicts.empty:

        conflict_resolution = (
            choose_conflict_resolution()
        )

        if conflict_resolution == "CANCEL":

            print(
                "\nMerge cancelled."
            )

            return

    # ========================================================
    # REMOVE / REPLACE CONFLICTING RESULTS
    # ========================================================

    (
        current_for_merge,
        other_for_merge
    ) = replace_conflicting_rows(
        current_database,
        other_database,
        conflicts,
        conflict_resolution
    )

    # ========================================================
    # REMOVE EQUIVALENT RESULTS FROM OTHER DATABASE
    #
    # Equivalent calculations already exist in the current
    # database. They must NOT be concatenated again.
    #
    # This uses the same calculation-key classification that
    # produced "Equivalent calculation rows" above.
    # ========================================================

    equivalent_keys = set(
        equivalent_rows[
            CALCULATION_KEY
        ]
        .apply(
            tuple,
            axis=1
        )
    )

    if equivalent_keys:

        other_keys_for_merge = (
            other_for_merge[
                CALCULATION_KEY
            ]
            .apply(
                tuple,
                axis=1
            )
        )

        equivalent_mask = (
            other_keys_for_merge
            .isin(
                equivalent_keys
            )
        )

        equivalent_rows_removed = int(
            equivalent_mask.sum()
        )

        other_for_merge = (
            other_for_merge[
                ~equivalent_mask
            ]
            .copy()
        )

    else:

        equivalent_rows_removed = 0

    # ========================================================
    # CONCATENATE
    # ========================================================

    combined_database = pd.concat(
        [
            current_for_merge,
            other_for_merge
        ],
        ignore_index=True
    )

    rows_before_exact_deduplication = (
        len(combined_database)
    )

    # ========================================================
    # REMOVE EXACT DUPLICATES
    #
    # This is now only a final safety check.
    # Equivalent rows have already been removed above.
    # ========================================================

    combined_database = (
        combined_database
        .drop_duplicates()
        .reset_index(
            drop=True
        )
    )

    exact_duplicates_removed = (
        rows_before_exact_deduplication
        - len(combined_database)
    )

    # ========================================================
    # SORT
    # ========================================================

    combined_database = (
        sort_database(
            combined_database
        )
    )

    # ========================================================
    # SAVE TEMPORARY MERGED DATABASE
    # ========================================================

    merged_database_path = (
        current_database_path.parent
        /
        (
            current_database_path.stem
            +
            "_merged.csv"
        )
    )

    save_database(
        combined_database,
        merged_database_path
    )

    # ========================================================
    # FINAL MERGE SUMMARY
    # ========================================================

    print()

    print_separator()
    print("MERGED DATABASE CREATED")
    print_separator()

    print(
        f"\nCurrent database rows       : "
        f"{current_rows}"
    )

    print(
        f"Other database rows         : "
        f"{other_rows}"
    )

    print(
        f"Rows before deduplication   : "
        f"{rows_before_exact_deduplication}"
    )

    print(
        f"Equivalent rows removed     : "
        f"{equivalent_rows_removed}"
    )

    print(
        f"Exact duplicate rows removed: "
        f"{exact_duplicates_removed}"
    )

    print(
        f"Final merged database rows  : "
        f"{len(combined_database)}"
    )

    print(
        f"\nMerged database:"
        f"\n    {merged_database_path}"
    )

    # ========================================================
    # INSPECTION
    # ========================================================

    print()

    print(
        "The merged database has been saved."
    )

    print(
        "Please inspect it before replacing "
        "the current master database."
    )

    print()

    print(
        "Open the following file:"
    )

    print(
        f"    {merged_database_path}"
    )

    print()

    input(
        "Press ENTER after you have inspected "
        "the merged database..."
    )

    # ========================================================
    # REPLACE CURRENT MASTER
    # ========================================================

    replace_master = (
        ask_replace_master()
    )

    if not replace_master:

        print()

        print(
            "Current master database was NOT "
            "modified."
        )

        print(
            f"\nThe merged database remains at:"
            f"\n    {merged_database_path}"
        )

        return

    # ========================================================
    # REPLACE MASTER DATABASE
    # ========================================================

    try:

        combined_database.to_csv(
            current_database_path,
            index=False
        )

    except Exception as error:

        print()

        print(
            "ERROR: Could not replace the "
            "current master database."
        )

        print(
            f"\n{error}"
        )

        print(
            f"\nYour merged database is still safely "
            f"available at:"
            f"\n    {merged_database_path}"
        )

        return

    # ========================================================
    # REMOVE TEMPORARY MERGED FILE
    # ========================================================

    try:

        merged_database_path.unlink()

    except Exception:

        print(
            "\nWarning: Could not remove temporary "
            "merged database:"
        )

        print(
            f"    {merged_database_path}"
        )

    # ========================================================
    # SUCCESS
    # ========================================================

    print()

    print_separator()
    print("MASTER SAF DATABASE UPDATED")
    print_separator()

    print()

    print(
        f"Updated database:"
    )

    print(
        f"    {current_database_path}"
    )

    print()

    print(
        f"Final rows:"
        f"    {len(combined_database)}"
    )

    print()

    try:

        from f_simulation_and_SAFs_further_analysis.e_create_and_update_publishable_saf_database import (
            create_publishable_saf_database
        )

        create_publishable_saf_database(
            {
                "simulation_code":
                    simulation_code
            }
        )

    except Exception as error:

        print()

        print(
            "The master SAF database was updated, "
            "but publishable SAF database generation "
            "failed."
        )

        print(
            f"\nError:"
            f"\n    {error}"
        )

# ============================================================
# MAIN
# ============================================================

def main():

    try:

        concatenate_master_saf_databases()

    except KeyboardInterrupt:

        print(
            "\n\nOperation cancelled by user."
        )

    except Exception as error:

        print()

        print_separator()
        print("ERROR")
        print_separator()

        print()

        print(error)

        print()


if __name__ == "__main__":

    main()