"""
Create publishable SAF tables from the master SAF databases.

The master databases contain long-format data:

    Phantom
    Source Organ
    Source Type
    Source Energy
    Target Region
    SAF
    Statistical Uncertainty
    Calculation Method
    ...

This script converts them into OpenDose-style wide tables.

For each phantom and source type, two files are produced:

    *_saf.csv
        Specific absorbed fractions (kg^-1)

    *_std.csv
        Statistical uncertainties (%)

The first two columns are:

    Source Organ
    Target Organ

and subsequent columns are source energies.

Example:

    Source Organ,Target Organ,0.1 MeV,0.2 MeV,0.5 MeV,...
    Liver,Brain,...
    Liver,Liver,...
    Liver,Thyroid,...
    Thyroid,Brain,...
    Thyroid,Liver,...
    ...

RBM and endosteum have two calculation methods in the master
database. For the publishable database, the method selected
below is used.
"""

import pandas as pd

import b_config.a_config as config


# ============================================================
# SETTINGS
# ============================================================

# ------------------------------------------------------------
# Calculation method used for RBM and endosteum in the
# publishable SAF database.
#
# Options:
#
#   "Fluence-to-dose response functions"
#   "Direct dose calculation"
# ------------------------------------------------------------

PUBLISHABLE_SKELETAL_METHODS = {
    "photon": "Fluence-to-dose response functions",
    "gamma": "Fluence-to-dose response functions",

    "electron": "Direct dose calculation",
    "e-": "Direct dose calculation",
}

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

publishable_dir = config.RESULTS_PUBLISHABLE_SAF_DATABASE_DIR

# ------------------------------------------------------------
# Master SAF databases and other configs
# ------------------------------------------------------------

MASTER_DATABASES = {

    "PHITS":
        config.RESULTS_SAF_DATABASE_DIR
        / "a_phits_all_safs_and_uncertainties.csv",

    "GEANT4":
        config.RESULTS_SAF_DATABASE_DIR
        / "b_geant4_all_safs_and_uncertainties.csv",
}

target_region_csv = config.TARGET_REGION_CSV

# ============================================================
# REQUIRED COLUMNS
# ============================================================

REQUIRED_COLUMNS = [

    "Phantom",

    "Source Organ ID",
    "Source Organ Name",

    "Source Type",
    "Source Energy (MeV)",

    "Target Region Name",

    "SAF (kg^-1)",

    "Statistical Uncertainty (%)",

    "Calculation Method",
]


# ============================================================
# PHANTOM DISPLAY NAMES
# ============================================================

PHANTOM_NAMES = {

    "Adult Female":
        "Adult Female Reference Computational Phantom",

    "Adult Male":
        "Adult Male Reference Computational Phantom",
}


# ============================================================
# SOURCE TYPE DISPLAY NAMES
# ============================================================

SOURCE_TYPE_NAMES = {

    "photon":
        "photons",

    "gamma":
        "photons",

    "electron":
        "electrons",

    "e-":
        "electrons",
}


# ============================================================
# FILE NAME COMPONENTS
# ============================================================

PHANTOM_FILE_NAMES = {

    "Adult Female":
        "adult_female",

    "Adult Male":
        "adult_male",
}


SOURCE_TYPE_FILE_NAMES = {

    "photon":
        "photons",

    "gamma":
        "photons",

    "electron":
        "electrons",

    "e-":
        "electrons",
}


# ============================================================
# VALIDATE MASTER DATABASE
# ============================================================

def validate_master_database(
    database,
    database_file
):

    missing = [

        column

        for column in REQUIRED_COLUMNS

        if column not in database.columns

    ]

    if missing:

        raise ValueError(
            f"\nMissing columns in "
            f"{database_file.name}:\n"
            +
            "\n".join(
                f"  {column}"
                for column in missing
            )
        )


# ============================================================
# CLEAN SOURCE TYPE
# ============================================================

def normalize_source_type(
    source_type
):

    source_type = str(
        source_type
    ).strip().lower()

    if source_type not in SOURCE_TYPE_NAMES:

        raise ValueError(
            f"Unsupported source type: "
            f"{source_type}"
        )

    return source_type


# ============================================================
# SELECT SKELETAL CALCULATION METHOD
# ============================================================

def select_skeletal_method(dataframe, source_type):

    source_type = normalize_source_type(source_type)

    publishable_method = (
        PUBLISHABLE_SKELETAL_METHODS[
            source_type
        ]
    )

    skeletal_regions = {
        "Red (active) marrow",
        "50-um endosteal region",
    }

    result = dataframe.copy()

    skeletal_mask = (
        result[
            "Target Region Name"
        ]
        .isin(skeletal_regions)
    )

    # --------------------------------------------------------
    # Keep all non-skeletal target regions.
    # --------------------------------------------------------

    non_skeletal = result[
        ~skeletal_mask
    ].copy()

    # --------------------------------------------------------
    # Select the requested method for RBM/endosteum.
    # --------------------------------------------------------

    skeletal = result[
        skeletal_mask
    ].copy()

    if skeletal.empty:

        return result

    skeletal = skeletal[
        skeletal["Calculation Method"].astype(str)
        == publishable_method
    ].copy()

    # --------------------------------------------------------
    # Make sure the requested skeletal method exists.
    # --------------------------------------------------------

    expected_skeletal_regions = (
        set(
            result.loc[
                skeletal_mask,
                "Target Region Name"
            ]
            .dropna()
            .unique()
        )
    )

    found_skeletal_regions = (
        set(
            skeletal[
                "Target Region Name"
            ]
            .dropna()
            .unique()
        )
    )

    missing = (
        expected_skeletal_regions
        - found_skeletal_regions
    )

    if missing:

        raise ValueError(
            "The selected publishable calculation "
            "method is missing for skeletal regions:\n"
            +
            "\n".join(
                f"  {region}"
                for region in sorted(missing)
            )
            +
            "\n\nSelected method:\n"
            f"  {publishable_method}"
        )

    return pd.concat(
        [
            non_skeletal,
            skeletal,
        ],
        ignore_index=True
    )


# ============================================================
# CREATE ONE PUBLISHABLE TABLE
# ============================================================

def create_publishable_table(
    dataframe,
    value_column,
):

    # --------------------------------------------------------
    # Columns identifying one SAF matrix entry
    # --------------------------------------------------------

    index_columns = [
        "Source Organ ID",
        "Source Organ Name",
        "Target Region Name",
    ]

    # --------------------------------------------------------
    # Check for duplicate entries.
    #
    # After skeletal-method filtering there should be only
    # one value for each:
    #
    # Source Organ + Target Region + Energy
    # --------------------------------------------------------

    duplicate_counts = (
        dataframe
        .groupby(
            [
                "Source Organ ID",
                "Target Region Name",
                "Source Energy (MeV)",
            ],
            dropna=False
        )
        .size()
    )

    duplicates = duplicate_counts[
        duplicate_counts > 1
    ]

    if not duplicates.empty:

        raise ValueError(
            "Duplicate SAF entries remain after "
            "selecting the skeletal calculation method.\n"
            "Each Source Organ + Target Region + "
            "Energy combination must have exactly "
            "one value."
        )

    # --------------------------------------------------------
    # Pivot:
    #
    # rows    = Source Organ + Target Region
    # columns = Source Energy
    # values  = SAF or uncertainty
    # --------------------------------------------------------

    table = dataframe.pivot(
        index=index_columns,
        columns="Source Energy (MeV)",
        values=value_column,
    )

    table = table.reset_index()

    # --------------------------------------------------------
    # Rename columns
    # --------------------------------------------------------

    table.rename(
        columns={
            "Source Organ Name":
                "Source Organ",

            "Target Region Name":
                "Target Organ",
        },
        inplace=True,
    )

    table.drop(
    columns=[
        "Source Organ ID"
    ], inplace=True,)
    
    # --------------------------------------------------------
    # Sort source organs by Source Organ ID
    # --------------------------------------------------------

    source_order = (
        dataframe[
            [
                "Source Organ ID",
                "Source Organ Name",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            "Source Organ ID"
        )
    )

    source_order = (
        source_order[
            "Source Organ Name"
        ]
        .tolist()
    )

    table["Source Organ"] = pd.Categorical(
        table["Source Organ"],
        categories=source_order,
        ordered=True,
    )

    # --------------------------------------------------------
    # Sort target regions according to the configured
    # Filipino target-region mapping.
    # --------------------------------------------------------

    target_mapping = pd.read_csv(target_region_csv)

    target_order = (
        target_mapping[
            "Target region"
        ]
        .dropna()
        .astype(str)
        .tolist()
    )

    table["Target Organ"] = pd.Categorical(
        table["Target Organ"],
        categories=target_order,
        ordered=True,
    )

    table.sort_values(
        [
            "Source Organ",
            "Target Organ",
        ],
        inplace=True,
        ignore_index=True,
    )

    # --------------------------------------------------------
    # Remove internal sorting columns
    # --------------------------------------------------------

    table["Source Organ"] = (
        table["Source Organ"]
        .astype(str)
    )

    table["Target Organ"] = (
        table["Target Organ"]
        .astype(str)
    )

    # --------------------------------------------------------
    # Sort energy columns numerically
    # --------------------------------------------------------

    fixed_columns = [
        "Source Organ",
        "Target Organ",
    ]

    energy_columns = [
        column
        for column in table.columns
        if column not in fixed_columns
    ]

    energy_columns = sorted(
        energy_columns,
        key=float
    )

    table = table[
        fixed_columns
        + energy_columns
    ]

    # --------------------------------------------------------
    # Convert energy column names into readable labels
    #
    # Example:
    #
    # 10.0
    #
    # becomes:
    #
    # 10.0 MeV
    # --------------------------------------------------------

    renamed_columns = {
        energy:
            f"{float(energy):g} MeV"
        for energy in energy_columns
    }

    table.rename(
        columns=renamed_columns,
        inplace=True,
    )

    return table


# ============================================================
# WRITE PUBLISHABLE CSV
# ============================================================

def write_publishable_csv(
    table,
    output_file,
    phantom,
    source_type,
    source_energies,
    source_organs,
    target_organs,
    document_name,
):

    phantom_name = PHANTOM_NAMES.get(
        phantom,
        phantom
    )

    source_type_name = SOURCE_TYPE_NAMES.get(
        source_type,
        source_type
    )

    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    header = [

        "# Header",

        f"# Document: {document_name}",

        f"# Phantom type: {phantom_name}",

        f"# Source type: {source_type_name}",

        f"# Number of source energies: "
        f"{len(source_energies)}",

        f"# Number of source organs: "
        f"{len(source_organs)}",

        f"# Number of target organs: "
        f"{len(target_organs)}",

        "# End Header",
    ]

    # --------------------------------------------------------
    # Write metadata header and CSV table
    # --------------------------------------------------------

    with open(
        output_file,
        "w",
        encoding="utf-8",
        newline="",
    ) as file:

        for line in header:

            file.write(
                line
                + "\n"
            )

        table.to_csv(
            file,
            index=False,
        )


# ============================================================
# CREATE PUBLISHABLE DATABASE
# ============================================================

def create_publishable_saf_database(params):

    simulation_code = (
        params[
            "simulation_code"
        ]
        .upper()
    )

    if simulation_code not in MASTER_DATABASES:

        raise ValueError(
            f"Unsupported simulation code: "
            f"{simulation_code}"
        )

    master_database_file = (
        MASTER_DATABASES[
            simulation_code
        ]
    )

    # --------------------------------------------------------
    # Check master database
    # --------------------------------------------------------

    if not master_database_file.exists():

        raise FileNotFoundError(
            "Could not find master SAF database:\n"
            f"{master_database_file}"
        )
    
    # --------------------------------------------------------
    # Read master database
    # --------------------------------------------------------

    database = pd.read_csv(
        master_database_file
    )

    validate_master_database(
        database,
        master_database_file
    )

    print()
    print("=" * 90)
    print(
        f"CREATING PUBLISHABLE "
        f"{simulation_code} SAF DATABASE"
    )
    print("=" * 90)

    print()
    print(
        f"Master database:"
    )

    print(
        f"  {master_database_file}"
    )

    print()
    print(
        f"Output directory:"
    )

    print(
        f"  {publishable_dir}"
    )

    # --------------------------------------------------------
    # Generate one SAF and one uncertainty table for each
    # phantom/source-type combination.
    # --------------------------------------------------------

    generated_files = []

    for phantom in sorted(
        database[
            "Phantom"
        ]
        .dropna()
        .unique()
    ):

        for source_type in sorted(
            database[
                "Source Type"
            ]
            .dropna()
            .unique()
        ):

            source_type = normalize_source_type(
                source_type
            )

            publishable_method = (PUBLISHABLE_SKELETAL_METHODS[source_type])

            subset = database[
                (
                    database[
                        "Phantom"
                    ]
                    == phantom
                )
                &
                (
                    database[
                        "Source Type"
                    ]
                    .str.lower()
                    == source_type
                )
            ].copy()

            if subset.empty:

                continue

            # ------------------------------------------------
            # Select the final skeletal values.
            # ------------------------------------------------

            subset = select_skeletal_method(subset, source_type)

            # ------------------------------------------------
            # Source energies
            # ------------------------------------------------

            source_energies = sorted(
                subset[
                    "Source Energy (MeV)"
                ]
                .dropna()
                .unique()
                .astype(float)
            )

            # ------------------------------------------------
            # Source organs
            # ------------------------------------------------

            source_organs = (
                subset[
                    "Source Organ ID"
                ]
                .dropna()
                .unique()
            )

            # ------------------------------------------------
            # Target organs
            # ------------------------------------------------

            target_mapping = pd.read_csv(target_region_csv)

            target_organs = (
                target_mapping[
                    "Target region"
                ]
                .dropna()
                .astype(str)
                .unique()
            )

            # ------------------------------------------------
            # Create SAF table
            # ------------------------------------------------

            saf_table = create_publishable_table(
                subset,
                "SAF (kg^-1)"
            )

            # ------------------------------------------------
            # Create uncertainty table
            # ------------------------------------------------

            std_table = create_publishable_table(
                subset,
                "Statistical Uncertainty (%)"
            )

            phantom_filename = (
                PHANTOM_FILE_NAMES[
                    phantom
                ]
            )

            source_filename = (
                SOURCE_TYPE_FILE_NAMES[
                    source_type
                ]
            )

            # ------------------------------------------------
            # SAF filename
            # ------------------------------------------------

            saf_file = (
                publishable_dir
                /
                f"{phantom_filename}_"
                f"{source_filename}_saf.csv"
            )

            # ------------------------------------------------
            # Standard uncertainty filename
            # ------------------------------------------------

            std_file = (
                publishable_dir
                /
                f"{phantom_filename}_"
                f"{source_filename}_std.csv"
            )

            # ------------------------------------------------
            # Write SAF file
            # ------------------------------------------------

            write_publishable_csv(
                saf_table,
                saf_file,
                phantom,
                source_type,
                source_energies,
                source_organs,
                target_organs,
                "Specific Absorbed Fractions (kg^-1)",
            )

            # ------------------------------------------------
            # Write uncertainty file
            # ------------------------------------------------

            write_publishable_csv(
                std_table,
                std_file,
                phantom,
                source_type,
                source_energies,
                source_organs,
                target_organs,
                "Statistical Uncertainties (%)",
            )

            generated_files.extend(
                [
                    saf_file,
                    std_file,
                ]
            )

            print()
            print(
                f"Generated:"
            )

            print(
                f"  {saf_file.name}"
            )

            print(
                f"  {std_file.name}"
            )

            print(f"    Skeletal method: {publishable_method}")
            
            print(
                f"    Source organs : "
                f"{len(source_organs)}"
            )

            print(
                f"    Target organs : "
                f"{len(target_organs)}"
            )

            print(
                f"    Energies      : "
                f"{len(source_energies)}"
            )

            print(
                f"    SAF rows      : "
                f"{len(saf_table)}"
            )

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("=" * 90)
    print(
        "PUBLISHABLE SAF DATABASE COMPLETE"
    )
    print("=" * 90)

    print()
    print(
        f"Files generated : "
        f"{len(generated_files)}"
    )

    print(
        f"Output directory:"
    )

    print(
        f"  {publishable_dir}"
    )

    for file in generated_files:

        print(
            f"  {file.name}"
        )

    print()
    print(
        "The files contain:"
    )

    print(
        "  Column 1 : Source Organ"
    )

    print(
        "  Column 2 : Target Organ"
    )

    print(
        "  Remaining : Source energy columns"
    )

    return generated_files