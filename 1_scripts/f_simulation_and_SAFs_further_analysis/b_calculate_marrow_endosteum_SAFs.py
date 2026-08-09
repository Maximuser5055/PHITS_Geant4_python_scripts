"""
Calculate active red bone marrow (RBM) dose and SAF using:

1. Geant4 energy-binned photon fluence in skeletal spongiosa.
2. ICRP Publication 116, Table D.1 AM response functions.

Geant4 CSV columns:
    Spongiosa ID
    Energy Low (MeV)
    Energy High (MeV)
    Energy Center (MeV)
    Fluence (photons/m2/source)

ICRP CSV:
    Photon energy (MeV)
    AM response-function columns for the ICRP skeletal IDs.

The Geant4 fluence values are histogram/bin fluences,
so they are NOT multiplied by Delta E again.

ICRP response functions are interpolated in log(E)-log(R)
space.

The final RBM dose is mass-weighted over the 13 skeletal
spongiosa regions.

For a 1 MeV photon source:

    SAF = RBM dose / emitted energy

where emitted energy = 1 MeV = 1.602176634e-13 J.
"""

# Import necessary libraries
import re

import numpy as np
import pandas as pd

import b_config.a_config as config

# ============================================================
# Constants
# ============================================================

MeV_to_J = config.MEV_TO_J

# ============================================================
# GEANT4 PHOTON-FLUENCE FILENAME PATTERN
# ============================================================

fluence_filename_pattern = re.compile(
    r"Geant4_deposit_MRCP_"
    r"(AM|AF)_source_"
    r"(.+?)_"
    r"(gamma|e-)_energy_"
    r"([0-9Ee.+-]+)"
    r"_photon_fluence\.csv$",
    re.IGNORECASE
)

# ============================================================
# ACTIVE MARROW MASSES (kg)
# ============================================================

marrow_mass_kg = config.MARROW_MASS_KG

# ============================================================
# FILES
# ============================================================

def find_fluence_files():

    root = config.GEANT4_GENERATED_INPUTS_DIR
    files = sorted(root.rglob("*photon_fluence*.csv"))

    return files

# ============================================================
# PARSE ICRP TABLE D.1
# ============================================================

def load_icrp_response_functions():

    icrp_file = config.SKELETAL_RESPONSE_FUNCTIONS_CSV

    if not icrp_file.exists():

        raise FileNotFoundError(
            f"Could not find ICRP 116 Table D.1:\n"
            f"{icrp_file}"
        )

    raw = pd.read_csv(
        icrp_file,
        header=None
    )

    # --------------------------------------------------------
    # Find all "Organ ID:" rows
    # --------------------------------------------------------

    organ_header_rows = []

    for row_index in range(len(raw)):

        for col_index, value in enumerate(
            raw.iloc[row_index]
        ):

            if pd.isna(value):
                continue

            match = re.search(
                r"Organ\s*ID\s*:\s*(\d+)",
                str(value),
                re.IGNORECASE
            )

            if match:

                organ_header_rows.append(
                    (
                        row_index,
                        col_index,
                        int(match.group(1))
                    )
                )

    if not organ_header_rows:

        raise RuntimeError(
            "No 'Organ ID:' entries were found "
            "in the ICRP 116 CSV."
        )

    # --------------------------------------------------------
    # Extract AM response functions
    # --------------------------------------------------------

    response_functions = {}

    for block_index, (
        header_row,
        _,
        _
    ) in enumerate(organ_header_rows):

        if block_index + 1 < len(
            organ_header_rows
        ):

            next_header_row = (
                organ_header_rows[
                    block_index + 1
                ][0]
            )

        else:

            next_header_row = len(raw)

        # ----------------------------------------------------
        # Find all organ IDs in this horizontal block
        # ----------------------------------------------------

        block_organs = []

        for col_index, value in enumerate(
            raw.iloc[header_row]
        ):

            if pd.isna(value):
                continue

            match = re.search(
                r"Organ\s*ID\s*:\s*(\d+)",
                str(value),
                re.IGNORECASE
            )

            if match:

                block_organs.append(
                    (
                        int(match.group(1)),
                        col_index
                    )
                )

        # ----------------------------------------------------
        # Find "(MeV)" row
        # ----------------------------------------------------

        data_start = header_row + 1

        for r in range(
            header_row,
            min(
                header_row + 10,
                next_header_row
            )
        ):

            first = raw.iloc[r, 0]

            if (
                pd.notna(first)
                and "(MeV)" in str(first)
            ):

                data_start = r + 1
                break

        # ----------------------------------------------------
        # Extract response for each organ
        # ----------------------------------------------------

        for organ_id, response_column in block_organs:

            rows = []

            for r in range(
                data_start,
                next_header_row
            ):

                energy_value = raw.iloc[r, 0]

                try:

                    energy = float(
                        energy_value
                    )

                except (
                    ValueError,
                    TypeError
                ):

                    continue

                response_value = raw.iloc[
                    r,
                    response_column
                ]

                try:

                    response = float(
                        response_value
                    )

                except (
                    ValueError,
                    TypeError
                ):

                    response = np.nan

                rows.append(
                    (
                        energy,
                        response
                    )
                )

            response_df = pd.DataFrame(
                rows,
                columns=[
                    "Energy_MeV",
                    "AM_Gy_m2"
                ]
            )

            response_df = response_df.dropna(
                subset=[
                    "Energy_MeV",
                    "AM_Gy_m2"
                ]
            )

            response_df = response_df[
                (response_df["Energy_MeV"] > 0)
                &
                (response_df["AM_Gy_m2"] > 0)
            ]

            response_df = response_df.sort_values(
                "Energy_MeV"
            )

            response_functions[
                organ_id
            ] = response_df.reset_index(
                drop=True
            )

    return response_functions

# ============================================================
# LOG-LOG INTERPOLATION
# ============================================================

def interpolate_response(
    energies,
    response_df
):

    table_E = response_df[
        "Energy_MeV"
    ].to_numpy(
        dtype=float
    )

    table_R = response_df[
        "AM_Gy_m2"
    ].to_numpy(
        dtype=float
    )

    if len(table_E) < 2:

        raise RuntimeError(
            "Fewer than two valid ICRP "
            "response points."
        )

    E = np.asarray(
        energies,
        dtype=float
    )

    if np.any(E <= 0):

        raise ValueError(
            "Requested energy must be > 0."
        )

    if np.any(
        E < table_E.min()
    ):

        raise ValueError(
            f"Energy below ICRP range: "
            f"{table_E.min()} MeV"
        )

    if np.any(
        E > table_E.max()
    ):

        raise ValueError(
            f"Energy above ICRP range: "
            f"{table_E.max()} MeV"
        )

    return np.exp(
        np.interp(
            np.log(E),
            np.log(table_E),
            np.log(table_R)
        )
    )


# ============================================================
# CALCULATE ONE FLUENCE FILE
# ============================================================

def calculate_from_fluence(
    fluence_file,
    response_functions,
    source_energy,
    source_organ,
    source_type,
    phantom_code
):

    print()
    print("=" * 90)
    print(
        f"Processing fluence file:\n"
        f"{fluence_file.name}"
    )
    print("=" * 90)

    fluence = pd.read_csv(
        fluence_file
    )

    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    required_columns = [

        "Spongiosa ID",

        "Energy Low (MeV)",
        "Energy High (MeV)",
        "Energy Center (MeV)",

        "Fluence (photons/m2/source)",

    ]

    missing = [

        column

        for column
        in required_columns

        if column
        not in fluence.columns

    ]

    if missing:

        raise ValueError(
            f"Missing columns in "
            f"{fluence_file.name}:\n"
            + "\n".join(
                f"  {column}"
                for column in missing
            )
        )

    # --------------------------------------------------------
    # Convert data types
    # --------------------------------------------------------

    for column in required_columns:

        fluence[column] = pd.to_numeric(
            fluence[column],
            errors="coerce"
        )

    fluence = fluence.dropna(
        subset=[
            "Spongiosa ID",
            "Energy Center (MeV)"
        ]
    )

    fluence["Spongiosa ID"] = (
        fluence["Spongiosa ID"]
        .astype(int)
    )

    # --------------------------------------------------------
    # Check all ICRP response functions
    # --------------------------------------------------------

    missing_icrp = [

        organ_id

        for organ_id in marrow_mass_kg

        if organ_id
        not in response_functions

    ]

    if missing_icrp:

        raise RuntimeError(
            "Missing ICRP response functions "
            f"for IDs: {missing_icrp}"
        )

    # --------------------------------------------------------
    # Calculate each skeletal site
    # --------------------------------------------------------

    site_results = []

    for organ_id in marrow_mass_kg:

        site = fluence[
            fluence["Spongiosa ID"] == organ_id
        ].copy()

        if site.empty:

            raise RuntimeError(
                f"No Geant4 fluence data "
                f"for spongiosa ID "
                f"{organ_id}."
            )

        response_df = (
            response_functions[
                organ_id
            ]
        )

        icrp_min_E = (
            response_df[
                "Energy_MeV"
            ].min()
        )

        icrp_max_E = (
            response_df[
                "Energy_MeV"
            ].max()
        )

        # ----------------------------------------------------
        # Total fluence
        # ----------------------------------------------------

        total_fluence = (
            site[
                "Fluence (photons/m2/source)"
            ]
            .sum()
        )

        # ----------------------------------------------------
        # Keep bins inside ICRP range
        # ----------------------------------------------------

        in_range = (

            site[
                "Energy Center (MeV)"
            ]
            >= icrp_min_E

        ) & (

            site[
                "Energy Center (MeV)"
            ]
            <= icrp_max_E

        )

        covered = site[
            in_range
        ].copy()

        excluded_fluence = (
            site.loc[
                ~in_range,
                "Fluence (photons/m2/source)"
            ]
            .sum()
        )

        if covered.empty:

            raise RuntimeError(
                f"No usable energy bins "
                f"for spongiosa ID "
                f"{organ_id}."
            )

        # ----------------------------------------------------
        # Dose calculation
        #
        # Phi is already bin-integrated.
        # Therefore:
        #
        # D = sum(Phi * R)
        #
        # NO Delta-E.
        # ----------------------------------------------------

        E = covered[
            "Energy Center (MeV)"
        ].to_numpy(
            dtype=float
        )

        Phi = covered[
            "Fluence (photons/m2/source)"
        ].to_numpy(
            dtype=float
        )

        R = interpolate_response(
            E,
            response_df
        )

        dose_contribution = (
            Phi * R
        )

        site_dose = np.sum(
            dose_contribution
        )

        site_mass_kg = marrow_mass_kg[organ_id]

        site_results.append({

            "Organ ID":
                organ_id,

            "AM mass (kg)":
                site_mass_kg,

            "Total fluence (photons/m2/source)":
                total_fluence,

            "Excluded fluence (photons/m2/source)":
                excluded_fluence,

            "Excluded fraction":
                (
                    excluded_fluence
                    / total_fluence
                    if total_fluence > 0
                    else 0.0
                ),

            "AM dose (Gy/source)":
                site_dose,

        })

    results = pd.DataFrame(
        site_results
    )

    # --------------------------------------------------------
    # Mass weighting
    # --------------------------------------------------------

    total_mass_kg = results["AM mass (kg)"].sum()

    results["AM mass fraction"] = (results["AM mass (kg)"] / total_mass_kg)

    results["Mass-weighted dose contribution (Gy/source)"] = (
        results["AM dose (Gy/source)"]

        *

        results["AM mass fraction"]
    )

    # --------------------------------------------------------
    # Total RBM dose
    # --------------------------------------------------------

    rbm_dose = (
        results[
            "Mass-weighted dose contribution (Gy/source)"
        ].sum()
    )

    # --------------------------------------------------------
    # SAF
    # --------------------------------------------------------

    emitted_energy_J = (source_energy * MeV_to_J)

    rbm_saf = (
        rbm_dose
        / emitted_energy_J
    )

    # --------------------------------------------------------
    # Add source information to every row
    # -------------------------------------------------------- 

    results.insert(
        0,
        "Phantom",
        config.PHANTOM_NAMES[phantom_code]
    )

    results.insert(
        1,
        "Source Organ",
        source_organ
    )

    results.insert(
        2,
        "Source Type",
        source_type
    )

    results.insert(
        3,
        "Source Energy (MeV)",
        source_energy
    )

    results["Total Active Marrow Mass (kg)"] = (
        total_mass_kg
    )

    results["RBM Dose (Gy/source)"] = (
        rbm_dose
    )

    results["RBM SAF (kg^-1)"] = (
        rbm_saf
    )

    # --------------------------------------------------------
    # Terminal output
    # --------------------------------------------------------

    print()

    print(
        results[
            [
                "Organ ID",
                "AM mass (kg)",
                "Total fluence (photons/m2/source)",
                "Excluded fraction",
                "AM dose (Gy/source)",
                "AM mass fraction",
            ]
        ].to_string(
            index=False,
            float_format=lambda x:
                f"{x:.6e}"
        )
    )

    print()
    print("-" * 90)

    print(
        f"Total active marrow mass : "
        f"{total_mass_kg:.1f} kg"
    )

    print(
        f"RBM dose                 : "
        f"{rbm_dose:.6e} Gy/source"
    )

    print(
        f"Source energy            : "
        f"{source_energy:.6g} MeV"
    )

    print(
        f"RBM SAF                  : "
        f"{rbm_saf:.6e} kg^-1"
    )

    print()
    print(
        f"RBM calculation completed for "
        f"{phantom_code}, {source_organ}, "
        f"{source_type}, {source_energy:g} MeV."
    )

    return results, rbm_dose, rbm_saf


# ============================================================
# MAIN PIPELINE FUNCTION
# ============================================================

def calculate_marrow_endosteum_SAFs(params):

    simulation = params[
        "simulation_code"
    ].upper()

    if simulation != "GEANT4":

        print(
            "\nSkipping ICRP 116 marrow "
            "fluence calculation."
        )

        return

    print()
    print("=" * 90)
    print(
        "ICRP 116 ACTIVE RED BONE MARROW "
        "CALCULATION"
    )
    print("=" * 90)

    # --------------------------------------------------------
    # Load ICRP response functions
    # --------------------------------------------------------

    response_functions = (
        load_icrp_response_functions()
    )

    # --------------------------------------------------------
    # Find Geant4 fluence files
    # --------------------------------------------------------

    fluence_files = find_fluence_files()

    if not fluence_files:

        raise FileNotFoundError(
            "No Geant4 photon-fluence CSV files "
            "were found in:\n"
            f"{config.GEANT4_BUILD_DIR}"
        )

    print(
        f"\nFound {len(fluence_files)} "
        "photon-fluence file(s)."
    )

    # --------------------------------------------------------
    # Calculate each source energy
    # --------------------------------------------------------

    results = []

    for fluence_file in fluence_files:

        match = fluence_filename_pattern.match(
            fluence_file.name
        )

        if not match:
            print(
                f"\n[WARNING] Could not parse "
                f"fluence filename:\n"
                f"  {fluence_file.name}"
            )
            continue

        phantom_code = match.group(1)

        source_organ = match.group(2)

        geant4_source_type = (
            match.group(3).lower()
        )

        source_energy = float(
            match.group(4)
        )

        source_type = (
            config.GEANT4_SOURCE_TYPE_MAP.get(
                geant4_source_type,
                geant4_source_type
            )
        )

        try:

            result = calculate_from_fluence(
                fluence_file,
                response_functions,
                source_energy,
                source_organ,
                source_type,
                phantom_code
            )

            results.append(result)

        except Exception as error:

            print(
                f"\n[WARNING] Could not process "
                f"{fluence_file.name}:"
            )

            print(
                f"  {error}"
            )

    # --------------------------------------------------------
    # Calculate every fluence file
    # --------------------------------------------------------

    all_results = []

    for fluence_file in fluence_files:

        match = fluence_filename_pattern.match(
            fluence_file.name
        )

        if not match:

            print(
                f"\n[WARNING] Could not parse "
                f"fluence filename:\n"
                f"  {fluence_file.name}"
            )

            continue

        phantom_code = match.group(1)

        source_organ = match.group(2)

        geant4_source_type = (
            match.group(3).lower()
        )

        source_energy = float(
            match.group(4)
        )

        source_type = (
            config.GEANT4_SOURCE_TYPE_MAP.get(
                geant4_source_type,
                geant4_source_type
            )
        )

        try:

            result, rbm_dose, rbm_saf = (
                calculate_from_fluence(
                    fluence_file,
                    response_functions,
                    source_energy,
                    source_organ,
                    source_type,
                    phantom_code
                )
            )

            all_results.append(result)

        except Exception as error:

            print(
                f"\n[WARNING] Could not process "
                f"{fluence_file.name}:"
            )

            print(
                f"  {error}"
            )

    if not all_results:

        raise RuntimeError(
            "No RBM results were successfully calculated."
        )

    combined_results = pd.concat(
        all_results,
        ignore_index=True
    )

    combined_results.sort_values(
        by=[
            "Phantom",
            "Source Organ",
            "Source Type",
            "Source Energy (MeV)",
            "Organ ID",
        ],
        inplace=True,
        ignore_index=True
    )

    output_file = (
        config.RESULTS_DIR /
        "o_geant4_RBM_ICRP116.csv"
    )

    combined_results.to_csv(
        output_file,
        index=False
    )

    print()
    print("=" * 90)
    print(
        "ICRP 116 RBM calculation complete."
    )
    print("=" * 90)

    print(
        f"Total result rows : "
        f"{len(combined_results)}"
    )

    print(
        f"Results saved to:\n"
        f"{output_file}"
    )

    return results