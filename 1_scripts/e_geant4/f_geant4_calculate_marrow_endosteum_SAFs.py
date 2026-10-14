"""
Calculate active red bone marrow (RBM) and endosteum dose and SAF using:

1. Geant4 energy-binned photon fluence in skeletal tissues.
2. phantom-specific photon dose-response functions (DRFs).

Geant4 CSV columns:
    Skeletal ID
    Energy Low (MeV)
    Energy High (MeV)
    Energy Center (MeV)
    Fluence (photons/m2/source)

Skeletal DRF CSV:
    Photon energy (MeV)
    RBM response function
    endosteum response function

The Geant4 fluence values are histogram/bin fluences,
so they are NOT multiplied by Delta E again.

ICRP response functions are interpolated in log(E)-log(R) space.

The final RBM and endosteum doses are mass-weighted
over the skeletal regions with valid
ICRP response functions.

For a photon source:

    SAF = dose / emitted energy

where:

    emitted energy = source energy (MeV) * 1.6021766339999e-13 J/MeV
"""

# ============================================================
# IMPORTS
# ============================================================

import re

import numpy as np
import pandas as pd

import b_config.a_config as config
from b_config.b_phantom_registry import (
    get_phantom,
    get_phantom_group,
    get_skeletal_ids,
    get_skeletal_mass_file,
    get_skeletal_response_functions_file,
)

# ============================================================
# CONSTANTS
# ============================================================

MeV_to_J = config.MEV_TO_J

# ============================================================
# Configs
# ============================================================

geant4_generated_inputs_dir = config.GEANT4_GENERATED_INPUTS_DIR
# The actual DRF file is selected from the phantom registry.
geant4_source_type_map = config.GEANT4_SOURCE_TYPE_MAP
geant4_results_dir = config.RESULTS_GEANT4_DIR
geant4_output_file = "d_geant4_rbm_endosteum_icrp116.csv"

# ============================================================
# GEANT4 PHOTON-FLUENCE FILENAME PATTERN
# ============================================================

fluence_filename_pattern = re.compile(
    r"geant4_deposit_(.+?)_source_"
    r"(.+?)_"
    r"(gamma|e-)_energy_"
    r"([0-9Ee.+-]+)"
    r"_photon_fluence\.csv$",
    re.IGNORECASE
)

# ============================================================
# LOAD SKELETAL TISSUE MASSES
# ============================================================

def load_skeletal_masses(mass_file):

    if not mass_file.exists():
        raise FileNotFoundError(f"Could not find skeletal tissue mass file:\n {mass_file}")

    masses = pd.read_csv(mass_file)

    required_columns = ["Organ ID",
                        "Ref_AM_Marrow_Mass(g)",
                        "Ref_AM_Endosteum_Mass(g)",
                        "Ref_AF_Marrow_Mass(g)",
                        "Ref_AF_Endosteum_Mass(g)",]

    missing = [
        column
        for column
        in required_columns
        if column
        not in masses.columns
    ]

    if missing:
        raise ValueError(
            "Missing columns in skeletal tissue mass CSV:\n"
            + "\n".join(
                f"  {column}"
                for column
                in missing
            )
        )

    masses["Organ ID"] = pd.to_numeric(masses["Organ ID"], errors="coerce")
    masses = masses.dropna(subset=["Organ ID"])
    masses["Organ ID"] = (masses["Organ ID"].astype(int))
    masses = masses.set_index("Organ ID")  

    return masses

# ============================================================
# FIND GEANT4 FLUENCE FILES
# ============================================================

def find_fluence_files():

    root = geant4_generated_inputs_dir

    return sorted(root.rglob("*_photon_fluence.csv"))

# ============================================================
# LOAD SKELETAL RESPONSE FUNCTIONS
# ============================================================

def load_response_functions(response_file):

    if not response_file.exists():

        raise FileNotFoundError(
            f"Could not find skeletal photon dose-response function file:\n"
            f"{response_file}"
        )

    raw = pd.read_csv(
        response_file,
        header=None,
        dtype=object
    )

    response_functions = {}

    # --------------------------------------------------------
    # Supported CSV block formats
    #
    # Each skeletal-site block is expected to contain:
    #
    #   row     = Organ ID: <ID>
    #   row + 1 = column headers
    #   row + 2 onward = energy / RBM / endosteum data
    #
    # The loader accepts either:
    #
    #   ICRP-116: AM / TM50
    #   Korean/Asian: RBM / BE
    #
    # and standardizes them internally to:
    #
    #   RBM_Gy_m2
    #   Endosteum_Gy_m2
    #
    # Organ ID may appear in any column of the Organ ID row.
    # --------------------------------------------------------

    def parse_float(value):
        if pd.isna(value):
            return np.nan

        try:
            return float(
                str(value).strip().replace("E", "e")
            )
        except (ValueError, TypeError):
            return np.nan

    def find_column(headers, keywords):
        for index, header in enumerate(headers):
            header = str(header).strip().lower()

            if any(keyword in header for keyword in keywords):
                return index

        return None

    for row in range(raw.shape[0]):

        # ----------------------------------------------------
        # Find "Organ ID: ####" anywhere in this row.
        # ----------------------------------------------------

        organ_id = None

        for column in range(raw.shape[1]):

            value = raw.iat[row, column]

            if pd.isna(value):
                continue

            match = re.search(
                r"Organ\s*ID\s*:\s*(\d+)",
                str(value),
                re.IGNORECASE
            )

            if match is not None:
                organ_id = int(match.group(1))
                break

        if organ_id is None:
            continue

        header_row = row + 1

        if header_row >= raw.shape[0]:
            raise RuntimeError(
                f"Missing response-function header after Organ ID {organ_id} "
                f"in {response_file.name}."
            )

        headers = [
            str(raw.iat[header_row, column]).strip()
            if not pd.isna(raw.iat[header_row, column])
            else ""
            for column in range(raw.shape[1])
        ]

        energy_column = find_column(
            headers,
            ("photon energy", "energy")
        )

        rbm_column = find_column(
            headers,
            ("rbm", "am response", "am")
        )

        endosteum_column = find_column(
            headers,
            ("tm50", "be", "endosteum")
        )

        if energy_column is None:
            energy_column = 0

        if rbm_column is None or endosteum_column is None:
            raise RuntimeError(
                f"Could not identify RBM/endosteum response columns "
                f"for Organ ID {organ_id} in {response_file.name}.\\n"
                f"Headers found: {headers}"
            )

        # ----------------------------------------------------
        # Read the data rows until the next block.
        # ----------------------------------------------------

        energies = []
        rbm_values = []
        endosteum_values = []

        current_row = row + 2

        while current_row < raw.shape[0]:

            # A new Organ ID marks the next block.
            row_contains_organ_id = False

            for column in range(raw.shape[1]):
                value = raw.iat[current_row, column]

                if pd.isna(value):
                    continue

                if re.search(
                    r"Organ\s*ID\s*:\s*\d+",
                    str(value),
                    re.IGNORECASE
                ):
                    row_contains_organ_id = True
                    break

            if row_contains_organ_id:
                break

            energy = parse_float(
                raw.iat[current_row, energy_column]
            )

            # Stop at non-data rows.
            if not np.isfinite(energy):
                break

            rbm = parse_float(
                raw.iat[current_row, rbm_column]
            )

            endosteum = parse_float(
                raw.iat[current_row, endosteum_column]
            )

            energies.append(energy)
            rbm_values.append(rbm)
            endosteum_values.append(endosteum)

            current_row += 1

        response_df = pd.DataFrame({
            "Energy_MeV": energies,
            "RBM_Gy_m2": rbm_values,
            "Endosteum_Gy_m2": endosteum_values,
        })

        response_df = response_df[
            response_df["RBM_Gy_m2"].notna()
            | response_df["Endosteum_Gy_m2"].notna()
        ].reset_index(drop=True)

        if not response_df.empty:
            response_functions[organ_id] = response_df

    if not response_functions:
        raise RuntimeError(
            f"No skeletal response-function blocks were found in:\\n"
            f"{response_file}"
        )

    return response_functions


# ============================================================
# LOG-LOG INTERPOLATION
# ============================================================

def interpolate_response(
    energies,
    response_df,
    response_column
):

    response_data = response_df[
        ["Energy_MeV",
         response_column]
         ].copy()

    response_data = response_data.dropna()
    response_data = response_data[
        (response_data["Energy_MeV"] > 0)
        &
        (response_data[response_column] > 0)]

    response_data = response_data.sort_values("Energy_MeV")

    # No response function for this target
    if response_data.empty:
        return None
    
    table_E = response_data[
        "Energy_MeV"].to_numpy(dtype=float)

    table_R = response_data[
        response_column].to_numpy(dtype=float)

    if len(table_E) < 2:
        raise RuntimeError(
            f"Fewer than two valid ICRP "
            f"response points for "
            f"{response_column}.")

    E = np.asarray(energies, dtype=float)

    if np.any(E <= 0):
        raise ValueError("Requested energy must be > 0.")

    if np.any(E < table_E.min()):
        raise ValueError(f"Energy below ICRP range: {table_E.min()} MeV")

    if np.any(E > table_E.max()):
        raise ValueError(f"Energy above ICRP range: {table_E.max()} MeV")
    
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
    skeletal_masses,
    rbm_ids,
    endosteum_ids,
    source_energy,
    source_organ,
    source_type,
    phantom_code,
    params
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

    # ========================================================
    # REQUIRED FLUENCE COLUMNS
    # ========================================================

    required_columns = [
        "Skeletal ID",
        "Energy Low (MeV)",
        "Energy High (MeV)",
        "Energy Center (MeV)",
        "Fluence (photons/m2/source)",
        "Relative Error",
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
            f"Missing columns in {fluence_file.name}:\n"
            + "\n".join(
                f"  {column}"
                for column
                in missing
            )
        )

    # ========================================================
    # CONVERT DATA TYPES
    # ========================================================

    for column in required_columns:

        fluence[column] = pd.to_numeric(fluence[column], errors="coerce")

    fluence = fluence.dropna(
        subset=["Skeletal ID",
                "Energy Center (MeV)"]
    )

    fluence["Skeletal ID"] = (fluence["Skeletal ID"].astype(int))

    # ========================================================
    # SELECT MASS COLUMNS
    # ========================================================

    phantom_spec = get_phantom(phantom_code)

    if phantom_spec.skeletal_mass_columns is None:
        raise ValueError(f"No skeletal mass columns defined for phantom '{phantom_code}'.")

    marrow_column, endosteum_column = phantom_spec.skeletal_mass_columns

    # ========================================================
    # Skeletal IDs from ICRP 116 Table D.1
    # ========================================================

    skeletal_ids = (rbm_ids | endosteum_ids)

    # ========================================================
    # CHECK THAT ALL MASSES EXIST
    # ========================================================

    missing_mass_ids = [
        organ_id
        for organ_id
        in skeletal_ids
        if organ_id
        not in skeletal_masses.index
    ]

    if missing_mass_ids:
        raise RuntimeError(f"Missing skeletal tissue masses for IDs: {missing_mass_ids}")

    # ========================================================
    # CALCULATE EACH SKELETAL SITE
    # ========================================================

    site_results = []

    for organ_id in sorted(skeletal_ids):

        # ----------------------------------------------------
        # Fluence for this skeletal site
        # ----------------------------------------------------

        site = fluence[
            fluence["Skeletal ID"]
            == organ_id
        ].copy()

        if site.empty:
            raise RuntimeError(
                f"No Geant4 fluence data "
                f"for skeletal ID "
                f"{organ_id}."
            )

        # ----------------------------------------------------
        # ICRP response function
        # ----------------------------------------------------

        if organ_id not in response_functions:

            raise RuntimeError(
                f"No ICRP response function "
                f"for skeletal ID "
                f"{organ_id}."
            )

        response_df = (response_functions[organ_id])

        # ----------------------------------------------------
        # ICRP energy range
        # ----------------------------------------------------

        icrp_min_E = (response_df[
                        "Energy_MeV"]
                        .min())

        icrp_max_E = (response_df[
                        "Energy_MeV"]
                        .max())

        # ----------------------------------------------------
        # Total fluence
        # ----------------------------------------------------

        total_fluence = (
            site["Fluence (photons/m2/source)"]
            .sum()
        )

        # ----------------------------------------------------
        # RBM energy range
        # ----------------------------------------------------

        rbm_table_E = response_df[
            response_df["RBM_Gy_m2"].notna()
        ]["Energy_MeV"]

        if rbm_table_E.empty:
            marrow_dose = np.nan
            marrow_relative_error = np.nan
            excluded_marrow_fluence = np.nan
            marrow_mass_kg = np.nan
        else:

            rbm_in_range = (
                (site["Energy Center (MeV)"] >= rbm_table_E.min())
                &
                (site["Energy Center (MeV)"] <= rbm_table_E.max())
            )

            rbm_covered = site[rbm_in_range].copy()

            if rbm_covered.empty:
                raise RuntimeError(
                    f"No usable Geant4 energy bins inside ICRP AM range "
                    f"for skeletal ID {organ_id}."
                )

            E = rbm_covered[
                "Energy Center (MeV)"
            ].to_numpy(dtype=float)

            Phi = rbm_covered[
                "Fluence (photons/m2/source)"
            ].to_numpy(dtype=float)

            relative_error = rbm_covered[
                "Relative Error"
            ].to_numpy(dtype=float)

            R_AM = interpolate_response(
                E,
                response_df,
                "RBM_Gy_m2"
            )

            dose_contribution = Phi * R_AM

            marrow_dose = np.sum(
                dose_contribution
            )

            marrow_sigma = np.sqrt(
                np.sum(
                    (
                        dose_contribution
                        * relative_error
                    ) ** 2
                )
            )

            marrow_relative_error = (
                marrow_sigma / marrow_dose
                if marrow_dose > 0
                else np.nan
            )

            excluded_marrow_fluence = (
                site.loc[
                    ~rbm_in_range,
                    "Fluence (photons/m2/source)"
                ].sum()
            )

            marrow_mass_kg = (
                float(
                    skeletal_masses.loc[
                        organ_id,
                        marrow_column
                    ]
                )
                / 1000.0
            )

        # ----------------------------------------------------
        # Endosteum energy range
        # ----------------------------------------------------

        endo_table_E = response_df[
            response_df["Endosteum_Gy_m2"].notna()
        ]["Energy_MeV"]

        if endo_table_E.empty:
            endosteum_dose = np.nan
            endosteum_relative_error = np.nan
            excluded_endosteum_fluence = np.nan
            endosteum_mass_kg = np.nan
        else:

            endo_in_range = (
                (site["Energy Center (MeV)"] >= endo_table_E.min())
                &
                (site["Energy Center (MeV)"] <= endo_table_E.max())
            )

            endo_covered = site[endo_in_range].copy()

            if endo_covered.empty:
                raise RuntimeError(
                    f"No usable Geant4 energy bins inside ICRP TM50 range "
                    f"for skeletal ID {organ_id}."
                )

            E = endo_covered[
                "Energy Center (MeV)"
            ].to_numpy(dtype=float)

            Phi = endo_covered[
                "Fluence (photons/m2/source)"
            ].to_numpy(dtype=float)

            relative_error = endo_covered[
                "Relative Error"
            ].to_numpy(dtype=float)

            R_TM50 = interpolate_response(
                E,
                response_df,
                "Endosteum_Gy_m2"
            )

            dose_contribution = Phi * R_TM50

            endosteum_dose = np.sum(
                dose_contribution
            )

            endosteum_sigma = np.sqrt(
                np.sum(
                    (
                        dose_contribution
                        * relative_error
                    ) ** 2
                )
            )

            endosteum_relative_error = (
                endosteum_sigma / endosteum_dose
                if endosteum_dose > 0
                else np.nan
            )

            excluded_endosteum_fluence = (
                site.loc[
                    ~endo_in_range,
                    "Fluence (photons/m2/source)"
                ].sum()
            )

            endosteum_mass_kg = (
                float(
                    skeletal_masses.loc[
                        organ_id,
                        endosteum_column
                    ]
                )
                / 1000.0
            )

        # ====================================================
        # STORE SITE RESULT
        # ====================================================

        site_results.append({

            "Organ ID":
                organ_id,

            "Marrow mass (kg)":
                marrow_mass_kg,

            "Endosteum mass (kg)":
                endosteum_mass_kg,

            "Total fluence (photons/m2/source)":
                total_fluence / 1.0e4,

            "Marrow excluded fluence (photons/m2/source)":
                excluded_marrow_fluence / 1.0e4,

            "Endosteum excluded fluence (photons/m2/source)":
                excluded_endosteum_fluence / 1.0e4,

            "Marrow dose (Gy/source)":
                marrow_dose,

            "Marrow Relative Error":
                marrow_relative_error,

            "Marrow Statistical Uncertainty (%)":
                (
                    marrow_relative_error * 100
                    if np.isfinite(marrow_relative_error)
                    else np.nan
                ),

            "Endosteum dose (Gy/source)":
                endosteum_dose,

            "Endosteum Relative Error":
                endosteum_relative_error,

            "Endosteum Statistical Uncertainty (%)":
                (
                    endosteum_relative_error * 100
                    if np.isfinite(endosteum_relative_error)
                    else np.nan
                ),
        })

    # ========================================================
    # RESULTS DATAFRAME
    # ========================================================

    results = pd.DataFrame(site_results)

    # ========================================================
    # RBM MASS WEIGHTING
    # ========================================================

    rbm_results = results[
        results["Marrow dose (Gy/source)"].notna()
    ].copy()

    rbm_results = rbm_results[
        rbm_results["Marrow mass (kg)"].notna()
    ].copy()

    if rbm_results.empty:

        total_marrow_mass_kg = np.nan
        rbm_dose = np.nan
        rbm_relative_error = np.nan
        rbm_statistical_uncertainty = np.nan

    else:

        total_marrow_mass_kg = (
            rbm_results["Marrow mass (kg)"].sum()
        )

        rbm_results["Marrow mass fraction"] = (
            rbm_results["Marrow mass (kg)"]
            / total_marrow_mass_kg
        )

        rbm_results[
            "Mass-weighted marrow dose contribution (Gy/source)"
        ] = (
            rbm_results["Marrow dose (Gy/source)"]
            * rbm_results["Marrow mass fraction"]
        )

        rbm_dose = (
            rbm_results[
                "Mass-weighted marrow dose contribution (Gy/source)"
            ].sum()
        )

        rbm_sigma = np.sqrt(
            np.sum(
                (
                    rbm_results["Marrow mass fraction"]
                    * rbm_results["Marrow dose (Gy/source)"]
                    * rbm_results["Marrow Relative Error"]
                ) ** 2
            )
        )

        rbm_relative_error = (
            rbm_sigma / rbm_dose
            if rbm_dose > 0
            else np.nan
        )

        rbm_statistical_uncertainty = (
            rbm_relative_error * 100
            if np.isfinite(rbm_relative_error)
            else np.nan
        )

    # ========================================================
    # ENDOSTEUM MASS WEIGHTING
    # ========================================================

    endosteum_results = results[
        results["Endosteum dose (Gy/source)"].notna()
    ].copy()

    endosteum_results = endosteum_results[
        endosteum_results["Endosteum mass (kg)"].notna()
    ].copy()

    if endosteum_results.empty:

        total_endosteum_mass_kg = np.nan
        endosteum_dose = np.nan
        endosteum_relative_error = np.nan
        endosteum_statistical_uncertainty = np.nan

    else:

        total_endosteum_mass_kg = (
            endosteum_results["Endosteum mass (kg)"].sum()
        )

        endosteum_results[
            "Endosteum mass fraction"
        ] = (
            endosteum_results["Endosteum mass (kg)"]
            / total_endosteum_mass_kg
        )

        endosteum_results[
            "Mass-weighted endosteum dose contribution (Gy/source)"
        ] = (
            endosteum_results["Endosteum dose (Gy/source)"]
            * endosteum_results["Endosteum mass fraction"]
        )

        endosteum_dose = (
            endosteum_results[
                "Mass-weighted endosteum dose contribution (Gy/source)"
            ].sum()
        )

        endosteum_sigma = np.sqrt(
            np.sum(
                (
                    endosteum_results["Endosteum mass fraction"]
                    * endosteum_results["Endosteum dose (Gy/source)"]
                    * endosteum_results["Endosteum Relative Error"]
                ) ** 2
            )
        )

        endosteum_relative_error = (
            endosteum_sigma / endosteum_dose
            if endosteum_dose > 0
            else np.nan
        )

        endosteum_statistical_uncertainty = (
            endosteum_relative_error * 100
            if np.isfinite(endosteum_relative_error)
            else np.nan
        )

    # ========================================================
    # PUT MASS FRACTIONS BACK INTO RESULTS
    # ========================================================

    results["Marrow mass fraction"] = np.nan
    results["Endosteum mass fraction"] = np.nan

    results.loc[
        rbm_results.index,
        "Marrow mass fraction"
    ] = rbm_results["Marrow mass fraction"]

    results.loc[
        endosteum_results.index,
        "Endosteum mass fraction"
    ] = endosteum_results["Endosteum mass fraction"]

    # ========================================================
    # SAF
    # ========================================================

    emitted_energy_J = (
        source_energy
        * MeV_to_J
    )

    rbm_saf = (
        rbm_dose
        / emitted_energy_J
    )

    endosteum_saf = (
        endosteum_dose
        / emitted_energy_J
    )

    # ========================================================
    # ADD SOURCE INFORMATION
    # ========================================================
    number_of_particles = params["nps"]

    results.insert(0, "Phantom", phantom_spec.display_name)
    results.insert(1, "Source Organ", source_organ)
    results.insert(2, "Source Type", source_type)
    results.insert(3, "Source Energy (MeV)",source_energy)
    results.insert(4, "Number of Particles", number_of_particles)

    # ========================================================
    # ADD TOTALS TO EVERY ROW
    # ========================================================

    results[
        "Total Active Marrow Mass (kg)"
    ] = (
        total_marrow_mass_kg
    )

    results[
        "Total Endosteum Mass (kg)"
    ] = (
        total_endosteum_mass_kg
    )

    results[
        "RBM Dose (Gy/source)"
    ] = (
        rbm_dose
    )

    results[
        "RBM Relative Error"
    ] = rbm_relative_error

    results[
        "RBM Statistical Uncertainty (%)"
    ] = rbm_statistical_uncertainty

    results[
        "RBM SAF (kg^-1)"
    ] = rbm_saf

    results[
        "RBM SAF Uncertainty (kg^-1)"
    ] = (
        rbm_relative_error * rbm_saf
        if np.isfinite(rbm_relative_error)
        else np.nan
    )

    results[
        "Endosteum Dose (Gy/source)"
    ] = (
        endosteum_dose
    )

    results[
        "Endosteum Relative Error"
    ] = endosteum_relative_error

    results[
        "Endosteum Statistical Uncertainty (%)"
    ] = endosteum_statistical_uncertainty

    results[
        "Endosteum SAF (kg^-1)"
    ] = endosteum_saf

    results[
        "Endosteum SAF Uncertainty (kg^-1)"
    ] = (
        endosteum_relative_error * endosteum_saf
        if np.isfinite(endosteum_relative_error)
        else np.nan
    )

    # ========================================================
    # TERMINAL OUTPUT
    # ========================================================

    print()
    print(
        results[
            [
                "Organ ID",
                "Marrow mass (kg)",
                "Endosteum mass (kg)",
                "Total fluence (photons/m2/source)",
                "Marrow dose (Gy/source)",
                "Endosteum dose (Gy/source)",
                "Marrow mass fraction",
                "Endosteum mass fraction",
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
        f"{total_marrow_mass_kg:.6f} kg"
    )
    print(
        f"Total endosteum mass     : "
        f"{total_endosteum_mass_kg:.6f} kg"
    )
    print(
        f"RBM dose                 : "
        f"{rbm_dose:.6e} Gy/source"
    )
    print(
        f"Endosteum dose           : "
        f"{endosteum_dose:.6e} Gy/source"
    )
    print(
        f"Source energy            : "
        f"{source_energy:.6g} MeV"
    )
    print(
        f"RBM SAF                  : "
        f"{rbm_saf:.6e} kg^-1"
    )
    print(
        f"Endosteum SAF            : "
        f"{endosteum_saf:.6e} kg^-1"
    )

    print()
    print(
        f"RBM and endosteum "
        f"calculation completed for "
        f"{phantom_code}, "
        f"{source_organ}, "
        f"{source_type}, "
        f"{source_energy:g} MeV."
    )

    return results

# ============================================================
# MAIN PIPELINE FUNCTION
# ============================================================

def geant4_calculate_marrow_endosteum_SAFs(params):

    # --------------------------------------------------------
    # Only calculate for photon sources
    # --------------------------------------------------------

    if (params["source_type"].lower() not in {"gamma", "photon"}):

        print(
            f"\nSkipping ICRP 116 skeletal "
            f"calculation for source type: "
            f"{params['source_type']}"
        )

        print("ICRP 116 fluence-to-dose "
            "response functions implemented "
            "here are for photon sources.")

        return None

    print()
    print("=" * 90)
    print(
        "ICRP 116 RBM AND ENDOSTEUM "
        "CALCULATION"
    )
    print("=" * 90)

    # ========================================================
    # LOAD DATA
    # ========================================================

    phantom_selection = params["phantom"]
    selected_phantoms = get_phantom_group(phantom_selection)

    # Load the skeletal tissue mass file required by each phantom.
    # ICRP 145 mesh -> ICRP 145 mass file;
    # Filipino / Filipino-sized -> Asian mass file.
    skeletal_mass_cache = {}
    phantom_mass_data = {}

    # --------------------------------------------------------
    # Load the correct DRF file for EACH phantom.
    #
    # ICRP 145 adult male/female:
    #     ICRP116-Table-D-1-Updated.csv
    #
    # Resized ICRP 145 + Filipino male:
    #     DRFs-Male-Korean-Asian-Updated.csv
    #
    # Resized ICRP 145 + Filipino female:
    #     DRFs-Female-Korean-Asian-Updated.csv
    #
    # This is important because a phantom group may contain
    # both male and female phantoms, but their DRF files are
    # not necessarily the same.
    # --------------------------------------------------------

    response_cache = {}
    phantom_response_data = {}

    for phantom_code in selected_phantoms:

        mass_file = get_skeletal_mass_file(phantom_code)
        mass_cache_key = str(mass_file.resolve())

        if mass_cache_key not in skeletal_mass_cache:
            skeletal_mass_cache[mass_cache_key] = load_skeletal_masses(
                mass_file
            )

        skeletal_masses = skeletal_mass_cache[mass_cache_key]
        phantom_mass_data[phantom_code] = skeletal_masses

        response_file = get_skeletal_response_functions_file(
            phantom_code
        )

        response_cache_key = str(response_file.resolve())

        if response_cache_key not in response_cache:

            response_cache[response_cache_key] = (
                load_response_functions(response_file)
            )

        response_functions = response_cache[response_cache_key]

        configured_skeletal_ids = set(
            get_skeletal_ids(phantom_code)
        )

        rbm_ids = {
            organ_id
            for organ_id in configured_skeletal_ids
            if (
                organ_id in response_functions
                and response_functions[organ_id][
                    "RBM_Gy_m2"
                ].notna().any()
            )
        }

        endosteum_ids = {
            organ_id
            for organ_id in configured_skeletal_ids
            if (
                organ_id in response_functions
                and response_functions[organ_id][
                    "Endosteum_Gy_m2"
                ].notna().any()
            )
        }

        if not rbm_ids:
            raise RuntimeError(
                f"No valid RBM skeletal IDs were found for "
                f"phantom '{phantom_code}' using:\\n"
                f"{response_file}"
            )

        if not endosteum_ids:
            raise RuntimeError(
                f"No valid endosteum skeletal IDs were found for "
                f"phantom '{phantom_code}' using:\\n"
                f"{response_file}"
            )

        phantom_response_data[phantom_code] = (
            response_functions,
            rbm_ids,
            endosteum_ids,
        )

        print()
        print(f"Phantom: {phantom_code}")
        print(f"  DRF file              : {response_file.name}")
        print(f"  Configured skeletal IDs: {len(configured_skeletal_ids)}")
        print(f"  RBM DRF regions       : {len(rbm_ids)}")
        print(f"  Endosteum DRF regions : {len(endosteum_ids)}")

    # ========================================================
    # FIND GEANT4 FLUENCE FILES
    # ========================================================

    selected_phantom_names = {
        phantom_code.casefold()
        for phantom_code in selected_phantoms
    }

    selected_energies = np.asarray(
        params["source_energies"],
        dtype=float,
    )

    fluence_files = []

    for fluence_file in find_fluence_files():

        match = fluence_filename_pattern.fullmatch(
            fluence_file.name
        )

        if match is None:
            continue

        file_phantom = match.group(1).casefold()
        file_energy = float(match.group(4))

        # Keep only the selected phantoms.
        if file_phantom not in selected_phantom_names:
            continue

        # Compare energies numerically, not as strings.
        if not np.any(
            np.isclose(
                file_energy,
                selected_energies,
                rtol=1e-9,
                atol=1e-12,
            )
        ):
            continue

        fluence_files.append(fluence_file)

    fluence_files.sort()

    if not fluence_files:
        print(f"\nFluence search directory: {geant4_generated_inputs_dir}")
        print(f"Selected phantoms: {selected_phantoms}")
        print(f"Selected energies: {selected_energies.tolist()}")
        print(
            "CSV files discovered:",
            len(find_fluence_files()),
        )

        raise FileNotFoundError(
            "No matching Geant4 photon-fluence CSV files were found.\n"
            f"Search directory: {geant4_generated_inputs_dir}\n"
            f"Selected phantoms: {selected_phantoms}\n"
            f"Selected energies: {selected_energies.tolist()}"
        )

    print(
        f"\nFound {len(fluence_files)} matching "
        "Geant4 photon-fluence CSV file(s)."
    )

    # ========================================================
    # CALCULATE EVERY FLUENCE FILE
    # ========================================================

    all_results = []

    for fluence_file in fluence_files:

        match = (
            fluence_filename_pattern.match(
                fluence_file.name
            )
        )

        if not match:

            print(
                f"\n[WARNING] Could not parse "
                f"fluence filename:\n"
                f"  {fluence_file.name}"
            )

            continue

        phantom_code = next(
            name
            for name in selected_phantoms
            if name.casefold() == match.group(1).casefold()
        )

        source_organ = match.group(2)

        geant4_source_type = match.group(3).lower()

        source_energy = float(match.group(4))

        # ----------------------------------------------------
        # Convert Geant4 particle syntax
        #
        # gamma -> photon
        # e-    -> electron
        # ----------------------------------------------------

        source_type = (
            geant4_source_type_map.get(
                geant4_source_type,
                geant4_source_type
            )
        )

        try:

            (
                response_functions,
                rbm_ids,
                endosteum_ids,
            ) = phantom_response_data[phantom_code]

            skeletal_masses = phantom_mass_data[phantom_code]

            result = calculate_from_fluence(
                fluence_file,
                response_functions,
                skeletal_masses,
                rbm_ids,
                endosteum_ids,
                source_energy,
                source_organ,
                source_type,
                phantom_code,
                params
            )

            all_results.append(
                result
            )

        except Exception as error:

            print(
                f"\n[WARNING] Could not process "
                f"{fluence_file.name}:"
            )

            print(
                f"  {error}"
            )

    # ========================================================
    # CHECK RESULTS
    # ========================================================

    if not all_results:

        raise RuntimeError(
            "No RBM/endosteum results "
            "were successfully calculated."
        )

    # ========================================================
    # COMBINE RESULTS
    # ========================================================

    combined_results = pd.concat(
        all_results,
        ignore_index=True
    )

    # ========================================================
    # SORT RESULTS
    # ========================================================

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


    # ========================================================
    # KEEP OVERALL RESULTS ONLY ON FIRST ROW
    # ========================================================

    summary_columns = [
        "Total Active Marrow Mass (kg)",
        "Total Endosteum Mass (kg)",

        "RBM Dose (Gy/source)",
        "RBM Relative Error",
        "RBM Statistical Uncertainty (%)",
        "RBM SAF (kg^-1)",
        "RBM SAF Uncertainty (kg^-1)",

        "Endosteum Dose (Gy/source)",
        "Endosteum Relative Error",
        "Endosteum Statistical Uncertainty (%)",
        "Endosteum SAF (kg^-1)",
        "Endosteum SAF Uncertainty (kg^-1)",
    ]

    simulation_columns = [
        "Phantom",
        "Source Organ",
        "Source Type",
        "Source Energy (MeV)",
    ]

    for _, group in combined_results.groupby(
        simulation_columns,
        sort=False
    ):

        # Keep the summary values only
        # on the first skeletal-region row.
        remaining_indices = group.index[1:]

        combined_results.loc[
            remaining_indices,
            summary_columns
        ] = np.nan

    # ========================================================
    # COLUMN ORDER
    # ========================================================

    column_order = [
        "Phantom",
        "Source Organ",
        "Source Type",
        "Source Energy (MeV)",
        "Number of Particles",
        "Organ ID",

        "Marrow mass (kg)",
        "Endosteum mass (kg)",

        "Total fluence (photons/m2/source)",

        "Marrow excluded fluence (photons/m2/source)",
        "Endosteum excluded fluence (photons/m2/source)",

        "Marrow dose (Gy/source)",
        "Marrow Relative Error",
        "Marrow Statistical Uncertainty (%)",

        "Endosteum dose (Gy/source)",
        "Endosteum Relative Error",
        "Endosteum Statistical Uncertainty (%)",

        "Marrow mass fraction",
        "Endosteum mass fraction",

        "Total Active Marrow Mass (kg)",
        "Total Endosteum Mass (kg)",

        "RBM Dose (Gy/source)",
        "RBM Relative Error",
        "RBM Statistical Uncertainty (%)",
        "RBM SAF (kg^-1)",
        "RBM SAF Uncertainty (kg^-1)",

        "Endosteum Dose (Gy/source)",
        "Endosteum SAF (kg^-1)",
        "Endosteum SAF Uncertainty (kg^-1)",
    ]

    combined_results = combined_results[column_order]

    # ========================================================
    # KEEP AND SAVE SELECTED PHANTOM FAMILY
    # ========================================================
    
    output_file = geant4_results_dir / geant4_output_file

    family_results = combined_results.copy()
    family_results.to_csv(output_file, index=False)

    # ========================================================
    # Final Print
    # ========================================================

    print()
    print("=" * 90)
    print(
        "ICRP 116 RBM AND ENDOSTEUM "
        "CALCULATION COMPLETE"
    )
    print("=" * 90)

    print(f"Result rows : {len(family_results)}")

    print(f"Results saved to:\n{output_file}")

    return family_results