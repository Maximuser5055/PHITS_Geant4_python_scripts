"""
Calculate active red bone marrow (RBM) and endosteum dose and SAF using:

1. Geant4 energy-binned photon fluence in skeletal tissues.
2. ICRP Publication 116, Table D.1 photon response functions.

Geant4 CSV columns:
    Skeletal ID
    Energy Low (MeV)
    Energy High (MeV)
    Energy Center (MeV)
    Fluence (photons/m2/source)

ICRP CSV:
    Photon energy (MeV)
    AM response function
    TM50 response function

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

# ============================================================
# CONSTANTS
# ============================================================

MeV_to_J = config.MEV_TO_J

# ============================================================
# GEANT4 PHOTON-FLUENCE FILENAME PATTERN
# ============================================================

fluence_filename_pattern = re.compile(
    r"geant4_deposit_MRCP_"
    r"(AM|AF)_source_"
    r"(.+?)_"
    r"(gamma|e-)_energy_"
    r"([0-9Ee.+-]+)"
    r"_photon_fluence\.csv$",
    re.IGNORECASE
)

# ============================================================
# LOAD SKELETAL TISSUE MASSES
# ============================================================

def load_skeletal_masses():

    mass_file = config.SKELETAL_MASSES_CSV

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

    root = config.GEANT4_GENERATED_INPUTS_DIR

    return sorted(root.rglob("*_gamma_energy_*_photon_fluence.csv"))

# ============================================================
# LOAD ICRP 116 RESPONSE FUNCTIONS
# ============================================================

def load_icrp_response_functions():

    icrp_file = config.SKELETAL_RESPONSE_FUNCTIONS_CSV

    if not icrp_file.exists():

        raise FileNotFoundError(f"Could not find ICRP 116 Table D.1:\n {icrp_file}")

    icrp_raw = pd.read_csv(icrp_file, header=None)

    response_functions = {}

    # --------------------------------------------------------
    # The ICRP CSV structure is:
    # Column 0 = Photon energy
    # Column 1 = AM
    # Column 2 = TM50
    # Organ ID is found in column 1.
    # --------------------------------------------------------

    for row in range(
        icrp_raw.shape[0]
    ):

        value = icrp_raw.iat[row, 1] 

        if pd.isna(value):
            continue

        match = re.search(
            r"Organ\s*ID\s*:\s*(\d+)",
            str(value),
            re.IGNORECASE
        )

        if match is None:
            continue

        organ_id = int(match.group(1))

        # ----------------------------------------------------
        # First data row
        # row     = Organ ID
        # row + 1 = Photon Energy / AM / TM50
        # row + 2 = first data point
        # ----------------------------------------------------

        data_start = row + 2

        energies = []
        am_values = []
        tm50_values = []

        current_row = data_start

        while (current_row < icrp_raw.shape[0]):

            energy_value = icrp_raw.iat[current_row,0]

            try:
                energy = float(energy_value)

            except (ValueError, TypeError):
                break

            # ------------------------------------------------
            # AM
            # ------------------------------------------------

            am_value = icrp_raw.iat[current_row, 1]

            try:
                am = float(
                    str(am_value)
                    .strip()
                    .replace("E", "e")
                )

            except (ValueError, TypeError):
                am = np.nan

            # ------------------------------------------------
            # TM50
            # ------------------------------------------------

            tm50_value = icrp_raw.iat[current_row, 2]

            try:
                tm50 = float(
                    str(tm50_value)
                    .strip()
                    .replace("E", "e")
                )

            except (ValueError, TypeError):
                tm50 = np.nan

            energies.append(energy)

            am_values.append(am)

            tm50_values.append(tm50)

            current_row += 1

        # ----------------------------------------------------
        # Store response function
        # ----------------------------------------------------

        response_df = pd.DataFrame({
            "Energy_MeV": energies,
            "AM_Gy_m2": am_values,
            "TM50_Gy_m2": tm50_values,
        })

        # Keep the region if it has at least one valid
        # response function.
        response_df = response_df[
            (response_df["AM_Gy_m2"].notna())
            |
            (response_df["TM50_Gy_m2"].notna())
        ].reset_index(drop=True)

        response_functions[organ_id] = response_df

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

    # ========================================================
    # REQUIRED FLUENCE COLUMNS
    # ========================================================

    required_columns = [
        "Skeletal ID",
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

    if phantom_code == "AM":

        marrow_column = ("Ref_AM_Marrow_Mass(g)")
        endosteum_column = ("Ref_AM_Endosteum_Mass(g)")

    elif phantom_code == "AF":

        marrow_column = ("Ref_AF_Marrow_Mass(g)")
        endosteum_column = ("Ref_AF_Endosteum_Mass(g)")

    else:
        raise ValueError(f"Unsupported phantom code: {phantom_code}")

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
                .sum())

        # ----------------------------------------------------
        # Keep bins inside ICRP range
        # ----------------------------------------------------

        in_range = (
            site["Energy Center (MeV)"]
            >= icrp_min_E
        ) & (
            site["Energy Center (MeV)"]
            <= icrp_max_E
        )

        covered = site[in_range].copy()

        excluded_fluence = (
            site.loc[
                ~in_range,
                "Fluence (photons/m2/source)"]
                .sum()
        )

        if covered.empty:
            raise RuntimeError(
                f"No usable energy bins "
                f"for skeletal ID "
                f"{organ_id}."
            )

        # ====================================================
        # ENERGY AND FLUENCE
        # ====================================================

        E = covered["Energy Center (MeV)"].to_numpy(dtype=float)

        Phi = covered["Fluence (photons/m2/source)"].to_numpy(dtype=float)

        # ====================================================
        # ICRP RESPONSE FUNCTIONS
        # ====================================================

        R_AM = interpolate_response(
            E,
            response_df,
            "AM_Gy_m2"
        )

        R_TM50 = interpolate_response(
            E,
            response_df,
            "TM50_Gy_m2"
        )

        # ====================================================
        # DOSE
        # Phi is already bin-integrated.
        # Therefore:
        # D = sum(Phi * R)
        # ====================================================

        # RBM Dose
        if R_AM is None:
            marrow_dose = np.nan
        else:
            marrow_dose = np.sum(Phi * R_AM)

        # Endosteum Dose
        if R_TM50 is None:
            endosteum_dose = np.nan
        else:
            endosteum_dose = np.sum(Phi * R_TM50)

        # ====================================================
        # MASS
        # ====================================================

        if R_AM is None:
            marrow_mass_kg = np.nan
        else:
            marrow_mass_kg = (
                float(
                    skeletal_masses.loc[
                        organ_id,
                        marrow_column
                    ]
                )
                / 1000.0
            )

        if R_TM50 is None:
            endosteum_mass_kg = np.nan
        else:
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

            "Marrow dose (Gy/source)":
                marrow_dose,

            "Endosteum dose (Gy/source)":
                endosteum_dose,
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

    if rbm_results.empty:

        total_marrow_mass_kg = np.nan
        rbm_dose = np.nan

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


    # ========================================================
    # ENDOSTEUM MASS WEIGHTING
    # ========================================================

    endosteum_results = results[
        results["Endosteum dose (Gy/source)"].notna()
    ].copy()

    if endosteum_results.empty:

        total_endosteum_mass_kg = np.nan
        endosteum_dose = np.nan

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

    results.insert(0, "Phantom", config.PHANTOM_NAMES[phantom_code])

    results.insert(1, "Source Organ", source_organ)

    results.insert(2, "Source Type", source_type)

    results.insert(3, "Source Energy (MeV)",source_energy)

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
        "RBM SAF (kg^-1)"
    ] = (
        rbm_saf
    )

    results[
        "Endosteum Dose (Gy/source)"
    ] = (
        endosteum_dose
    )

    results[
        "Endosteum SAF (kg^-1)"
    ] = (
        endosteum_saf
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
                "Excluded fraction",
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

    return (
        results,
        rbm_dose,
        rbm_saf,
        endosteum_dose,
        endosteum_saf
    )


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

    response_functions = (load_icrp_response_functions())

    skeletal_masses = (load_skeletal_masses())

    # ========================================================
    # VALID ICRP RESPONSE-FUNCTION IDs
    # ========================================================

    rbm_ids = {
        organ_id
        for organ_id, response_df
        in response_functions.items()
        if response_df["AM_Gy_m2"].notna().any()
    }

    endosteum_ids = {
        organ_id
        for organ_id, response_df
        in response_functions.items()
        if response_df["TM50_Gy_m2"].notna().any()
    }

    # ========================================================
    # FIND FLUENCE FILES
    # ========================================================

    fluence_files = [
        f for f in find_fluence_files()
        if any(
            f"energy_{energy}_photon_fluence.csv" in f.name
            for energy in params["source_energies"]
        )
    ]

    if not fluence_files:
        raise FileNotFoundError(
            "No Geant4 photon-fluence CSV "
            "files were found in:\n"
            f"{config.GEANT4_GENERATED_INPUTS_DIR}"
        )

    print(
        f"\nFound {len(fluence_files)} "
        "photon-fluence file(s)."
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

        phantom_code = (match.group(1))

        source_organ = (match.group(2))

        geant4_source_type = (match.group(3).lower())

        source_energy = float(match.group(4))

        # ----------------------------------------------------
        # Convert Geant4 particle syntax
        #
        # gamma -> photon
        # e-    -> electron
        # ----------------------------------------------------

        source_type = (
            config.GEANT4_SOURCE_TYPE_MAP.get(
                geant4_source_type,
                geant4_source_type
            )
        )

        try:
            (
                result,
                rbm_dose,
                rbm_saf,
                endosteum_dose,
                endosteum_saf

            ) = calculate_from_fluence(

                fluence_file,

                response_functions,

                skeletal_masses,

                rbm_ids,

                endosteum_ids,

                source_energy,

                source_organ,

                source_type,

                phantom_code

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
    # OUTPUT FILE
    # ========================================================

    if (params["simulation_code"]== "GEANT4"):
        output_filename = ("j_geant4_rbm_endosteum_icrp116.csv")

    else:

        raise ValueError(
            "Unsupported simulation code: "
            f"{params['simulation_code']}"
        )

    output_file = (
        config.RESULTS_DIR
        / output_filename
    )

    combined_results.to_csv(
        output_file,
        index=False
    )

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

    print(
        f"Total result rows : "
        f"{len(combined_results)}"
    )

    print(
        f"Results saved to:\n"
        f"{output_file}"
    )

    return combined_results