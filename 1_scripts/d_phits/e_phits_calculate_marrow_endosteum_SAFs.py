"""
Calculate active red bone marrow (RBM) and endosteum
dose and SAF from PHITS T-Track photon fluence.

1. PHITS energy-binned photon fluence in skeletal tissues.
2. ICRP Publication 116, Table D.1 photon response functions.

PHITS T-Track:
    unit = 1
    Flux = 1/cm^2/source

Therefore:

    Phi [1/m^2/source]
        =
    Phi [1/cm^2/source] * 1e4

The PHITS fluence values are histogram/bin fluences,
so they are NOT multiplied by Delta E again.

ICRP response functions are interpolated in
log(E)-log(R) space.

RBM and endosteum are mass-weighted independently
over the skeletal regions for which the corresponding
ICRP response function exists.

For a photon source:

    SAF = dose / emitted energy

where:

    emitted energy =
        source energy (MeV)
        * 1.6021766339999e-13 J/MeV
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

# PHITS T-Track unit=1 gives 1/cm^2/source.
# ICRP response functions use m^2.
CM2_TO_M2 = 1.0e4

# ============================================================
# Configs
# ============================================================

mass_file = config.SKELETAL_MASSES_CSV
phits_generated_inputs_dir = config.GENERATED_INPUTS_DIR
fluence_to_dose_response_functions = config.SKELETAL_RESPONSE_FUNCTIONS_CSV
phantom_names = config.PHANTOM_NAMES
phits_results_dir = config.RESULTS_PHITS_DIR
phits_output_fluence = "e_phits_rbm_endosteum_icrp116.csv"

# ============================================================
# PHITS FLUENCE FILENAME PATTERN
# ============================================================

fluence_filename_pattern = re.compile(
    r"phits_fluence_MRCP_"
    r"(AM|AF)_source_"
    r"(.+?)_"
    r"(photon|electron)_energy_"
    r"([0-9Ee.+-]+)"
    r"\.out$",
    re.IGNORECASE
)


# ============================================================
# LOAD SKELETAL TISSUE MASSES
# ============================================================

def load_skeletal_masses():

    if not mass_file.exists():

        raise FileNotFoundError(
            f"Could not find skeletal tissue mass file:\n"
            f"{mass_file}"
        )

    masses = pd.read_csv(
        mass_file
    )

    required_columns = [
        "Organ ID",

        "Ref_AM_Marrow_Mass(g)",
        "Ref_AM_Endosteum_Mass(g)",

        "Ref_AF_Marrow_Mass(g)",
        "Ref_AF_Endosteum_Mass(g)",
    ]

    missing = [
        column
        for column in required_columns
        if column not in masses.columns
    ]

    if missing:

        raise ValueError(
            "Missing columns in skeletal tissue mass CSV:\n"
            + "\n".join(
                f"  {column}"
                for column in missing
            )
        )

    masses["Organ ID"] = pd.to_numeric(
        masses["Organ ID"],
        errors="coerce"
    )

    masses = masses.dropna(
        subset=["Organ ID"]
    )

    masses["Organ ID"] = (
        masses["Organ ID"].astype(int)
    )

    masses = masses.set_index(
        "Organ ID"
    )

    return masses


# ============================================================
# FIND PHITS FLUENCE FILES
# ============================================================

def find_fluence_files():

    root = phits_generated_inputs_dir

    return sorted(
        root.rglob(
            "phits_fluence_MRCP_*_photon_energy_*.out"
        )
    )


# ============================================================
# LOAD ICRP 116 RESPONSE FUNCTIONS
# ============================================================

def load_icrp_response_functions():

    icrp_file = fluence_to_dose_response_functions

    if not icrp_file.exists():

        raise FileNotFoundError(
            f"Could not find ICRP 116 Table D.1:\n"
            f"{icrp_file}"
        )

    icrp_raw = pd.read_csv(
        icrp_file,
        header=None
    )

    response_functions = {}

    # --------------------------------------------------------
    # The ICRP CSV structure is:
    #
    # Column 0 = Photon energy
    # Column 1 = AM response
    # Column 2 = TM50 response
    #
    # Organ ID is identified in column 1.
    # --------------------------------------------------------

    for row in range(
        icrp_raw.shape[0]
    ):

        value = icrp_raw.iat[
            row,
            1
        ]

        if pd.isna(value):
            continue

        match = re.search(
            r"Organ\s*ID\s*:\s*(\d+)",
            str(value),
            re.IGNORECASE
        )

        if match is None:
            continue

        organ_id = int(
            match.group(1)
        )

        # ----------------------------------------------------
        # First data row
        #
        # row     = Organ ID
        # row + 1 = column headers
        # row + 2 = first data point
        # ----------------------------------------------------

        data_start = row + 2

        energies = []
        am_values = []
        tm50_values = []

        current_row = data_start

        while current_row < icrp_raw.shape[0]:

            energy_value = icrp_raw.iat[
                current_row,
                0
            ]

            try:

                energy = float(
                    energy_value
                )

            except (
                ValueError,
                TypeError
            ):

                break

            # ------------------------------------------------
            # AM
            # ------------------------------------------------

            am_value = icrp_raw.iat[
                current_row,
                1
            ]

            try:

                am = float(
                    str(am_value)
                    .strip()
                    .replace("E", "e")
                )

            except (
                ValueError,
                TypeError
            ):

                am = np.nan

            # ------------------------------------------------
            # TM50
            # ------------------------------------------------

            tm50_value = icrp_raw.iat[
                current_row,
                2
            ]

            try:

                tm50 = float(
                    str(tm50_value)
                    .strip()
                    .replace("E", "e")
                )

            except (
                ValueError,
                TypeError
            ):

                tm50 = np.nan

            energies.append(
                energy
            )

            am_values.append(
                am
            )

            tm50_values.append(
                tm50
            )

            current_row += 1

        # ----------------------------------------------------
        # Store response function
        # ----------------------------------------------------

        response_df = pd.DataFrame({

            "Energy_MeV":
                energies,

            "AM_Gy_m2":
                am_values,

            "TM50_Gy_m2":
                tm50_values,

        })

        # Keep the region if at least one response function
        # exists.
        response_df = response_df[
            (
                response_df["AM_Gy_m2"].notna()
            )
            |
            (
                response_df["TM50_Gy_m2"].notna()
            )
        ].reset_index(
            drop=True
        )

        response_functions[
            organ_id
        ] = response_df

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
        [
            "Energy_MeV",
            response_column
        ]
    ].copy()

    # Remove NA values for THIS response function only.
    response_data = response_data.dropna()

    response_data = response_data[
        (
            response_data["Energy_MeV"] > 0
        )
        &
        (
            response_data[response_column] > 0
        )
    ]

    response_data = response_data.sort_values(
        "Energy_MeV"
    )

    # No response function for this target.
    if response_data.empty:

        return None

    table_E = response_data[
        "Energy_MeV"
    ].to_numpy(
        dtype=float
    )

    table_R = response_data[
        response_column
    ].to_numpy(
        dtype=float
    )

    if len(table_E) < 2:

        raise RuntimeError(
            f"Fewer than two valid ICRP "
            f"response points for "
            f"{response_column}."
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
# PARSE ONE PHITS T-TRACK REGION
# ============================================================

def parse_phits_region(
    text,
    region_id
):

    # Find the region section.
    region_pattern = re.compile(
        rf"newpage:\s*"
        rf"\n#\s*no\.\s*=\s*\d+"
        rf"\s+reg\s*=\s*{region_id}\b"
        rf"(.*?)(?=\n#\s*-{{10,}}|\Z)",
        re.IGNORECASE |
        re.DOTALL
    )

    match = region_pattern.search(
        text
    )

    if match is None:

        raise RuntimeError(
            f"Could not find PHITS T-Track "
            f"section for region {region_id}."
        )

    section = match.group(1)

    # --------------------------------------------------------
    # Find the actual energy table.
    #
    # Data begin after:
    #
    # # e-lower      e-upper      photon    r.err
    #
    # and end at:
    #
    # # sum over
    # --------------------------------------------------------

    table_match = re.search(
        r"#\s*e-lower\s+"
        r"e-upper\s+"
        r"photon\s+"
        r"r\.err"
        r"(.*?)(?=\n\s*#\s*sum\s+over)",
        section,
        re.IGNORECASE |
        re.DOTALL
    )

    if table_match is None:

        raise RuntimeError(
            f"Could not find T-Track energy "
            f"table for region {region_id}."
        )

    table_text = table_match.group(1)

    rows = []

    # --------------------------------------------------------
    # Each data row has:
    #
    # E_low E_high fluence relative_error
    # --------------------------------------------------------

    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"

    row_pattern = re.compile(
        rf"^\s*"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s*$",
        re.MULTILINE
    )

    for match in row_pattern.finditer(
        table_text
    ):

        energy_low = float(
            match.group(1)
        )

        energy_high = float(
            match.group(2)
        )

        fluence_cm2 = float(
            match.group(3)
        )

        relative_error = float(
            match.group(4)
        )

        rows.append({

            "Skeletal ID":
                region_id,

            "Energy Low (MeV)":
                energy_low,

            "Energy High (MeV)":
                energy_high,

            "Energy Center (MeV)":
                np.sqrt(
                    energy_low
                    * energy_high
                ),

            "Fluence (1/cm2/source)":
                fluence_cm2,

            "Relative Error":
                relative_error,

        })

    if not rows:

        raise RuntimeError(
            f"No T-Track energy bins found "
            f"for region {region_id}."
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# LOAD PHITS FLUENCE FILE
# ============================================================

def load_phits_fluence(
    fluence_file
):

    text = fluence_file.read_text(
        encoding="utf-8",
        errors="ignore"
    )

    # Find regions from the T-Track file itself.
    region_matches = re.findall(
        r"newpage:\s*"
        r"\n#\s*no\.\s*=\s*\d+"
        r"\s+reg\s*=\s*(\d+)",
        text,
        re.IGNORECASE
    )

    if not region_matches:

        raise RuntimeError(
            f"No PHITS T-Track regions found in:\n"
            f"{fluence_file}"
        )

    region_ids = [
        int(region_id)
        for region_id in region_matches
    ]

    all_regions = []

    for region_id in region_ids:

        region_data = parse_phits_region(
            text,
            region_id
        )

        all_regions.append(
            region_data
        )

    return pd.concat(
        all_regions,
        ignore_index=True
    )


# ============================================================
# CALCULATE ONE PHITS FLUENCE FILE
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
        f"Processing PHITS fluence file:\n"
        f"{fluence_file.name}"
    )
    print("=" * 90)

    # ========================================================
    # LOAD PHITS T-TRACK DATA
    # ========================================================

    fluence = load_phits_fluence(
        fluence_file
    )

    # ========================================================
    # SELECT MASS COLUMNS
    # ========================================================

    if phantom_code == "AM":

        marrow_column = (
            "Ref_AM_Marrow_Mass(g)"
        )

        endosteum_column = (
            "Ref_AM_Endosteum_Mass(g)"
        )

    elif phantom_code == "AF":

        marrow_column = (
            "Ref_AF_Marrow_Mass(g)"
        )

        endosteum_column = (
            "Ref_AF_Endosteum_Mass(g)"
        )

    else:

        raise ValueError(
            f"Unsupported phantom code: "
            f"{phantom_code}"
        )

    # ========================================================
    # SKELETAL IDS
    # ========================================================

    skeletal_ids = (
        rbm_ids
        |
        endosteum_ids
    )

    # ========================================================
    # CHECK MASSES
    # ========================================================

    missing_mass_ids = [

        organ_id

        for organ_id
        in skeletal_ids

        if organ_id
        not in skeletal_masses.index

    ]

    if missing_mass_ids:

        raise RuntimeError(
            "Missing skeletal tissue masses "
            f"for IDs: {missing_mass_ids}"
        )

    # ========================================================
    # CALCULATE EACH SKELETAL SITE
    # ========================================================

    site_results = []

    for organ_id in sorted(
        skeletal_ids
    ):

        # ----------------------------------------------------
        # Fluence for this skeletal site
        # ----------------------------------------------------

        site = fluence[
            fluence["Skeletal ID"]
            == organ_id
        ].copy()

        if site.empty:

            raise RuntimeError(
                f"No PHITS T-Track fluence "
                f"data for skeletal ID "
                f"{organ_id}."
            )

        # ----------------------------------------------------
        # ICRP response
        # ----------------------------------------------------

        if organ_id not in response_functions:

            raise RuntimeError(
                f"No ICRP response function "
                f"for skeletal ID "
                f"{organ_id}."
            )

        response_df = (
            response_functions[
                organ_id
            ]
        )

        # ====================================================
        # TOTAL FLUENCE
        # ====================================================

        total_fluence_cm2 = (site["Fluence (1/cm2/source)"].sum())

        total_fluence_m2 = (total_fluence_cm2 * CM2_TO_M2)
        
        # ====================================================
        # CONVERT cm^-2 -> m^-2
        # ====================================================

        site["Fluence (1/m2/source)"] = (
            site["Fluence (1/cm2/source)"] * CM2_TO_M2)

        # ====================================================
        # RBM
        # ====================================================

        rbm_site = site.copy()

        table_E = response_df[
            response_df["AM_Gy_m2"].notna()
        ]["Energy_MeV"]

        if table_E.empty:

            marrow_dose = np.nan
            marrow_mass_kg = np.nan
            excluded_marrow_fluence_m2 = np.nan
            marrow_relative_error = np.nan
            marrow_statistical_uncertainty = np.nan

        else:

            rbm_in_range = (
                (rbm_site["Energy Center (MeV)"] >= table_E.min())
                &
                (rbm_site["Energy Center (MeV)"] <= table_E.max())
            )

            rbm_covered = rbm_site[
                rbm_in_range
            ]

            if rbm_covered.empty:

                raise RuntimeError(
                    f"No PHITS energy bins inside ICRP AM range "
                    f"for region {organ_id}."
                )

            R_AM = interpolate_response(
                rbm_covered[
                    "Energy Center (MeV)"
                ].to_numpy(dtype=float),
                response_df,
                "AM_Gy_m2"
            )

            Phi = rbm_covered[
                "Fluence (1/m2/source)"
            ].to_numpy(dtype=float)

            relative_error = rbm_covered[
                "Relative Error"
            ].to_numpy(dtype=float)

            dose_contribution = (
                Phi * R_AM
            )

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

            marrow_statistical_uncertainty = (
                marrow_relative_error * 100
                if np.isfinite(marrow_relative_error)
                else np.nan
            )

            excluded_marrow_fluence_m2 = (
                rbm_site.loc[
                    ~rbm_in_range,
                    "Fluence (1/cm2/source)"
                ].sum()
            ) * CM2_TO_M2

            marrow_mass_kg = (
                float(
                    skeletal_masses.loc[
                        organ_id,
                        marrow_column
                    ]
                )
                / 1000.0
            )

        # ====================================================
        # ENDOSTEUM
        # ====================================================

        endo_site = site.copy()

        table_E = response_df[
            response_df["TM50_Gy_m2"].notna()
        ]["Energy_MeV"]

        if table_E.empty:

            endosteum_dose = np.nan
            endosteum_mass_kg = np.nan
            excluded_endosteum_fluence_m2 = np.nan
            endosteum_relative_error = np.nan
            endosteum_statistical_uncertainty = np.nan

        else:

            endo_in_range = (
                (endo_site["Energy Center (MeV)"] >= table_E.min())
                &
                (endo_site["Energy Center (MeV)"] <= table_E.max())
            )

            endo_covered = endo_site[
                endo_in_range
            ]

            if endo_covered.empty:

                raise RuntimeError(
                    f"No PHITS energy bins inside ICRP TM50 range "
                    f"for region {organ_id}."
                )

            R_TM50 = interpolate_response(
                endo_covered[
                    "Energy Center (MeV)"
                ].to_numpy(dtype=float),
                response_df,
                "TM50_Gy_m2"
            )

            Phi = endo_covered[
                "Fluence (1/m2/source)"
            ].to_numpy(dtype=float)

            relative_error = endo_covered[
                "Relative Error"
            ].to_numpy(dtype=float)

            dose_contribution = (
                Phi * R_TM50
            )

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

            endosteum_statistical_uncertainty = (
                endosteum_relative_error * 100
                if np.isfinite(endosteum_relative_error)
                else np.nan
            )

            excluded_endosteum_fluence_m2 = (
                endo_site.loc[
                    ~endo_in_range,
                    "Fluence (1/cm2/source)"
                ].sum()
            ) * CM2_TO_M2

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
                total_fluence_m2,

            "Marrow excluded fluence (photons/m2/source)":
                excluded_marrow_fluence_m2,

            "Endosteum excluded fluence (photons/m2/source)":
                excluded_endosteum_fluence_m2,

            "Marrow dose (Gy/source)":
                marrow_dose,

            "Marrow Relative Error":
                marrow_relative_error,

            "Marrow Statistical Uncertainty (%)":
                marrow_statistical_uncertainty,

            "Endosteum dose (Gy/source)":
                endosteum_dose,

            "Endosteum Relative Error":
                endosteum_relative_error,

            "Endosteum Statistical Uncertainty (%)":
                endosteum_statistical_uncertainty,

        })

    # ========================================================
    # RESULTS
    # ========================================================

    results = pd.DataFrame(
        site_results
    )

    # ========================================================
    # RBM MASS WEIGHTING
    # ========================================================

    rbm_results = results[
        results[
            "Marrow dose (Gy/source)"
        ].notna()
    ].copy()

    if rbm_results.empty:

        total_marrow_mass_kg = np.nan
        rbm_dose = np.nan
        rbm_sigma = np.nan
        rbm_relative_error = np.nan
        rbm_statistical_uncertainty = np.nan

    else:

        rbm_results = rbm_results[rbm_results["Marrow mass (kg)"].notna()].copy()
        total_marrow_mass_kg = (
            rbm_results[
                "Marrow mass (kg)"
            ].sum()
        )

        rbm_results[
            "Marrow mass fraction"
        ] = (
            rbm_results[
                "Marrow mass (kg)"
            ]
            / total_marrow_mass_kg
        )

        rbm_results[
            "Mass-weighted marrow dose contribution (Gy/source)"
        ] = (
            rbm_results[
                "Marrow dose (Gy/source)"
            ]
            *
            rbm_results[
                "Marrow mass fraction"
            ]
        )

        # Total RBM dose FIRST
        rbm_dose = (
            rbm_results[
                "Mass-weighted marrow dose contribution (Gy/source)"
            ].sum()
        )

        # Then calculate uncertainty
        rbm_sigma = np.sqrt(
            np.sum(
                (
                    rbm_results[
                        "Marrow mass fraction"
                    ]
                    *
                    rbm_results[
                        "Marrow dose (Gy/source)"
                    ]
                    *
                    rbm_results[
                        "Marrow Relative Error"
                    ]
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
        results[
            "Endosteum dose (Gy/source)"
        ].notna()
    ].copy()

    if endosteum_results.empty:

        total_endosteum_mass_kg = np.nan
        endosteum_dose = np.nan
        endosteum_sigma = np.nan
        endosteum_relative_error = np.nan
        endosteum_statistical_uncertainty = np.nan

    else:

        endosteum_results = endosteum_results[endosteum_results["Endosteum mass (kg)"].notna()].copy()
        total_endosteum_mass_kg = (
            endosteum_results[
                "Endosteum mass (kg)"
            ].sum()
        )

        endosteum_results[
            "Endosteum mass fraction"
        ] = (
            endosteum_results[
                "Endosteum mass (kg)"
            ]
            /
            total_endosteum_mass_kg
        )

        endosteum_results[
            "Mass-weighted endosteum dose contribution (Gy/source)"
        ] = (
            endosteum_results[
                "Endosteum dose (Gy/source)"
            ]
            *
            endosteum_results[
                "Endosteum mass fraction"
            ]
        )

        # Total endosteum dose FIRST
        endosteum_dose = (
            endosteum_results[
                "Mass-weighted endosteum dose contribution (Gy/source)"
            ].sum()
        )

        # Then calculate uncertainty
        endosteum_sigma = np.sqrt(
            np.sum(
                (
                    endosteum_results[
                        "Endosteum mass fraction"
                    ]
                    *
                    endosteum_results[
                        "Endosteum dose (Gy/source)"
                    ]
                    *
                    endosteum_results[
                        "Endosteum Relative Error"
                    ]
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
    # MASS FRACTIONS BACK INTO RESULTS
    # ========================================================

    results[
        "Marrow mass fraction"
    ] = np.nan

    results[
        "Endosteum mass fraction"
    ] = np.nan

    results.loc[
        rbm_results.index,
        "Marrow mass fraction"
    ] = rbm_results[
        "Marrow mass fraction"
    ]

    results.loc[
        endosteum_results.index,
        "Endosteum mass fraction"
    ] = endosteum_results[
        "Endosteum mass fraction"
    ]

    # ========================================================
    # SAF
    # ========================================================

    emitted_energy_J = (
        source_energy
        *
        MeV_to_J
    )

    rbm_saf = (
        rbm_dose
        /
        emitted_energy_J
    )

    endosteum_saf = (
        endosteum_dose
        /
        emitted_energy_J
    )

    rbm_saf_sigma = (
        rbm_sigma
        /
        emitted_energy_J
    )

    endosteum_saf_sigma = (
        endosteum_sigma
        /
        emitted_energy_J
    )

    # ========================================================
    # SOURCE INFORMATION
    # ========================================================

    results.insert(0, "Phantom", phantom_names[phantom_code])
    results.insert(1, "Source Organ", source_organ)
    results.insert(2, "Source Type", source_type)
    results.insert(3, "Source Energy (MeV)", source_energy)

    # ========================================================
    # TOTALS
    # ========================================================

    results[
        "Total Active Marrow Mass (kg)"
    ] = total_marrow_mass_kg

    results[
        "Total Endosteum Mass (kg)"
    ] = total_endosteum_mass_kg

    results[
        "RBM Dose (Gy/source)"
    ] = rbm_dose

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
        "RBM SAF Statistical Uncertainty (kg^-1)"
    ] = rbm_saf_sigma

    results[
        "Endosteum Dose (Gy/source)"
    ] = endosteum_dose

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
        "Endosteum SAF Statistical Uncertainty (kg^-1)"
    ] = endosteum_saf_sigma

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

def phits_calculate_marrow_endosteum_SAFs(
    params
):

    # --------------------------------------------------------
    # Only calculate for photon sources
    # --------------------------------------------------------

    if (
        params["source_type"].lower()
        not in {"photon", "gamma"}
    ):

        print(
            f"\nSkipping ICRP 116 skeletal "
            f"calculation for source type: "
            f"{params['source_type']}"
        )

        print(
            "ICRP 116 fluence-to-dose "
            "response functions implemented "
            "here are for photon sources."
        )

        return None

    print()
    print("=" * 90)

    print(
        "ICRP 116 RBM AND ENDOSTEUM "
        "CALCULATION FROM PHITS T-TRACK"
    )

    print("=" * 90)

    # ========================================================
    # LOAD DATA
    # ========================================================

    response_functions = (
        load_icrp_response_functions()
    )

    skeletal_masses = (
        load_skeletal_masses()
    )

    # ========================================================
    # VALID RESPONSE-FUNCTION IDs
    # ========================================================

    rbm_ids = {

        organ_id

        for organ_id, response_df
        in response_functions.items()

        if response_df[
            "AM_Gy_m2"
        ].notna().any()

    }

    endosteum_ids = {

        organ_id

        for organ_id, response_df
        in response_functions.items()

        if response_df[
            "TM50_Gy_m2"
        ].notna().any()

    }

    print(
        f"\nICRP RBM regions        : "
        f"{len(rbm_ids)}"
    )

    print(
        f"ICRP endosteum regions : "
        f"{len(endosteum_ids)}"
    )

    # ========================================================
    # FIND PHITS FLUENCE FILES
    # ========================================================

    fluence_files = [
        f
        for f in find_fluence_files()

        if any(
            f"energy_{energy}.out"
            in f.name

            for energy
            in params["source_energies"]
        )
    ]

    if not fluence_files:

        raise FileNotFoundError(
            "No PHITS photon-fluence "
            "output files were found in:\n"
            f"{phits_generated_inputs_dir}"
        )

    print(
        f"\nFound {len(fluence_files)} "
        "PHITS photon-fluence file(s)."
    )

    # ========================================================
    # PROCESS EVERY PHITS FLUENCE FILE
    # ========================================================

    all_results = []

    for fluence_file in fluence_files:

        match = (
            fluence_filename_pattern.match(
                fluence_file.name
            )
        )

        if match is None:

            print(
                f"\n[WARNING] Could not parse "
                f"fluence filename:\n"
                f"  {fluence_file.name}"
            )

            continue

        phantom_code = (
            match.group(1)
        )

        source_organ = (
            match.group(2)
        )

        source_type = (
            match.group(3).lower()
        )

        source_energy = float(
            match.group(4)
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
    # COMBINE
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

    # ============================================================
    # Keep overall simulation results only on the first row
    # ============================================================

    summary_columns = [
        "Total Active Marrow Mass (kg)",
        "Total Endosteum Mass (kg)",

        "RBM Dose (Gy/source)",
        "RBM Relative Error",
        "RBM Statistical Uncertainty (%)",
        "RBM SAF (kg^-1)",
        "RBM SAF Statistical Uncertainty (kg^-1)",

        "Endosteum Dose (Gy/source)",
        "Endosteum Relative Error",
        "Endosteum Statistical Uncertainty (%)",
        "Endosteum SAF (kg^-1)",
        "Endosteum SAF Statistical Uncertainty (kg^-1)",
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

        # Clear repeated values from remaining rows
        remaining_indices = group.index[1:]

        combined_results.loc[
            remaining_indices,
            summary_columns
        ] = np.nan

    # ========================================================
    # OUTPUT FILE
    # ========================================================

    if (params["simulation_code"] == "PHITS"):
        output_filename = phits_output_fluence
        

    else:

        raise ValueError(
            "Unsupported simulation code: "
            f"{params['simulation_code']}"
        )

    output_file = (phits_results_dir / output_filename)

    combined_results.to_csv(
        output_file,
        index=False
    )

    # ========================================================
    # FINAL MESSAGE
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