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

Geant4 IDs are mapped to ICRP IDs as:

    1400 -> 14
    2500 -> 25
    2700 -> 27
    2900 -> 29
    4000 -> 40
    4200 -> 42
    4400 -> 44
    4600 -> 46
    4800 -> 48
    5000 -> 50
    5200 -> 52
    5400 -> 54
    5600 -> 56

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

from pathlib import Path
import re

import numpy as np
import pandas as pd


# ============================================================
# FILES
# ============================================================

FLUENCE_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/3_geant4/Internal/build/example_AM_2_photon_fluence.csv")

ICRP_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/ICRP116-Table-D-1-Updated.csv")

OUTPUT_FILE = Path(
    "RBM_ICRP116_results.csv"
)


# ============================================================
# GEANT4 ID -> ICRP ID
# ============================================================

SKELETAL_IDS = [
    1400,
    2500,
    2700,
    2900,
    4000,
    4200,
    4400,
    4600,
    4800,
    5000,
    5200,
    5400,
    5600,
]

# ============================================================
# ACTIVE MARROW MASSES
#
# These are in kg.
#
# IMPORTANT:
# Put the correct ICRP active-marrow masses here before
# using the final mass-weighted RBM SAF.
# ============================================================

AM_MASS_KG = {

    1400: .0269,
    2500: .0093,
    2700: .0889,
    2900: .0784,
    4000: .0094,
    4200: .2052,
    4400: .1888,
    4600: .0328,
    4800: .0456,
    5000: .1888,
    5200: .1439,
    5400: .1159,
    5600: .0363,

}


# ============================================================
# NORMALIZE COLUMN NAMES
# ============================================================

def normalize(text):

    return (
        str(text)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
        .replace("(", "")
        .replace(")", "")
    )


# ============================================================
# FIND COLUMN
# ============================================================

def find_column(
    dataframe,
    candidates
):

    normalized_columns = {
        normalize(column): column
        for column in dataframe.columns
    }

    for candidate in candidates:

        key = normalize(candidate)

        if key in normalized_columns:

            return normalized_columns[key]

    return None


# ============================================================
# LOAD GEANT4 FLUENCE
# ============================================================

print()
print("=" * 80)
print("LOADING GEANT4 PHOTON FLUENCE")
print("=" * 80)

fluence = pd.read_csv(
    FLUENCE_FILE
)


# ------------------------------------------------------------
# Find columns
# ------------------------------------------------------------

id_col = find_column(
    fluence,
    [
        "Spongiosa ID",
        "SpongiosaID"
    ]
)


energy_center_col = find_column(
    fluence,
    [
        "Energy Center (MeV)",
        "Energy Center"
    ]
)


fluence_col = find_column(
    fluence,
    [
        "Fluence (photons/m2/source)",
        "Fluence"
    ]
)


if id_col is None:

    raise ValueError(
        "Could not find the Spongiosa ID column."
    )


if energy_center_col is None:

    raise ValueError(
        "Could not find the Energy Center column."
    )


if fluence_col is None:

    raise ValueError(
        "Could not find the Fluence column."
    )


# ------------------------------------------------------------
# Convert data types
# ------------------------------------------------------------

fluence[id_col] = pd.to_numeric(
    fluence[id_col],
    errors="coerce"
)


fluence[energy_center_col] = pd.to_numeric(
    fluence[energy_center_col],
    errors="coerce"
)


fluence[fluence_col] = pd.to_numeric(
    fluence[fluence_col],
    errors="coerce"
).fillna(0.0)


fluence = fluence.dropna(
    subset=[
        id_col,
        energy_center_col
    ]
)


# Convert IDs to integer

fluence[id_col] = (
    fluence[id_col]
    .astype(int)
)


print(
    f"Loaded {len(fluence)} fluence rows."
)

# ============================================================
# LOAD ICRP 116 TABLE D.1
# ============================================================

print()
print("=" * 80)
print("LOADING ICRP 116 TABLE D.1")
print("=" * 80)

icrp_raw = pd.read_csv(
    ICRP_FILE,
    header=None
)


# ============================================================
# PARSE UPDATED ICRP CSV
#
# The CSV has this structure:
#
# Column 0                  Column 1                    Column 2
# Photon Energy (MeV)       AM                          TM50
#
# Row:
# [blank]                   Organ ID: 1400              Humeri...
#
# Then:
# Photon Energy (MeV)       AM                          TM50
# 0.01                      6.13E-16                    5.36E-16
# 0.015                     ...
#
# The Organ ID is in COLUMN 1.
# ============================================================

response_functions = {}


for row in range(
    icrp_raw.shape[0]
):

    # Organ ID is in column 1
    value = icrp_raw.iat[
        row,
        1
    ]


    if pd.isna(value):
        continue


    text = str(
        value
    ).strip()


    # --------------------------------------------------------
    # Look for:
    #
    # Organ ID: 1400
    # --------------------------------------------------------

    match = re.search(
        r"Organ\s*ID\s*:\s*(\d+)",
        text,
        re.IGNORECASE
    )


    if match is None:
        continue


    organ_id = int(
        match.group(1)
    )


    # --------------------------------------------------------
    # Organ name is in column 2
    # --------------------------------------------------------

    organ_name = str(
        icrp_raw.iat[
            row,
            2
        ]
    ).strip()


    # --------------------------------------------------------
    # Data structure:
    #
    # row     = Organ ID
    # row + 1 = Photon Energy / AM / TM50
    # row + 2 = first data point
    # --------------------------------------------------------

    data_start = row + 2


    energies = []
    am_values = []
    tm50_values = []


    current_row = data_start


    while (
        current_row
        < icrp_raw.shape[0]
    ):

        # ----------------------------------------------------
        # Photon energy is in column 0
        # ----------------------------------------------------

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


        # ----------------------------------------------------
        # AM is column 1
        # ----------------------------------------------------

        am_value = icrp_raw.iat[
            current_row,
            1
        ]


        try:

            am = float(
                str(am_value)
                .strip()
                .replace(
                    "E",
                    "e"
                )
            )

        except (
            ValueError,
            TypeError
        ):

            break


        # ----------------------------------------------------
        # TM50 is column 2
        # ----------------------------------------------------

        tm50_value = icrp_raw.iat[
            current_row,
            2
        ]


        try:

            tm50 = float(
                str(tm50_value)
                .strip()
                .replace(
                    "E",
                    "e"
                )
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


    # --------------------------------------------------------
    # Store response function
    # --------------------------------------------------------

    response_functions[
        organ_id
    ] = {

        "name":
            organ_name,

        "energy":
            np.array(
                energies,
                dtype=float
            ),

        "AM":
            np.array(
                am_values,
                dtype=float
            ),

        "TM50":
            np.array(
                tm50_values,
                dtype=float
            )

    }


# ============================================================
# PRINT WHAT WAS FOUND
# ============================================================

print()
print(
    "Spongiosa AM response functions found:"
)


for organ_id in sorted(
    response_functions
):

    data = response_functions[
        organ_id
    ]


    print(
        f"  Organ ID {organ_id:4d}: "
        f"{len(data['energy']):2d} energy points  "
        f"{data['name']}"
    )


# ============================================================
# CHECK REQUIRED ORGAN IDS
# ============================================================

missing_icrp_ids = [

    organ_id

    for organ_id
    in SKELETAL_IDS

    if organ_id
    not in response_functions

]


if missing_icrp_ids:

    raise RuntimeError(

        "The following required ICRP "
        "response functions were not found:\n"
        f"{missing_icrp_ids}"

    )


# ============================================================
# CHECK NUMBER OF ENERGY POINTS
# ============================================================

bad_energy_counts = {

    organ_id:

    len(
        response_functions[
            organ_id
        ]["energy"]
    )

    for organ_id
    in SKELETAL_IDS

    if len(
        response_functions[
            organ_id
        ]["energy"]
    ) != 25

}


if bad_energy_counts:

    raise RuntimeError(

        "Unexpected number of energy points "
        "for these organ IDs:\n"
        f"{bad_energy_counts}"

    )


print()
print(
    "All 13 required ICRP spongiosa "
    "AM response functions were found."
)

# ============================================================
# LOG-LOG INTERPOLATION
# ============================================================

def make_response_interpolator(
    energy,
    response
):

    energy = np.asarray(
        energy,
        dtype=float
    )

    response = np.asarray(
        response,
        dtype=float
    )


    # --------------------------------------------------------
    # Keep only valid positive values
    # --------------------------------------------------------

    mask = (
        np.isfinite(energy)
        &
        np.isfinite(response)
        &
        (energy > 0)
        &
        (response > 0)
    )


    energy = energy[mask]

    response = response[mask]


    # --------------------------------------------------------
    # Sort by energy
    # --------------------------------------------------------

    order = np.argsort(
        energy
    )

    energy = energy[order]

    response = response[order]


    if len(energy) < 2:

        raise ValueError(
            "Not enough valid ICRP response points."
        )


    # --------------------------------------------------------
    # Log-log interpolation
    # --------------------------------------------------------

    log_energy = np.log(
        energy
    )

    log_response = np.log(
        response
    )


    def interpolate(E):

        E = np.asarray(
            E,
            dtype=float
        )


        if (
            np.any(
                E < energy.min()
            )
            or
            np.any(
                E > energy.max()
            )
        ):

            raise ValueError(

                "Energy outside ICRP "
                "Table D.1 range: "

                f"{energy.min()} - "
                f"{energy.max()} MeV."

            )


        return np.exp(

            np.interp(
                np.log(E),
                log_energy,
                log_response
            )

        )


    return (
        interpolate,
        energy.min(),
        energy.max()
    )

# ============================================================
# CALCULATE SITE DOSES
# ============================================================

print()
print("=" * 80)
print("CALCULATING SKELETAL SITE DOSES")
print("=" * 80)


site_results = []


for organ_id in SKELETAL_IDS:

    print(f"Processing Organ ID {organ_id}...")

    # --------------------------------------------------------
    # Select this skeletal region
    # --------------------------------------------------------

    site_data = fluence[fluence[id_col] == organ_id].copy()

    if site_data.empty:

        print(
            "  WARNING: no fluence data."
        )

        continue


    # --------------------------------------------------------
    # Get ICRP response function for this skeletal site
    # --------------------------------------------------------

    if organ_id not in response_functions:

        print(
            f"WARNING: no ICRP AM response "
            f"found for ID {organ_id}."
        )

        continue


    icrp_energy = (
        response_functions[
            organ_id
        ]["energy"]
    )

    response = (
        response_functions[
            organ_id
        ]["AM"]
    )

    # --------------------------------------------------------
    # Build interpolation function
    # --------------------------------------------------------

    (
        interpolate_R,
        icrp_Emin,
        icrp_Emax
    ) = make_response_interpolator(
        icrp_energy,
        response
    )


    print(
        f"  ICRP energy range: "
        f"{icrp_Emin:.3f} - "
        f"{icrp_Emax:.3f} MeV"
    )
    # --------------------------------------------------------
    # Keep only Geant4 bins within the ICRP range
    # --------------------------------------------------------

    in_range = (

        (
            site_data[
                energy_center_col
            ]
            >= icrp_Emin
        )

        &

        (
            site_data[
                energy_center_col
            ]
            <= icrp_Emax
        )

    )


    # Fluence outside ICRP range

    excluded_fluence = (
        site_data.loc[
            ~in_range,
            fluence_col
        ]
        .sum()
    )


    # Usable bins

    used = site_data.loc[
        in_range
    ].copy()


    if used.empty:

        print(
            "  WARNING: no usable "
            "energy bins."
        )

        continue


    # --------------------------------------------------------
    # Energy and fluence arrays
    # --------------------------------------------------------

    E = (
        used[
            energy_center_col
        ]
        .to_numpy(
            dtype=float
        )
    )


    Phi = (
        used[
            fluence_col
        ]
        .to_numpy(
            dtype=float
        )
    )


    # --------------------------------------------------------
    # Interpolate ICRP response
    # --------------------------------------------------------

    R = interpolate_R(
        E
    )


    # --------------------------------------------------------
    # Dose contribution
    #
    # Phi:
    #
    # photons / m2 / source
    #
    # R:
    #
    # Gy m2 / photon
    #
    # Therefore:
    #
    # Phi * R:
    #
    # Gy / source
    #
    # --------------------------------------------------------

    dose_contribution = (
        Phi * R
    )


    dose_Gy_per_source = (
        np.sum(
            dose_contribution
        )
    )


    total_fluence = (
        np.sum(
            Phi
        )
    )


    # --------------------------------------------------------
    # Save site result
    # --------------------------------------------------------

    site_results.append(

        {
            "Organ ID":
                organ_id,

            "Skeletal site":
                response_functions[organ_id]["name"],

            "Total fluence (photons/m2/source)":
                total_fluence,

            "Excluded fluence outside Table D.1 range":
                excluded_fluence,

            "AM dose (Gy/source)":
                dose_Gy_per_source,

            "AM mass (kg)":
                AM_MASS_KG.get(organ_id,np.nan),
        }

    )


    print(
        f"  Total fluence: "
        f"{total_fluence:.6e}"
    )


    print(
        f"  AM dose: "
        f"{dose_Gy_per_source:.6e} Gy/source"
    )


# ============================================================
# CREATE RESULTS DATAFRAME
# ============================================================

results = pd.DataFrame(
    site_results
)


if results.empty:

    raise RuntimeError(
        "No skeletal-site results "
        "were calculated."
    )


# ============================================================
# CHECK ACTIVE MARROW MASSES
# ============================================================

missing_masses = (
    results.loc[
        results["AM mass (kg)"].isna(),
        "Organ ID"
    ]
    .tolist()
)

if missing_masses:

    print()
    print("=" * 80)
    print("ACTIVE MARROW MASSES ARE MISSING")
    print("=" * 80)

    print()

    print(
        "The site-specific AM doses "
        "were successfully calculated."
    )

    print()

    print(
        "However, the final mass-weighted "
        "RBM dose cannot yet be calculated."
    )

    print()

    print(
        "Add the active-marrow masses "
        "to AM_MASS_KG near the top "
        "of this script."
    )

    print()

    print(
        "Missing ICRP IDs:"
    )

    print(
        missing_masses
    )

    print()

    print(
        "Site results:"
    )

    print()

    print(

        results[
            [
                "Organ ID",
                "Skeletal site",
                "Total fluence (photons/m2/source)",
                "AM dose (Gy/source)",
                "AM mass (kg)"
            ]
        ]
        .to_string(
            index=False,
            float_format=lambda x:
                f"{x:.6e}"
        )

    )

    # Save site-level results anyway.

    results.to_csv(
        OUTPUT_FILE,
        index=False
    )

    print()

    print(
        f"Site results saved to: "
        f"{OUTPUT_FILE}"
    )

    raise SystemExit


# ============================================================
# MASS-WEIGHTED RBM DOSE
# ============================================================

total_AM_mass = (
    results[
        "AM mass (kg)"
    ]
    .sum()
)


results[
    "Mass fraction"
] = (

    results[
        "AM mass (kg)"
    ]

    /

    total_AM_mass

)


results[
    "Mass-weighted dose contribution (Gy/source)"
] = (

    results[
        "AM dose (Gy/source)"
    ]

    *

    results[
        "Mass fraction"
    ]

)


# ------------------------------------------------------------
# Total RBM dose
# ------------------------------------------------------------

RBM_dose = (
    results[
        "Mass-weighted dose contribution (Gy/source)"
    ]
    .sum()
)


# ============================================================
# RBM SAF
# ============================================================

# 1 MeV photon source

MeV_to_J = (
    1.602176634e-13
)

RBM_SAF = (
    RBM_dose
    /
    MeV_to_J
)


# ============================================================
# SAVE RESULTS
# ============================================================

results.to_csv(
    OUTPUT_FILE,
    index=False
)

# ============================================================
# PRINT FINAL RESULTS
# ============================================================

print()
print()
print("=" * 100)
print("ICRP 116 TABLE D.1 RBM CALCULATION")
print("=" * 100)

print()

print(
    results[
        [
            "Organ ID",
            "Skeletal site",
            "Total fluence (photons/m2/source)",
            "AM dose (Gy/source)",
            "AM mass (kg)",
            "Mass fraction"
        ]
    ]
    .to_string(
        index=False,
        float_format=lambda x:
            f"{x:.6e}"
    )
)

print()
print("-" * 100)

print(
    f"Total active marrow mass : "
    f"{total_AM_mass:.6e} kg"
)

print(
    f"RBM dose                 : "
    f"{RBM_dose:.6e} Gy/source"
)

print(
    f"RBM SAF                  : "
    f"{RBM_SAF:.6e} kg^-1"
)

print(
    "-" * 100
)

print()

print(
    f"Results saved to: "
    f"{OUTPUT_FILE}"
)

print()