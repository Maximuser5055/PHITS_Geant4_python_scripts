# This script compares the results of this simulation against the literature

# Import sub-Python files
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import b_config.a_config as config

# ============================================================
# FILE LOCATIONS
# ============================================================

LITERATURE_FILE = config.OTHER_INPUT_FILES_DIR / "literature_1_comparison_source_liver.csv"
AM_FILE = config.RESULTS_DIR / "k_geant4_target_regions_dose_SAFs_AM.csv"
AF_FILE = config.RESULTS_DIR / "l_geant4_target_regions_dose_SAFs_AF.csv"
OUTPUT_COMPARISON_FILE = config.RESULTS_DIR / "q_liver_SAF_1_MeV_comparison_literature_1.csv"
OUTPUT_PNG = config.RESULTS_DIR/ "r_liver_SAF_1_MeV_comparison.png"

ENERGY_MEV = 1.0
SOURCE_ORGAN = "Liver"

# ============================================================
# TARGET ORGAN MAPPING
# ============================================================

# Left side:
#   Target region name used in the literature
# Right side:
#   Target region name used in the results

TARGETS = {
    "Red bone marrow": "Red (active) marrow",
    "Colon": "Colon wall",
    "Lungs": "RLung + LLung",
    "Stomach": "Stomach wall",
    "Breast": "Breast-a + Breast-g",
    "Esophagus": "Oesophagus wall",
    "Liver": "Liver",
    "Thyroid": "Thyroid",
    "Brain": "Brain",
    "Salivary gland": "Salivary glands",
    "Skin": "Skin",
}

# ============================================================
# READ LITERATURE CSV
# ============================================================

def read_literature(path):

    # The literature CSV has multiple header rows,
    # so read it without assigning a header.
    raw = pd.read_csv(path, header=None)

    records = []

    current_target = None

    # The actual data begins around row 5.
    for i in range(5, len(raw)):

        target = raw.iloc[i, 0]
        energy = raw.iloc[i, 1]

        # Skip rows without energy
        if pd.isna(energy):
            continue

        try:
            energy = float(energy)
        except (TypeError, ValueError):
            continue

        # The target name only appears on the first
        # row of each 3-energy group.
        if pd.notna(target):
            current_target = str(target).strip()

        # We only want 1 MeV
        if energy != ENERGY_MEV:
            continue

        # Literature CSV structure:
        #
        # Column 0 = Target organ
        # Column 1 = Energy
        #
        # Female:
        # Column 2 = ICRP 133
        # Column 3 = Geant4
        #
        # Male:
        # Column 7 = ICRP 133
        # Column 8 = Geant4

        female_geant4 = pd.to_numeric(
            raw.iloc[i, 3],
            errors="coerce"
        )

        male_geant4 = pd.to_numeric(
            raw.iloc[i, 8],
            errors="coerce"
        )

        records.append({
            "Target organ": current_target,
            "Literature Female Geant4": female_geant4,
            "Literature Male Geant4": male_geant4,
        })

    return pd.DataFrame(records)


# ============================================================
# READ YOUR GEANT4 SIMULATION CSV
# ============================================================

def read_simulation(path, phantom_name):

    df = pd.read_csv(path)

    # Keep only:
    #   Thyroid source
    #   1 MeV photon
    #   Correct phantom
    df = df[
        (
            df["Source Organ Name"]
            .astype(str)
            .str.strip()
            .str.lower()
            == SOURCE_ORGAN.lower()
        )
        &
        (
            pd.to_numeric(
                df["Source Energy (MeV)"],
                errors="coerce"
            )
            == ENERGY_MEV
        )
        &
        (
            df["Phantom"]
            .astype(str)
            .str.strip()
            == phantom_name
        )
    ].copy()

    df["Target Region Name"] = (
        df["Target Region Name"]
        .astype(str)
        .str.strip()
    )

    # Keep only the target regions we want
    df = df[
        df["Target Region Name"].isin(TARGETS.values())
    ].copy()

    return df[
        [
            "Target Region Name",
            "SAF (kg^-1)",
            "Statistical Uncertainty (%)"
        ]
    ]


# ============================================================
# MAIN COMPARISON
# ============================================================

def main():

    # --------------------------------------------------------
    # Read literature
    # --------------------------------------------------------

    literature = read_literature(LITERATURE_FILE)

    # --------------------------------------------------------
    # Read your Adult Male and Adult Female results
    # --------------------------------------------------------

    male = read_simulation(
        AM_FILE,
        "Adult Male"
    )

    female = read_simulation(
        AF_FILE,
        "Adult Female"
    )

    # Rename columns so we can distinguish AM and AF
    male = male.rename(
        columns={
            "SAF (kg^-1)": "AM SAF",
            "Statistical Uncertainty (%)":
                "AM Statistical Uncertainty (%)"
        }
    )

    female = female.rename(
        columns={
            "SAF (kg^-1)": "AF SAF",
            "Statistical Uncertainty (%)":
                "AF Statistical Uncertainty (%)"
        }
    )

    # --------------------------------------------------------
    # Build comparison table
    # --------------------------------------------------------

    rows = []

    for target, region in TARGETS.items():

        # ----------------------------------------------------
        # Literature value
        # ----------------------------------------------------

        lit_row = literature[
            literature["Target organ"]
            .astype(str)
            .str.strip()
            == target
        ]

        if lit_row.empty:
            raise ValueError(
                f"Literature target not found: {target}"
            )

        lit_row = lit_row.iloc[0]

        literature_male = float(
            lit_row["Literature Male Geant4"]
        )

        literature_female = float(
            lit_row["Literature Female Geant4"]
        )

        # ----------------------------------------------------
        # Adult Male Geant4 value
        # ----------------------------------------------------

        am_row = male[
            male["Target Region Name"]
            == region
        ]

        if am_row.empty:
            raise ValueError(
                f"Adult Male target region not found: {region}"
            )

        am_row = am_row.iloc[0]

        am_saf = float(
            am_row["AM SAF"]
        )

        am_uncertainty = float(
            am_row["AM Statistical Uncertainty (%)"]
        )

        # ----------------------------------------------------
        # Adult Female Geant4 value
        # ----------------------------------------------------

        af_row = female[
            female["Target Region Name"]
            == region
        ]

        if af_row.empty:
            raise ValueError(
                f"Adult Female target region not found: {region}"
            )

        af_row = af_row.iloc[0]

        af_saf = float(
            af_row["AF SAF"]
        )

        af_uncertainty = float(
            af_row["AF Statistical Uncertainty (%)"]
        )

        # ----------------------------------------------------
        # Calculate percentage differences
        #
        # Difference (%) =
        # (Your value - Literature value)
        # / Literature value * 100
        # ----------------------------------------------------

        male_difference = (
            (am_saf - literature_male)
            / literature_male
            * 100
        )

        female_difference = (
            (af_saf - literature_female)
            / literature_female
            * 100
        )

        # ----------------------------------------------------
        # Store result
        # ----------------------------------------------------

        rows.append({

            "Target organ": target,

            "Target region": region,

            "Literature Male Geant4 (kg^-1)":
                literature_male,

            "Adult Male SAF (kg^-1)":
                am_saf,

            "Male difference (%)":
                male_difference,

            "Literature Female Geant4 (kg^-1)":
                literature_female,

            "Adult Female SAF (kg^-1)":
                af_saf,

            "Female difference (%)":
                female_difference,

            "AM statistical uncertainty (%)":
                am_uncertainty,

            "AF statistical uncertainty (%)":
                af_uncertainty,
        })

    # Convert to DataFrame
    comparison = pd.DataFrame(rows)


    # ========================================================
    # SAVE COMPARISON TABLE
    # ========================================================

    output_csv = OUTPUT_COMPARISON_FILE

    comparison.to_csv(
        output_csv,
        index=False
    )


    # ========================================================
    # PRINT RESULTS
    # ========================================================

    print()
    print("=" * 170)
    print(
        "1 MeV LIVER-SOURCE SAF COMPARISON"
    )
    print("=" * 170)

    print()

    display_comparison = comparison[
        [
            "Target organ",

            "Literature Male Geant4 (kg^-1)",

            "Adult Male SAF (kg^-1)",

            "Male difference (%)",

            "Literature Female Geant4 (kg^-1)",

            "Adult Female SAF (kg^-1)",

            "Female difference (%)",
        ]
    ].copy()


    # Format SAF values in scientific notation
    display_comparison[
        "Literature Male Geant4 (kg^-1)"
    ] = display_comparison[
        "Literature Male Geant4 (kg^-1)"
    ].map(lambda x: f"{x:.2e}")

    display_comparison[
        "Adult Male SAF (kg^-1)"
    ] = display_comparison[
        "Adult Male SAF (kg^-1)"
    ].map(lambda x: f"{x:.2e}")

    display_comparison[
        "Literature Female Geant4 (kg^-1)"
    ] = display_comparison[
        "Literature Female Geant4 (kg^-1)"
    ].map(lambda x: f"{x:.2e}")

    display_comparison[
        "Adult Female SAF (kg^-1)"
    ] = display_comparison[
        "Adult Female SAF (kg^-1)"
    ].map(lambda x: f"{x:.2e}")


    # Format percentage differences as normal decimal
    display_comparison[
        "Male difference (%)"
    ] = display_comparison[
        "Male difference (%)"
    ].map(lambda x: f"{x:.2f}")

    display_comparison[
        "Female difference (%)"
    ] = display_comparison[
        "Female difference (%)"
    ].map(lambda x: f"{x:.2f}")


    # Print formatted table
    print(
        display_comparison.to_string(
            index=False
        )
    )

    print()

    print(f"Comparison table saved to:\n{output_csv}")

    # ========================================================
    # PLOT
    # ========================================================

    x = range(len(comparison))

    width = 0.22

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(14, 11),
        sharex=True
    )


    # ========================================================
    # ADULT MALE PLOT
    # ========================================================

    ax = axes[0]

    ax.bar(
        [
            i - width / 2
            for i in x
        ],

        comparison[
            "Literature Male Geant4 (kg^-1)"
        ],

        width=width,

        label="Literature Geant4"
    )

    ax.bar(
        [
            i + width / 2
            for i in x
        ],

        comparison[
            "Adult Male SAF (kg^-1)"
        ],

        width=width,

        label="This Geant4 simulation"
    )

    # SAF spans several orders of magnitude
    ax.set_yscale("log")

    ax.set_ylabel(
        "SAF (kg$^{-1}$)"
    )

    ax.set_title(
        "Adult Male — Liver source, 1 MeV photon"
    )

    ax.grid(
        axis="y",
        which="both",
        alpha=0.25
    )

    ax.legend()


    # ========================================================
    # ADULT FEMALE PLOT
    # ========================================================

    ax = axes[1]

    ax.bar(
        [
            i - width / 2
            for i in x
        ],

        comparison[
            "Literature Female Geant4 (kg^-1)"
        ],

        width=width,

        label="Literature Geant4"
    )

    ax.bar(
        [
            i + width / 2
            for i in x
        ],

        comparison[
            "Adult Female SAF (kg^-1)"
        ],

        width=width,

        label="This Geant4 simulation"
    )

    ax.set_yscale("log")

    ax.set_ylabel(
        "SAF (kg$^{-1}$)"
    )

    ax.set_title(
        "Adult Female — Liver source, 1 MeV photon"
    )

    ax.grid(
        axis="y",
        which="both",
        alpha=0.25
    )

    ax.legend()


    # ========================================================
    # X-AXIS
    # ========================================================

    axes[1].set_xticks(
        list(x)
    )

    axes[1].set_xticklabels(
        comparison["Target organ"],
        rotation=45,
        ha="right"
    )

    axes[1].set_xlabel(
        "Target organ"
    )


    # ========================================================
    # FIGURE SETTINGS
    # ========================================================

    fig.suptitle(
        "Liver-source SAF comparison at 1 MeV",
        fontsize=15
    )

    fig.tight_layout()


    # ========================================================
    # SAVE FIGURE
    # ========================================================

    output_png = OUTPUT_PNG

    fig.savefig(
        output_png,
        dpi=300,
        bbox_inches="tight"
    )

    print()
    print(
        f"Comparison plot saved to:\n{output_png}"
    )

    plt.show()


# ============================================================
# RUN PROGRAM
# ============================================================

if __name__ == "__main__":
    main()