# This script combines the target organs defined in the mapping file and calculates the true dose and SAFs for each target region.

# temporary
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import necessary libraries
import numpy as np
import pandas as pd

from b_config import a_config as config


def combine_target_organs_and_calculate_true_dose_and_SAFs():

    print("\nConfiguring target regions...")

    # ==========================================================
    # Determine input/output files
    # ==========================================================

    simulation = config.SIMULATION_CODE

    # Process both phantoms
    for phantom in ["AM", "AF"]:
        
        if simulation == "PHITS":
            input_prefix = "e" if phantom == "AM" else "f"
            output_prefix = "i" if phantom == "AM" else "j"

        elif simulation == "GEANT4":
            input_prefix = "g" if phantom == "AM" else "h"
            output_prefix = "k" if phantom == "AM" else "l"

        else:
            raise ValueError(
                f"Unsupported simulation code: {config.SIMULATION_CODE}"
            )

        # ==========================================================
        # Read files
        # ==========================================================
        input_csv = (
            config.RESULTS_DIR /
            f"{input_prefix}_{simulation.lower()}_all_dose_and_SAFs_{phantom}.csv"
        )

        output_csv = (
            config.RESULTS_DIR /
            f"{output_prefix}_{simulation.lower()}_target_regions_dose_SAFs_{phantom}.csv"
        )

        print(f"\nProcessing {phantom}...")
        df = pd.read_csv(input_csv)
        mapping = pd.read_csv(config.TARGET_REGION_CSV)

        df["Target Organ ID"] = df["Target Organ ID"].astype(int)

        output_rows = []

        grouping_columns = [
            "Phantom",
            "Source Organ ID",
            "Source Organ Name",
            "Source Type",
            "Source Energy (MeV)"
        ]

        grouped = df.groupby(grouping_columns)

        # ==========================================================
        # Loop over every source organ / energy
        # ==========================================================

        for group_key, source_df in grouped:

            # ------------------------------------------------------
            # Process every target region listed in mapping file
            # ------------------------------------------------------

            for _, region in mapping.iterrows():

                target_name = region["Target region"]

                ids = [
                    int(x)
                    for x in str(region["ID number(s)"]).split("_")
                ]

                rows = source_df[
                    source_df["Target Organ ID"].isin(ids)
                ].copy()

                if rows.empty:
                    continue

                masses = rows["Target Organ Mass (g)"].to_numpy(float)
                doses = rows["Dose (Gy/source)"].to_numpy(float)
                safs = rows["SAF (kg^-1)"].to_numpy(float)
                rel_errors = rows["Relative Error"].to_numpy(float)

                total_mass = masses.sum()

                if total_mass == 0:
                    continue

                # ==================================================
                # Combined dose
                # ==================================================

                combined_dose = np.sum(
                    masses * doses
                ) / total_mass

                # ==================================================
                # Combined SAF
                # ==================================================

                combined_saf = np.sum(
                    masses * safs
                ) / total_mass

                # ==================================================
                # Combined uncertainty
                # ==================================================

                absolute_sigma = doses * rel_errors

                combined_sigma = (
                    np.sqrt(
                        np.sum(
                            (masses * absolute_sigma) ** 2
                        )
                    )
                    / total_mass
                )

                if combined_dose > 0:

                    combined_relative_error = (
                        combined_sigma /
                        combined_dose
                    )

                else:

                    combined_relative_error = 0.0

                statistical_uncertainty = (
                    combined_relative_error * 100
                )

                output_rows.append({

                    "Phantom":
                        group_key[0],

                    "Source Organ ID":
                        group_key[1],

                    "Source Organ Name":
                        group_key[2],

                    "Source Type":
                        group_key[3],

                    "Source Energy (MeV)":
                        group_key[4],

                    "Target Organ ID":
                        "_".join(map(str, ids)),

                    "Target Organ Name":
                        target_name,

                    "Target Organ Mass (g)":
                        total_mass,

                    "Dose (Gy/source)":
                        combined_dose,

                    "SAF (kg^-1)":
                        combined_saf,

                    "Relative Error":
                        combined_relative_error,

                    "Statistical Uncertainty (%)":
                        statistical_uncertainty

                })

        # ==========================================================
        # Save output
        # ==========================================================

        output = pd.DataFrame(output_rows)

        # ==========================================================
        # Sort according to the order in target_regions_ICRP_145.csv
        # ==========================================================

        target_region_order = mapping["Target region"].tolist()

        output["Target Organ Name"] = pd.Categorical(
            output["Target Organ Name"],
            categories=target_region_order,
            ordered=True,
        )

        output.sort_values(
            [
                "Source Organ ID",
                "Source Energy (MeV)",
                "Target Organ Name",
            ],
            inplace=True,
            ignore_index=True,
        )

        output.to_csv(output_csv, index=False)

        print(f"Target regions configured successfully.\n"f"Saved to:\n{output_csv}")

combine_target_organs_and_calculate_true_dose_and_SAFs()