# This script combines the target organs defined in the mapping file and calculates the true dose and SAFs for each target region.

# Import necessary libraries
import numpy as np
import pandas as pd

from b_config import a_config as config


def combine_target_organs_and_calculate_true_dose_and_SAFs(params):

    print("\nConfiguring target regions...")

    # ==========================================================
    # Determine input/output files
    # ==========================================================

    simulation = params["simulation_code"].upper()
    mapping = pd.read_csv(config.TARGET_REGION_CSV)
    skeletal = pd.read_csv(config.SKELETAL_MASSES_CSV)
    results_phits_dir = config.RESULTS_PHITS_DIR
    results_geant4_dir = config.RESULTS_GEANT4_DIR
    
    RBM_REGION = "Red (active) marrow"
    ENDOSTEUM_REGION = "50-um endosteal region"

    skeletal = skeletal.set_index("Organ ID")

    am_marrow = skeletal["Ref_AM_Marrow_Mass(g)"]
    af_marrow = skeletal["Ref_AF_Marrow_Mass(g)"]

    am_endosteum = skeletal["Ref_AM_Endosteum_Mass(g)"]
    af_endosteum = skeletal["Ref_AF_Endosteum_Mass(g)"]

    # Process both phantoms
    for phantom in ["AM", "AF"]:
        
        if simulation == "PHITS":
            input_prefix = "c" if phantom == "AM" else "d"
            output_prefix = "f" if phantom == "AM" else "g"
            fluence_csv = results_phits_dir / "e_phits_rbm_endosteum_icrp116.csv"
            input_csv = results_phits_dir /f"{input_prefix}_{simulation.lower()}_all_dose_and_SAFs_{phantom}.csv"
            output_csv = results_phits_dir / f"{output_prefix}_{simulation.lower()}_target_regions_dose_SAFs_{phantom}.csv"

        elif simulation == "GEANT4":
            input_prefix = "c" if phantom == "AM" else "d"
            output_prefix = "f" if phantom == "AM" else "g"
            fluence_csv = results_geant4_dir / "e_geant4_rbm_endosteum_icrp116.csv"
            input_csv = results_geant4_dir /f"{input_prefix}_{simulation.lower()}_all_dose_and_SAFs_{phantom}.csv"
            output_csv = results_geant4_dir / f"{output_prefix}_{simulation.lower()}_target_regions_dose_SAFs_{phantom}.csv"

        else:
            raise ValueError(f"Unsupported simulation code: {simulation}")

        # ==========================================================
        # Skip phantom if input file does not exist
        # ==========================================================

        if not input_csv.is_file():

            print(f"\n[SKIP] {phantom} input file not found:")
            print(f"       {input_csv}")
            print(f"       Skipping {phantom}...")

            continue

        print(f"\nProcessing {phantom}...")

        df = pd.read_csv(input_csv)

        df["Target Organ ID"] = df["Target Organ ID"].astype(int)

        fluence_df = pd.read_csv(fluence_csv) if fluence_csv.is_file() else pd.DataFrame()
        
        output_rows = []

        grouping_columns = [
            "Phantom",
            "Source Organ ID",
            "Source Organ Name",
            "Source Type",
            "Source Energy (MeV)"
        ]

        grouped = df.groupby(grouping_columns)

        if phantom == "AM":
            marrow_lookup = am_marrow
            endosteum_lookup = am_endosteum
        else:
            marrow_lookup = af_marrow
            endosteum_lookup = af_endosteum

        # ==========================================================
        # Loop over every source organ / energy
        # ==========================================================

        for group_key, source_df in grouped:

            # ------------------------------------------------------
            # Process every target region listed in mapping file
            # ------------------------------------------------------

            for _, region in mapping.iterrows():

                target_name = region["Target region"]

                ids = [int(x.strip()) for x in str(region["ID number(s)"]).split("_")]

                rows = source_df[
                    source_df["Target Organ ID"].isin(ids)
                ].copy()

                if rows.empty:
                    continue

                expected_ids = set(ids)
                found_ids = set(rows["Target Organ ID"])

                missing_ids = sorted(expected_ids - found_ids)

                if missing_ids:
                    print(
                        f"Warning: {target_name} "
                        f"(Source Organ ID {group_key[1]}, "
                        f"Energy {group_key[4]} MeV) "
                        f"is missing Target Organ IDs: {missing_ids}"
                    )

                # ======================================================
                # Determine masses used for weighting
                # ======================================================

                if target_name == RBM_REGION:

                    mapped_masses = rows["Target Organ ID"].map(marrow_lookup)

                    if mapped_masses.isna().any():
                        missing_ids = rows.loc[
                            mapped_masses.isna(),
                            "Target Organ ID"
                        ].tolist()

                        raise ValueError(
                            f"Missing marrow masses for Organ IDs: {missing_ids}"
                        )

                    masses = mapped_masses.to_numpy(float)

                elif target_name == ENDOSTEUM_REGION:

                    mapped_masses = rows["Target Organ ID"].map(endosteum_lookup)

                    if mapped_masses.isna().any():
                        missing_ids = rows.loc[
                            mapped_masses.isna(),
                            "Target Organ ID"
                        ].tolist()

                        raise ValueError(
                            f"Missing endosteum masses for Organ IDs: {missing_ids}"
                        )

                    masses = mapped_masses.to_numpy(float)

                else:

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

                row = {

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

                    "Target Organ IDs":
                        "_".join(map(str, ids)),

                    "Target Region Name":
                        target_name,

                    "Target Region Mass (g)":
                        total_mass,

                    "Dose (Gy/source)":
                        combined_dose,

                    "SAF (kg^-1)":
                        combined_saf,

                    "Relative Error":
                        combined_relative_error,

                    "Statistical Uncertainty (%)":
                        statistical_uncertainty,
                }

                # Add calculation method ONLY for RBM and endosteum
                if target_name in {
                    RBM_REGION,
                    ENDOSTEUM_REGION,
                }:
                    row["Calculation Method"] = (
                        "Direct dose calculation"
                    )

                output_rows.append(row)

            # ==========================================================
            # Append fluence-to-dose response-function results
            # ==========================================================

            if not fluence_df.empty:

                # Select the same source simulation
                # as the current direct-dose group.
                source_fluence = fluence_df[
                    (fluence_df["Phantom"] == group_key[0])
                    &
                    (fluence_df["Source Organ"] == group_key[2]
                    )
                    &
                    (fluence_df["Source Type"] == group_key[3]
                    )
                    &
                    (fluence_df["Source Energy (MeV)"] == group_key[4]
                    )
                ]

                # ------------------------------------------------------
                # Keep the complete skeletal-region table.
                # ------------------------------------------------------

                if not source_fluence.empty:

                    fluence_result = source_fluence.iloc[0]

                    # --------------------------------------------------
                    # RBM
                    # --------------------------------------------------

                    if pd.notna(fluence_result["RBM Dose (Gy/source)"]):

                        output_rows.append({

                            "Phantom":
                                fluence_result["Phantom"],

                            "Source Organ ID":
                                group_key[1],

                            "Source Organ Name":
                                fluence_result["Source Organ"],

                            "Source Type":
                                fluence_result["Source Type"],

                            "Source Energy (MeV)":
                                fluence_result[
                                    "Source Energy (MeV)"
                                ],

                            "Target Organ IDs":
                                "_".join(
                                    source_fluence.loc[
                                        source_fluence["Marrow mass (kg)"].notna(),
                                        "Organ ID"
                                    ]
                                    .astype(int)
                                    .astype(str)
                                    .tolist()
                                ),

                            "Target Region Name":
                                RBM_REGION,

                            "Target Region Mass (g)":
                                fluence_result[
                                    "Total Active Marrow Mass (kg)"
                                ] * 1000.0,

                            "Dose (Gy/source)":
                                fluence_result[
                                    "RBM Dose (Gy/source)"
                                ],

                            "SAF (kg^-1)":
                                fluence_result[
                                    "RBM SAF (kg^-1)"
                                ],

                            "Relative Error":
                                fluence_result[
                                    "RBM Relative Error"
                                ],

                            "Statistical Uncertainty (%)":
                                fluence_result[
                                    "RBM Statistical Uncertainty (%)"
                                ],

                            "Calculation Method":
                                "Fluence-to-dose response functions",
                        })

                    # --------------------------------------------------
                    # Endosteum
                    # --------------------------------------------------

                    if pd.notna(fluence_result["Endosteum Dose (Gy/source)"]):

                        output_rows.append({

                            "Phantom":
                                fluence_result["Phantom"],

                            "Source Organ ID":
                                group_key[1],

                            "Source Organ Name":
                                fluence_result["Source Organ"],

                            "Source Type":
                                fluence_result["Source Type"],

                            "Source Energy (MeV)":
                                fluence_result[
                                    "Source Energy (MeV)"
                                ],

                            "Target Organ IDs":
                                "_".join(
                                    source_fluence.loc[
                                        source_fluence["Endosteum mass (kg)"].notna(),
                                        "Organ ID"
                                    ]
                                    .astype(int)
                                    .astype(str)
                                    .tolist()
                                ),

                            "Target Region Name":
                                ENDOSTEUM_REGION,

                            "Target Region Mass (g)":
                                fluence_result[
                                    "Total Endosteum Mass (kg)"
                                ] * 1000.0,

                            "Dose (Gy/source)":
                                fluence_result[
                                    "Endosteum Dose (Gy/source)"
                                ],

                            "SAF (kg^-1)":
                                fluence_result[
                                    "Endosteum SAF (kg^-1)"
                                ],

                            "Relative Error":
                                fluence_result[
                                    "Endosteum Relative Error"
                                ],

                            "Statistical Uncertainty (%)":
                                fluence_result[
                                    "Endosteum Statistical Uncertainty (%)"
                                ],

                            "Calculation Method":
                                "Fluence-to-dose response functions",
                        })

        # ==========================================================
        # Save output
        # ==========================================================
        if not output_rows:
            print(f"\n[WARNING] No target-region data generated for {phantom}.")
            continue

        output = pd.DataFrame(output_rows)

        # ==========================================================
        # Sort according to the order in target_regions_ICRP_145.csv
        # ==========================================================

        target_region_order = mapping["Target region"].tolist()

        output["Target Region Name"] = pd.Categorical(
            output["Target Region Name"],
            categories=target_region_order,
            ordered=True,
        )

        output["Calculation Method"] = pd.Categorical(
            output["Calculation Method"],
            categories=[
                "Direct dose calculation",
                "Fluence-to-dose response functions",
            ],
            ordered=True,
        )

        output.sort_values(
            [
                "Phantom",
                "Source Organ ID",
                "Source Type",
                "Source Energy (MeV)",
                "Target Region Name",
            ],
            inplace=True,
            ignore_index=True,
        )

        output.to_csv(output_csv, index=False)

        print(f"Target regions configured successfully.\n"f"Saved to:\n{output_csv}")