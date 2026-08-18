# This script checks the statistical uncertainty of each PHITS/Geant4 simulation.
# If any source-target organ pair has a statistical uncertainty >= 5%,
# the corresponding PHITS input file is updated with new maxcas/maxbch values
# provided by the user.

# Import necessary libraries
from pathlib import Path
import pandas as pd
import re
import b_config.a_config as config

def check_uncertainty(params):

    # -------------------------------------------------------------------------
    # Parameters
    # -------------------------------------------------------------------------

    uncertainty_limit = params["uncertainty_limit"]
    simulation_code = params["simulation_code"].upper()

    # -------------------------------------------------------------------------
    # Directories and configs
    # -------------------------------------------------------------------------

    results_phits_dir = config.RESULTS_PHITS_DIR
    results_geant4_dir = config.RESULTS_GEANT4_DIR
    
    if simulation_code == "PHITS":

        input_dir = config.GENERATED_INPUTS_DIR

        csv_files = [
            results_phits_dir / "f_phits_target_regions_dose_SAFs_AM.csv",
            results_phits_dir / "g_phits_target_regions_dose_SAFs_AF.csv",
        ]

        rerun_file = config.PHITS_RERUN_CSV_FILE

        input_extension = ".inp"

    elif simulation_code == "GEANT4":

        input_dir = config.GEANT4_GENERATED_INPUTS_DIR

        csv_files = [
            results_geant4_dir / "f_geant4_target_regions_dose_SAFs_AM.csv",
            results_geant4_dir / "g_geant4_target_regions_dose_SAFs_AF.csv",
        ]

        rerun_file = config.GEANT4_RERUN_CSV_FILE

        input_extension = ".in"

    else:
        raise ValueError(f"Unsupported simulation code: {simulation_code}")

    # -------------------------------------------------------------------------
    # Read CSV files
    # -------------------------------------------------------------------------

    dfs = []

    for csv_file in csv_files:

        if not csv_file.exists():

            print(
                f"\n[SKIP] Phantom result file not found:"
            )

            print(
                f"       {csv_file}"
            )

            continue

        print(
            f"\nReading: {csv_file.name}"
        )

        dfs.append(
            pd.read_csv(csv_file)
        )

    if not dfs:

        raise FileNotFoundError(
            "No phantom tally CSV files found."
        )

    df = pd.concat(
        dfs,
        ignore_index=True
    )

    # -------------------------------------------------------------------------
    # Ignore rows with zero dose
    # -------------------------------------------------------------------------

    df = df[df["Dose (Gy/source)"] > 0]

    # -------------------------------------------------------------------------
    # Determine which simulations require rerunning
    # -------------------------------------------------------------------------

    group_columns = [
        "Phantom",
        "Source Organ ID",
        "Source Organ Name",
        "Source Type",
        "Source Energy (MeV)"
    ]

    summary = (
        df.groupby(group_columns)["Statistical Uncertainty (%)"]
        .max()
        .reset_index()
        .rename(columns={"Statistical Uncertainty (%)":
                        "Maximum Statistical Uncertainty (%)"})
    )

    rerun = summary[summary["Maximum Statistical Uncertainty (%)"] >= uncertainty_limit].copy()

    rerun = rerun.sort_values(
        by=[
            "Phantom",
            "Source Organ ID",
            "Source Type",
            "Source Energy (MeV)",
        ]
    )

    # -------------------------------------------------------------------------
    # Determine corresponding input files and job directories
    # -------------------------------------------------------------------------

    rerun["Input File"] = ""
    rerun["Job Directory"] = ""

    suffix = re.escape(input_extension)

    for index, row in rerun.iterrows():

        phantom = "AM" if row["Phantom"] == "Adult Male" else "AF"

        source = row["Source Organ Name"]

        particle = row["Source Type"]

        if simulation_code == "GEANT4":
            particle = next(
                k for k, v in config.GEANT4_SOURCE_TYPE_MAP.items()
                if v == particle
            )

        energy = row["Source Energy (MeV)"]

        input_pattern = (f"*MRCP_{phantom}_source_{source}_{particle}_energy_*{input_extension}")

        for file in input_dir.rglob(input_pattern):

            match = re.search(rf"_energy_([0-9Ee.+-]+){suffix}$",file.name,)

            if not match:
                continue

            file_energy = float(match.group(1))

            if abs(file_energy - energy) < 1e-12:

                rerun.at[index, "Input File"] = str(
                    file.relative_to(input_dir)
                )

                rerun.at[index, "Job Directory"] = str(
                    file.parent.relative_to(input_dir)
                )

                break

    # -------------------------------------------------------------------------
    # Save rerun report
    # -------------------------------------------------------------------------

    rerun.to_csv(rerun_file,index=False)
    
    # -------------------------------------------------------------------------
    # Report
    # -------------------------------------------------------------------------

    print("\n========================================")
    print("Uncertainty Check")
    print("========================================")

    if rerun.empty:

        print(f"\nAll simulations satisfy the {uncertainty_limit:.1f}% criterion.")

        raise SystemExit

    # -------------------------------------------------------------------------
    # Add PASS / FAIL column
    # -------------------------------------------------------------------------

    summary["Status"] = summary[
        "Maximum Statistical Uncertainty (%)"
    ].apply(
        lambda x: "PASS" if x < uncertainty_limit else "FAIL"
    )

    print(f"\nThreshold         : {uncertainty_limit:.1f}%")
    print(f"Total simulations : {len(summary)}")
    print(f"PASS              : {(summary['Status'] == 'PASS').sum()}")
    print(f"FAIL              : {(summary['Status'] == 'FAIL').sum()}\n")

    display = (
        summary.rename(
            columns={
                "Source Organ Name": "Source Organ",
                "Source Type": "Source Type",
                "Source Energy (MeV)": "Energy (MeV)",
                "Maximum Statistical Uncertainty (%)": "Max Uncertainty (%)"
            }
        )
        .sort_values(
            by=[
                "Phantom",
                "Source Organ ID",
                "Source Type",
                "Energy (MeV)",
            ]
        )
    )

    print(
        display.to_string(
            index=False,
            justify="left",
            col_space={
                "Phantom": 13,
                "Source Organ": 15,
                "Source Type": 12,
                "Energy (MeV)": 12,
                "Max Uncertainty (%)": 20,
                "Status": 8,
            },
            formatters={
                "Energy (MeV)": "{:.3f}".format,
                "Max Uncertainty (%)": "{:.2f}".format,
            },
        )
    )

    # -------------------------------------------------------------------------
    # Ask for confirmation in updating input files
    # -------------------------------------------------------------------------

    answer = input(
        "\nUpdate the input files for all FAILED simulations? (y/n): "
    ).strip().lower()

    if answer not in ("y", "yes"):
        print("\nNo input files were modified.")
        raise SystemExit

    # -------------------------------------------------------------------------
    # Ask for new simulation parameters
    # -------------------------------------------------------------------------

    if simulation_code == "PHITS":

        print("\nEnter new PHITS parameter/s:")

        new_maxcas = input(f"New maxbch (no. of batches) [Current = {params["maxcas"]}]: ").strip()
        new_maxcas = int(new_maxcas) if new_maxcas else params["maxcas"]
        config.update_config("MAXCAS", new_maxcas)
        config.MAXCAS = new_maxcas

        new_maxbch = input(f"New maxbch (no. of batches) [Current = {params["maxbch"]}]: ").strip()
        new_maxbch = int(new_maxbch) if new_maxbch else params["maxbch"]
        config.update_config("MAXBCH", new_maxbch)
        config.MAXBCH = new_maxbch

        def update_input(text):

            text, n1 = re.subn(
                r"maxcas\s*=\s*\d+",
                f"maxcas = {new_maxcas}",
                text,
            )

            text, n2 = re.subn(
                r"maxbch\s*=\s*\d+",
                f"maxbch = {new_maxbch}",
                text,
            )

            return text, (n1 > 0 and n2 > 0)

    elif simulation_code == "GEANT4":

        print("\nEnter new Geant4 parameter/s:")

        new_nps = input(f"New nps (no. of particle histories) [Current = {params["nps"]}]: ").strip()
        new_nps = int(new_nps) if new_nps else params["nps"]
        config.update_config("NPS", new_nps)
        config.NPS = new_nps
        
        def update_input(text):

            text, n1 = re.subn(
                r"/run/beamOn\s+\d+",
                f"/run/beamOn {new_nps}",
                text,
            )

            return text, (n1 > 0)

    else:
        raise ValueError(f"Unsupported simulation code: {simulation_code}")

    # -------------------------------------------------------------------------
    # Update corresponding input files
    # -------------------------------------------------------------------------

    updated = 0

    for _, row in rerun.iterrows():

        if not row["Input File"]:
            print(f"No matching input file found for:")
            print(row)
            continue

        input_file = input_dir / Path(row["Input File"])

        if not input_file.exists():
            print(f"Missing input file: {input_file}")
            continue

        text = input_file.read_text(encoding="utf-8")

        text, success = update_input(text)

        if not success:
            print(f"Could not update {input_file.name}")
            continue

        print(f"Updating {input_file.name}")

        input_file.write_text(text, encoding="utf-8")

        updated += 1

    print("\n----------------------------------------")
    print(f"{updated} input file(s) updated.")
    print("----------------------------------------")

    print("\nYou can now rerun the failed simulations.")