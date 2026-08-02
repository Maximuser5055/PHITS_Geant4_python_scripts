# This script checks the statistical uncertainty of each PHITS simulation.
# If any source-target organ pair has a statistical uncertainty >= 5%,
# the corresponding PHITS input file is updated with new maxcas/maxbch values
# provided by the user.

# Import necessary libraries
from pathlib import Path
import pandas as pd
import re
import c_config as config

def check_uncertainty():
    # -------------------------------------------------------------------------
    # Directories
    # -------------------------------------------------------------------------

    results_dir = config.RESULTS_DIR

    input_dir = config.GENERATED_INPUTS_DIR

    csv_files = [
        results_dir / "3_all_dose_and_SAFs_AM.csv",
        results_dir / "4_all_dose_and_SAFs_AF.csv",
    ]

    # -------------------------------------------------------------------------
    # Read CSV files
    # -------------------------------------------------------------------------

    dfs = []

    for csv_file in csv_files:

        if csv_file.exists():
            dfs.append(pd.read_csv(csv_file))

    if not dfs:
        raise FileNotFoundError("No tally CSV files found.")

    df = pd.concat(dfs, ignore_index=True)

    # -------------------------------------------------------------------------
    # Ignore rows with zero dose
    # -------------------------------------------------------------------------

    df = df[df["Dose (Gy/source)"] > 0]

    # -------------------------------------------------------------------------
    # Determine which simulations require rerunning
    # -------------------------------------------------------------------------

    group_columns = [
        "Phantom",
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

    rerun = summary[
        summary["Maximum Statistical Uncertainty (%)"] >= config.UNCERTAINTY_LIMIT
    ].copy()

    rerun = rerun.sort_values(
        by=[
            "Phantom",
            "Source Type",
            "Source Energy (MeV)",
            "Source Organ Name",
        ]
    )

    # -------------------------------------------------------------------------
    # Determine corresponding input files and job directories
    # -------------------------------------------------------------------------

    rerun["Input File"] = ""
    rerun["Job Directory"] = ""

    for index, row in rerun.iterrows():

        phantom = "AM" if row["Phantom"] == "Adult Male" else "AF"

        source = row["Source Organ Name"]

        particle = row["Source Type"]

        energy = row["Source Energy (MeV)"]

        input_pattern = (
            f"MRCP_{phantom}_source_{source}_{particle}_energy_*.inp"
        )

        for file in input_dir.rglob(input_pattern):

            match = re.search(
                r"_energy_([0-9Ee.+-]+)\.inp$",
                file.name
            )

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

    rerun.to_csv(
        results_dir / "5_rerun_required.csv",
        index=False
    )

    # -------------------------------------------------------------------------
    # Report
    # -------------------------------------------------------------------------

    print("\n========================================")
    print("Uncertainty Check")
    print("========================================")

    if rerun.empty:

        print(f"\nAll simulations satisfy the {config.UNCERTAINTY_LIMIT:.1f}% criterion.")

        raise SystemExit

    # -------------------------------------------------------------------------
    # Add PASS / FAIL column
    # -------------------------------------------------------------------------

    summary["Status"] = summary[
        "Maximum Statistical Uncertainty (%)"
    ].apply(
        lambda x: "PASS" if x < config.UNCERTAINTY_LIMIT else "FAIL"
    )

    print(f"\nThreshold         : {config.UNCERTAINTY_LIMIT:.1f}%")
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
                "Source Type",
                "Energy (MeV)",
                "Source Organ",
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
    # Ask for new maxcas/maxbch
    # -------------------------------------------------------------------------

    print("\nEnter the new PHITS parameters.\n")

    new_maxcas = int(input("New maxcas : "))
    new_maxbch = int(input("New maxbch : "))

    # -------------------------------------------------------------------------
    # Update corresponding input files
    # -------------------------------------------------------------------------

    updated = 0

    for _, row in rerun.iterrows():

        input_file = input_dir / Path(row["Input File"])

        if not input_file.exists():
            print(f"Missing input file: {input_file}")
            continue

        text = input_file.read_text(encoding="utf-8")

        text, n1 = re.subn(
            r"maxcas\s*=\s*\d+",
            f"maxcas = {new_maxcas}",
        text
        )

        text, n2 = re.subn(
            r"maxbch\s*=\s*\d+",
            f"maxbch = {new_maxbch}",
            text
        )

        if n1 == 0 or n2 == 0:
            print(f"Could not update {input_file.name}")
            continue

        print(f"Updating {input_file.name}")

        input_file.write_text(text, encoding="utf-8")

        updated += 1

    print("\n----------------------------------------")
    print(f"{updated} input file(s) updated.")
    print("----------------------------------------")

    print("\nYou can now rerun those simulations using d_running_inputs.py.")