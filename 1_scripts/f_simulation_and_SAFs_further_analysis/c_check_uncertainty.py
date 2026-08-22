# This script checks the statistical uncertainty of each PHITS/Geant4 simulation.
# If any source-target organ pair has a statistical uncertainty >= 5%,
# the corresponding PHITS input file is updated with new maxcas/maxbch values
# provided by the user.

# Import necessary libraries
from pathlib import Path
import pandas as pd
import re
import b_config.a_config as config

# ============================================================
# EXISTING PUBLISHABLE SAF DATABASE CHECK
# ============================================================

PHANTOM_FILE_NAMES = {
    "MRCP_AF": "mrcp_af",
    "MRCP_AM": "mrcp_am",
    "MFCP_AF": "mfcp_af",
    "MFCP_AM": "mfcp_am",
}

def get_required_saf_database_files(phantom, simulation_code):
    """
    Return all required publishable SAF and STD files for
    the selected phantom.

    Both photon and electron databases are always required.
    """

    simulation_code = simulation_code.upper()

    if simulation_code == "PHITS":

        publishable_dir = config.RESULTS_PHITS_PUBLISHABLE_SAF_DATABASE_DIR

    elif simulation_code == "GEANT4":

        publishable_dir = config.RESULTS_GEANT4_PUBLISHABLE_SAF_DATABASE_DIR

    else:

        raise ValueError(f"Unsupported simulation code: {simulation_code}")
    
    phantom_map = {
        "MRCP_AF": ["MRCP_AF"],
        "MRCP_AM": ["MRCP_AM"],
        "MFCP_AF": ["MFCP_AF"],
        "MFCP_AM": ["MFCP_AM"],
        "MRCP_AF_AM": ["MRCP_AF", "MRCP_AM"],
        "MFCP_AF_AM": ["MFCP_AF", "MFCP_AM"],
    }

    if phantom not in phantom_map:
        raise ValueError(
            f"Unsupported phantom selection: {phantom}"
        )

    required_files = []

    for phantom_code in phantom_map[phantom]:

        phantom_name = PHANTOM_FILE_NAMES[phantom_code]

        # ----------------------------------------------------
        # Photon SAF database
        # ----------------------------------------------------

        required_files.extend([
            publishable_dir
            / f"{phantom_name}_photons_saf.csv",

            publishable_dir
            / f"{phantom_name}_photons_std.csv",
        ])

        # ----------------------------------------------------
        # Electron SAF database
        # ----------------------------------------------------

        required_files.extend([
            publishable_dir
            / f"{phantom_name}_electrons_saf.csv",

            publishable_dir
            / f"{phantom_name}_electrons_std.csv",
        ])

    return required_files

def check_existing_saf_database(phantom, simulation_code, uncertainty_limit,):
    """
    Check whether the complete publishable SAF database required
    by the selected phantom and simulation code exists.

    Both photon and electron SAF databases are required.

    The uncertainty is read directly from the corresponding
    publishable *_std.csv files.

    This function is report-only. It does not modify input files
    and does not request a rerun. The user interface in
    b_input_user_parameters.py decides whether to reuse the
    existing database or redo the SAF calculations.
    """

    required_files = get_required_saf_database_files(phantom, simulation_code)

    existing_files = [
        file
        for file in required_files
        if file.is_file()
    ]

    missing_files = [
        file
        for file in required_files
        if not file.is_file()
    ]

    if missing_files:

        return {
            "exists": bool(existing_files),
            "complete": False,
            "required_files": required_files,
            "existing_files": existing_files,
            "missing_files": missing_files,
            "max_uncertainty": None,
            "uncertainty_by_file": {},
            "uncertainty_pass": False,
        }

    uncertainty_by_file = {}
    overall_max_uncertainty = None

    for std_file in existing_files:

        dataframe = pd.read_csv(
            std_file,
            comment="#",
        )

        if dataframe.shape[1] < 3:
            raise ValueError(
                f"Uncertainty database {std_file.name} "
                "does not contain source-energy columns."
            )

        energy_columns = list(dataframe.columns[2:])

        numeric_values = (
            dataframe[energy_columns]
            .apply(pd.to_numeric, errors="coerce")
            .stack()
            .dropna()
        )

        if numeric_values.empty:
            raise ValueError(
                f"No numerical uncertainty values were found "
                f"in {std_file.name}."
            )

        file_max = float(numeric_values.max())

        uncertainty_by_file[std_file.name] = file_max

        if (
            overall_max_uncertainty is None
            or file_max > overall_max_uncertainty
        ):
            overall_max_uncertainty = file_max

    return {
        "exists": True,
        "complete": True,
        "required_files": required_files,
        "existing_files": existing_files,
        "missing_files": [],
        "max_uncertainty": overall_max_uncertainty,
        "uncertainty_by_file": uncertainty_by_file,
        "uncertainty_pass": (
            overall_max_uncertainty < uncertainty_limit
        ),
    }


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

        return

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
        return

    # -------------------------------------------------------------------------
    # Ask for new simulation parameters
    # -------------------------------------------------------------------------

    if simulation_code == "PHITS":

        print("\nEnter new PHITS parameter/s:")

        new_maxcas = input(f"New maxcas (no. of histories per batch) [Current = {params["maxcas"]}]: ").strip()
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