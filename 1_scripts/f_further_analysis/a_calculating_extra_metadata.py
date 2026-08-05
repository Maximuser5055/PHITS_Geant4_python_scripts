#This script calculates additional metadata based on the extracted metadata from PHITS .out and .inp files. 
# It reads the extracted metadata, performs calculations, and saves the results to a CSV file.

# temp
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import necessary libraries
import pandas as pd
import b_config.a_config as config

def calculate_extra_metadata():
    if config.SIMULATION_CODE.upper() == "PHITS":
        input_prefix = "1_phits"
        output_prefix = "3_phits"

    elif config.SIMULATION_CODE.upper() == "GEANT4":
        input_prefix = "2_geant4"
        output_prefix = "4_geant4"

    else:
        raise ValueError(
            f"Unknown SIMULATION_CODE: {config.SIMULATION_CODE}"
        )

    # Input metadata log
    input_file = config.RESULTS_DIR / f"{input_prefix}_all_simulations_log.csv"

    # Output statistics
    output_file = config.RESULTS_DIR / f"{output_prefix}_extra_metadata_statistics.csv"

    # Output directory
    results_dir = input_file.parent

    df = pd.read_csv(input_file)

    # Format seconds to days, hours, minutes, and seconds
    def format_duration(seconds):
        if pd.isna(seconds):
            return ""

        seconds = int(round(seconds))

        days, rem = divmod(seconds, 86400)
        hours, rem = divmod(rem, 3600)
        minutes, seconds = divmod(rem, 60)

        return f"{days}d {hours:02}h {minutes:02}m {seconds:02}s"

    # Return summary statistics for a grouping.
    def summarize(df, group_cols, category):

        summary = (
            df.groupby(group_cols)["Individual Wall Time (s)"]
            .agg(
                Count="count",
                Average="mean",
                StdDev="std",
                Minimum="min",
                Maximum="max"
            )
            .reset_index()
            .round(2)
        )

        summary.rename(
            columns={
                "Average": "Average (s)",
                "StdDev": "StdDev (s)",
                "Minimum": "Minimum (s)",
                "Maximum": "Maximum (s)",
            },
            inplace=True,
        )

        summary.insert(0, "Category", category)

        # Human-readable time columns
        for col in ["Average (s)"]:
            summary[f"{col.replace('(s)', '(dhms)')}"] = (summary[col].apply(format_duration))

        # Ensure all grouping columns exist
        for col in ["Source Energy (MeV)", "Phantom", "Source Organ", "Source Type"]:
            if col not in summary.columns:
                summary[col] = ""

        # Consistent column order
        summary = summary[
            [
                "Category",
                "Source Type",
                "Source Energy (MeV)",
                "Phantom",
                "Source Organ",
                "Count",
                "Average (s)",
                "Average (dhms)",
                "StdDev (s)",
                "Minimum (s)",
                "Maximum (s)",
            ]
        ]

        return summary

    summary_energy = summarize(
        df,
        ["Source Type", "Source Energy (MeV)"],
        "Type + Energy"
    )

    summary_energy_phantom = summarize(
        df,
        ["Phantom", "Source Type", "Source Energy (MeV)"],
        "Phantom + Type + Energy"
    )

    summary_energy_organ = summarize(
        df,
        ["Source Organ", "Source Type", "Source Energy (MeV)"],
        "Source Organ + Type + Energy"
    )

    # Starting and termination datetime to calculate for the total elapsed time
    df["Starting Datetime"] = pd.to_datetime(df["Starting Datetime"])
    df["Termination Datetime"] = pd.to_datetime(df["Termination Datetime"])

    earliest_start = df["Starting Datetime"].min()
    latest_finish = df["Termination Datetime"].max()
    total_elapsed = latest_finish - earliest_start
    total_elapsed_seconds = total_elapsed.total_seconds()

    summary_all = pd.DataFrame({
        "Category": ["All simulations"],
        "Source Type": [""],
        "Source Energy (MeV)": [""],
        "Phantom": [""],
        "Source Organ": [""],
        "Count": [df["Individual Wall Time (s)"].count()],
        "Average (s)": [round(df["Individual Wall Time (s)"].mean(), 2)],
        "Average (dhms)": [format_duration(round(df["Individual Wall Time (s)"].mean(), 2))],
        "StdDev (s)": [round(df["Individual Wall Time (s)"].std(), 2)],
        "Minimum (s)": [round(df["Individual Wall Time (s)"].min(), 2)],
        "Maximum (s)": [round(df["Individual Wall Time (s)"].max(), 2)],
        "Overall Wall Time (s)": [total_elapsed_seconds],
        "Overall Wall Time (dhms)": [format_duration(total_elapsed_seconds)],
    })

    summary = pd.concat(
        [
            summary_energy,
            summary_energy_phantom,
            summary_energy_organ,
            summary_all
        ],
        ignore_index=True
    )

    summary.to_csv(output_file,index=False)

    print(f"Extra metadata statistics written to:\n"f"{output_file}")

calculate_metadata()
