#This script calculates additional metadata based on the extracted metadata from PHITS .out and .inp files. 
# It reads the extracted metadata, performs calculations, and saves the results to a CSV file.

# Import necessary libraries
import pandas as pd
import c_config as config

def calculate_metadata():
    # Master log produced by e_extracting_metadata.py
    log_file = config.RESULTS_DIR / "1_all_simulations_log.csv"

    # Output directory
    results_dir = log_file.parent

    df = pd.read_csv(log_file)

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
            df.groupby(group_cols)["Elapsed Time (s)"]
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

        summary.insert(0, "Category", category)

        # Human-readable time columns
        for col in ["Average"]:
            summary[f"{col} (dhms)"] = summary[col].apply(format_duration)

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
                "Average",
                "Average (dhms)",
                "StdDev",
                "Minimum",
                "Maximum",
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
        "Count": [df["Elapsed Time (s)"].count()],
        "Average": [round(df["Elapsed Time (s)"].mean(), 2)],
        "Average (dhms)": [format_duration(round(df["Elapsed Time (s)"].mean(), 2))],
        "StdDev": [round(df["Elapsed Time (s)"].std(), 2)],
        "Minimum": [round(df["Elapsed Time (s)"].min(), 2)],
        "Maximum": [round(df["Elapsed Time (s)"].max(), 2)],
        "Total Elapsed Time (s)": [total_elapsed_seconds],
        "Total Elapsed Time (dhms)": [format_duration(total_elapsed_seconds)],
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

    summary.to_csv(
        results_dir / "2_extra_metadata_statistics.csv",
        index=False
    )

    print("Extra metadata statistics written successfully.")