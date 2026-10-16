#This script calculates additional metadata based on the extracted metadata from PHITS .out and .inp files. 
# It reads the extracted metadata, performs calculations, and saves the results to a CSV file.

# Import necessary libraries
import pandas as pd
import b_config.a_config as config

def calculate_extra_metadata(params):
    if params["simulation_code"].upper() == "PHITS":
        input_prefix = "a_phits"
        output_prefix = "b_phits"
        input_file = config.RESULTS_PHITS_DIR / f"{input_prefix}_all_simulations_log.csv"
        output_file = config.RESULTS_PHITS_DIR / f"{output_prefix}_extra_metadata_statistics.csv"
        
    elif params["simulation_code"].upper() == "GEANT4":
        input_prefix = "a_geant4"
        output_prefix = "b_geant4"
        input_file = config.RESULTS_GEANT4_DIR / f"{input_prefix}_all_simulations_log.csv"
        output_file = config.RESULTS_GEANT4_DIR / f"{output_prefix}_extra_metadata_statistics.csv"

    else:
        raise ValueError(
            f"Unknown SIMULATION_CODE: {params["simulation_code"].upper()}"
        )

    df = pd.read_csv(input_file)

    # Geant4 now records only total elapsed wall time. PHITS retains its
    # existing column name, so support both schemas.
    wall_time_candidates = [
        "Total Elapsed Wall Time (s)",
        "Individual Wall Time (s)",
    ]
    wall_time_column = next(
        (column for column in wall_time_candidates if column in df.columns),
        None,
    )
    if wall_time_column is None:
        raise KeyError(
            "No wall-time column found in metadata CSV. Expected one of: "
            + ", ".join(wall_time_candidates)
            + f". Available columns: {df.columns.tolist()}"
        )

    df[wall_time_column] = pd.to_numeric(df[wall_time_column], errors="coerce")

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
            df.groupby(group_cols)[wall_time_column]
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

    # With Geant4's new timing format, start/termination timestamps are no
    # longer stored. Report the sum of per-job wall times instead; this is
    # cumulative simulation time, not calendar elapsed time across parallel jobs.
    mean_wall_time = df[wall_time_column].mean()
    summary_all = pd.DataFrame({
        "Category": ["All simulations"],
        "Source Type": [""],
        "Source Energy (MeV)": [""],
        "Phantom": [""],
        "Source Organ": [""],
        "Count": [df[wall_time_column].count()],
        "Average (s)": [round(mean_wall_time, 2) if pd.notna(mean_wall_time) else float("nan")],
        "Average (dhms)": [format_duration(mean_wall_time)],
        "StdDev (s)": [round(df[wall_time_column].std(), 2)],
        "Minimum (s)": [round(df[wall_time_column].min(), 2)],
        "Maximum (s)": [round(df[wall_time_column].max(), 2)],
        "Total Simulation Wall Time (sum of jobs, s)": [df[wall_time_column].sum()],
        "Total Simulation Wall Time (sum of jobs, dhms)": [format_duration(df[wall_time_column].sum())],
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
