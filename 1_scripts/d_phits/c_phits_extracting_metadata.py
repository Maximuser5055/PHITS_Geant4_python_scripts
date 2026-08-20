# This script extracts useful metadata from PHITS .out and .inp files using keyword matching.

# import necessary libraries
import re
from pathlib import Path
import csv
import b_config.a_config as config

def phits_extract_metadata_stats():
    # PHITS output file
    root = config.GENERATED_INPUTS_DIR 
    output_files = [file for file in sorted(root.rglob("phits_*.out"))
                    if "_deposit_" not in file.name
                    and "_fluence_" not in file.name
                    ]

    # Define metadata output file
    metadata_output_file = config.RESULTS_PHITS_DIR / "a_phits_all_simulations_log.csv"

    # Configs
    phantom_names = config.PHANTOM_NAMES

    # Regex patterns
    patterns = {

        "starting_date":
            re.compile(r"Starting Date\s*=\s*(.+)"),

        "starting_time":
            re.compile(r"Starting Time\s*=\s*(.+)"),

        "file(6)":
            re.compile(r"file\(6\)\s*=\s*phits_(MRCP|MFCP)_(AM|AF)_source_(.+?)_.+?_energy_.+?\.out",re.IGNORECASE),

        "maxcas":
            re.compile(r"maxcas\s*=\s*(\d+)"),

        "maxbch":
            re.compile(r"maxbch\s*=\s*(\d+)"),

        "source_energy":
            re.compile(r"e0\s*=\s*([0-9Ee.+-]+)"),

        "source_type":
        re.compile(r"proj\s*=\s*(\S+)", re.IGNORECASE),

        "parallelization_pattern":
            re.compile(r"\$(OMP|MPI)\s*=\s*(\d+)",re.IGNORECASE),

        "cpu_time_pattern":
            re.compile(r"cpu\s+time\s*=\s*(.+)$",re.IGNORECASE),

        "elapsed_time":
            re.compile(r"total cpu time\s*=\s*([0-9.]+)")
    }

    def parse_cpu_time(cpu_string):
        """
        Convert PHITS cpu time strings into seconds.

        Examples
        --------
        5.20 s.
        1 m. 5.20 s.
        2 h. 1 m. 5.20 s.
        1 d. 2 h. 1 m. 5.20 s.
        """

        total = 0.0

        day = re.search(r"([0-9.]+)\s*d\.", cpu_string)
        hour = re.search(r"([0-9.]+)\s*h\.", cpu_string)
        minute = re.search(r"([0-9.]+)\s*m\.", cpu_string)
        second = re.search(r"([0-9.]+)\s*s\.?", cpu_string)

        if day:
            total += float(day.group(1)) * 86400

        if hour:
            total += float(hour.group(1)) * 3600

        if minute:
            total += float(minute.group(1)) * 60

        if second:
            total += float(second.group(1))

        return total

    # Function for extracting metadata from a PHITS output and input file
    def extract_metadata(output_file: Path):

        # Determine the corresponding PHITS input file   
        input_files = output_file.with_suffix(".inp")

        # For extracting the particle transport times per batch
        batch_cpu_times = []

        start_collecting = False

        # Read the files
        results = {}

        lines = output_file.read_text(encoding="utf-8",errors="ignore").splitlines()

        # Search every line
        for i, line in enumerate(lines):

            for key, pattern in patterns.items():

                if key == "cpu_time_pattern":
                    continue
                if key in results:
                    continue
                if key in ["cpu_time_pattern", "parallelization_pattern"]:
                    continue

                match = pattern.search(line)

                if match:

                    value = match.group(1).strip()

                    if key in ["source_energy", "elapsed_time"]:
                        value = float(value)
                    elif key in ["maxcas", "maxbch"]:
                        value = int(value)

                    results[key] = value

            # Phantom, source organ and source type
            if "file(6)" in line.lower():

                match = patterns["file(6)"].search(line)

                if match:
                    
                    phantom_prefix = match.group(1).upper()
                    sex = match.group(2).upper()

                    phantom_key = f"{phantom_prefix}_{sex}"

                    results["phantom"] = phantom_names[phantom_key]
                    results["source_organ"] = match.group(3)

            # Particle transport times per batch
            if "bat[" in line:
                start_collecting = True

            if not start_collecting:
                continue

            match = patterns["cpu_time_pattern"].search(line)

            if match:

                batch_cpu_times.append(
                    parse_cpu_time(match.group(1))
                )
            # Job termination date/time
            if "job termination date" in line.lower():

                date_match = re.search(
                    r"job termination date\s*:\s*(.+)",
                    line,
                    re.IGNORECASE
                )

                if date_match:

                    results["termination_date"] = date_match.group(1).strip()

                # look ahead for the next line containing "time"
                for j in range(i + 1, min(i + 6, len(lines))):

                    if "time" in lines[j].lower():

                        time_match = re.search(
                            r"time\s*:\s*(.+)",
                            lines[j],
                            re.IGNORECASE
                        )

                        if time_match:

                            results["termination_time"] = time_match.group(1).strip()

                        break
                
        with open(input_files, encoding="utf-8", errors="ignore") as f:

            for line in f:

                match = patterns["parallelization_pattern"].search(line)

                if match:
                    results["parallelization"] = (
                        f"{match.group(1).upper()}-{match.group(2)}"
                    )
                    break

        # Format starting datetime (ISO 8601)
        start_time = results["starting_time"]

        start_time = (
            start_time.replace("h", ":")
                    .replace("m", ":")
                    .replace(" ", "")
        )

        results["starting_datetime"] = (
            f'{results["starting_date"]}T{start_time}'
        )

        results["termination_datetime"] = (
            f'{results["termination_date"].replace("/", "-")}T'
            f'{results["termination_time"]}'
        )

        # Remove old entries
        del results["starting_date"]
        del results["starting_time"]
        del results["termination_date"]
        del results["termination_time"]

        # Sum the batch particle transport times and store the total particle transport time
        results["particle_transport_time"] = sum(batch_cpu_times)

        # Calculate the initialization time. This represents the time it takes for the mesh generation, source definition, and 
        # other setup processes before the actual particle transport begins.
        results["initialization_time"] = results["elapsed_time"] - results["particle_transport_time"]

        # Order the results for better readability
        ordered_results = {
            "Starting Datetime": results["starting_datetime"],
            "Termination Datetime": results["termination_datetime"],

            "Phantom": results.get("phantom"),
            "Source Organ": results.get("source_organ"),
            "Source Type": results.get("source_type"),
            "Source Energy (MeV)": results.get("source_energy"),

            "Parallel Mode": results.get("parallelization"),
            "maxcas (nps if Geant4)": results.get("maxcas"),
            "maxbch": results.get("maxbch"),

            "Setup Time (s)": results.get("initialization_time"),
            "Execution Time (s)": results.get("particle_transport_time"),
            "Individual Wall Time (s)": results.get("elapsed_time"),

        }

        print("\nProcessing:")
        print(output_file.name)
        print(input_files.name)

        return ordered_results

    all_results = []

    if not output_files:
        print("\n[WARNING] No PHITS output files were found.")
        print(f"Search directory: {root}")
        return

    for output_file in output_files:
        all_results.append(
            extract_metadata(output_file)
        )

    with open(metadata_output_file,
              "w",
              newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=all_results[0].keys()
        )

        writer.writeheader()
        writer.writerows(all_results)

    print("\nMetadata extraction complete.")
    print(f"Processed {len(all_results)} simulation(s).")
    print(f"Output written to:")
    print(f"{metadata_output_file}")