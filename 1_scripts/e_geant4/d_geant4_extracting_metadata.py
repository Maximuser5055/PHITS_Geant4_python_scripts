# This script extracts useful metadata from Geant4 output files using keyword matching.

# Import necessary libraries
import re
from pathlib import Path
import csv

import b_config.a_config as config


def geant4_extract_metadata_stats():

    root = config.GEANT4_GENERATED_INPUTS_DIR

    timing_files = sorted(root.rglob("geant4_timing_*.txt"))

    metadata_output_file = (
        config.RESULTS_DIR /
        "b_geant4_all_simulations_log.csv"
    )

    PHANTOM_NAMES = {
        "AM": "Adult Male",
        "AF": "Adult Female",
    }

    timing_patterns = {

        "starting_datetime":
            re.compile(r"Starting\s+Datetime\s*:\s*(.+)"),

        "termination_datetime":
            re.compile(r"Termination\s+Datetime\s*:\s*(.+)"),

        "setup_time":
            re.compile(r"Setup\s+Time.*?:\s*([0-9.]+)"),

        "execution_time":
            re.compile(r"Execution\s+Time.*?:\s*([0-9.]+)"),

        "total_wall_time":
            re.compile(r"Total\s+Wall\s+Time.*?:\s*([0-9.]+)",
                       re.IGNORECASE),
    }

    filename_pattern = re.compile(
        r"geant4_MRCP_(AM|AF)_source_(.+?)_(.+?)_energy_([0-9Ee.+-]+)\.in",
        re.IGNORECASE
    )

    thread_pattern = re.compile(
        r"^\s*/run/numberOfThreads\s+(\d+)",
        re.IGNORECASE
    )

    beamon_pattern = re.compile(
        r"/run/beamOn\s+(\d+)",
        re.IGNORECASE
    )

    def extract_metadata(timing_file: Path):

        results = {}

        # --------------------------
        # Timing file
        # --------------------------

        with open(timing_file, encoding="utf-8",
                  errors="ignore") as f:

            for line in f:

                for key, pattern in timing_patterns.items():

                    if key in results:
                        continue

                    match = pattern.search(line)

                    if match:

                        value = match.group(1).strip()

                        if key in [
                            "setup_time",
                            "execution_time",
                            "total_wall_time",
                        ]:
                            value = float(value)

                        results[key] = value

        # --------------------------
        # Corresponding .in file
        # --------------------------

        input_files = timing_file.with_name(
            timing_file.name
                .replace("geant4_timing_", "geant4_")
                .replace(".txt", ".in")
        )

        match = filename_pattern.match(input_files.name)

        if not match:
            raise RuntimeError(
                f"Cannot parse filename:\n{input_files.name}"
            )

        phantom = PHANTOM_NAMES[match.group(1)]
        source_organ = match.group(2)
        source_type = match.group(3)
        source_energy = float(match.group(4))

        results["phantom"] = phantom
        results["source_organ"] = source_organ
        results["source_type"] = source_type
        results["source_energy"] = source_energy

        parallel_mode = "SEQ"

        nps = None

        with open(input_files,
                  encoding="utf-8",
                  errors="ignore") as f:

            for line in f:

                stripped = line.strip()

                # Ignore commented commands

                if stripped.startswith("#"):
                    continue

                thread_match = thread_pattern.search(line)

                if thread_match:

                    nthreads = int(thread_match.group(1))

                    if nthreads > 1:
                        parallel_mode = f"MT-{nthreads}"
                    else:
                        parallel_mode = "SEQ"

                beamon_match = beamon_pattern.search(line)

                if beamon_match:
                    nps = int(beamon_match.group(1))

        results["parallel_mode"] = parallel_mode
        results["nps"] = nps
        results["maxbch"] = 1

        ordered_results = {

            "Starting Datetime":results.get("starting_datetime"),
            "Termination Datetime":results.get("termination_datetime"),

            "Phantom":results.get("phantom"),
            "Source Organ":results.get("source_organ"),
            "Source Type":results.get("source_type"),
            "Source Energy (MeV)":results.get("source_energy"),

            "Parallel Mode":results.get("parallel_mode"),
            "maxcas (nps if Geant4)":results.get("nps"),
            "maxbch (1 if Geant4)":results.get("maxbch"),

            "Setup Time (s)":results.get("setup_time"),
            "Execution Time (s)":results.get("execution_time"),
            "Individual Wall Time (s)":results.get("total_wall_time"),
        }

        print(f"\nProcessing:")
        print(f"{timing_file.name}")
        print(f"{input_files.name}")

        return ordered_results

    all_results = []

    if not timing_files:
        print("\n[WARNING] No Geant4 timing files were found.")
        print(f"Search directory: {root}")
        return

    for timing_file in timing_files:
        all_results.append(
            extract_metadata(timing_file)
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
