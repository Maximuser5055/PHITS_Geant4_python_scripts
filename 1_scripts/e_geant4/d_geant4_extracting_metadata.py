# This script extracts useful metadata from Geant4 output files using keyword matching.

# Import necessary libraries
import re
from pathlib import Path
import csv

import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom

def geant4_extract_metadata_stats():

    # Configs
    root = config.GEANT4_GENERATED_INPUTS_DIR
    timing_files = sorted(root.rglob("geant4_timing_*.txt"))
    metadata_output_file = (config.RESULTS_GEANT4_DIR / "a_geant4_all_simulations_log.csv")

    # Python runner writes only total elapsed wall time.
    # Accept both the current label and the older C++ timing-file label.
    total_wall_time_pattern = re.compile(
        r"Total\s+(?:Elapsed\s+)?Wall\s+Time\s*\(s\)\s*:\s*([0-9.eE+-]+)",
        re.IGNORECASE,
    )

    # Match the source particle type explicitly so organ names containing
    # underscores (e.g. Cortical_bone) are not split into Source Organ/Type.
    filename_pattern = re.compile(
        r"^geant4_(?P<phantom>.+?)_source_(?P<source_organ>.+)_"
        r"(?P<source_type>gamma|photon|electron|positron|proton|neutron)_"
        r"energy_(?P<energy>[0-9Ee.+-]+)\.in$",
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
        # Timing file: total elapsed wall time only
        # --------------------------

        with open(timing_file, encoding="utf-8", errors="ignore") as f:
            for line in f:
                match = total_wall_time_pattern.search(line)
                if match:
                    results["total_wall_time"] = float(match.group(1))
                    break

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

        phantom_code = match.group("phantom").upper()
        phantom_spec = get_phantom(phantom_code)

        source_organ = match.group("source_organ")
        source_type = match.group("source_type")
        source_energy = float(match.group("energy"))

        results["phantom"] = phantom_spec.display_name
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

            "Phantom": results.get("phantom"),
            "Source Organ": results.get("source_organ"),
            "Source Type": results.get("source_type"),
            "Source Energy (MeV)": results.get("source_energy"),

            "Parallel Mode": results.get("parallel_mode"),
            "maxcas (nps if Geant4)": results.get("nps"),
            "maxbch (1 if Geant4)": results.get("maxbch"),

            "Total Elapsed Wall Time (s)": results.get("total_wall_time"),
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
