# This script automates the running of Geant4 simulations using the generated input files.

# Import necessary libraries
#from pathlib import Path
import subprocess
from datetime import datetime
from time import perf_counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import pandas as pd

import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS
from b_config.b_phantom_registry import get_phantom, get_phantom_group
from e_geant4.a_geant4_setup_and_build_executable import find_geant4make

def run_geant4(params):

    geant4make_path = find_geant4make()

    # Configs
    executable = config.GEANT4_EXECUTABLE_FILE
    input_root = config.GEANT4_GENERATED_INPUTS_DIR
    current_working_dir = config.GEANT4_BUILD_DIR
    source_csv = config.SOURCE_CSV

    #########################################
    # User Interface
    #########################################

    print("\n")
    print("=" * 90)
    print("Geant4 SAF Simulation Launcher")
    print("=" * 90)

    phantom_selection = params["phantom"]

    selected_phantoms = get_phantom_group(phantom_selection)

    # ==========================================================
    # Source Organ Selection
    # ==========================================================

    if not source_csv.is_file():
        raise FileNotFoundError(f"Source organ file not found:\\n{source_csv}")

    # Current source_organs.csv schema:
    # Source region | Acronym | ID number(s)
    source_df = pd.read_csv(source_csv, encoding="utf-8-sig")
    source_df.columns = [str(column).strip() for column in source_df.columns]

    print("\n")
    print("=" * 90)
    print("Source Organ Selection")
    print("=" * 90)
    print(f"Source organ file:\n{source_csv}")

    required_columns = {"Source region", "Acronym", "ID number(s)"}
    missing_columns = required_columns.difference(source_df.columns)
    if missing_columns:
        raise ValueError(
            "The source-organ CSV is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}. "
            f"Found columns: {source_df.columns.tolist()}\n"
            f"File: {source_csv}"
        )

    for column in required_columns:
        source_df[column] = source_df[column].astype("string").str.strip()

    source_df = source_df.dropna(subset=["Source region", "Acronym", "ID number(s)"])
    source_df = source_df[
        source_df["Source region"].ne("")
        & source_df["Acronym"].ne("")
        & source_df["ID number(s)"].ne("")
    ].drop_duplicates(subset=["Acronym"], keep="first").reset_index(drop=True)

    if source_df.empty:
        raise ValueError(f"No valid source-organ rows were found in {source_csv}")

    def normalize_source_ids(value):
        """Normalize grouped source IDs such as 9700_9900 or 9700,9900."""
        text = str(value).strip()
        # Tolerate legacy source_id.txt files containing literal newline escapes.
        text = text.replace("\\n", "").replace("\n", "").strip()
        if not text or text.lower() == "nan":
            return None
        parts = [part.strip() for part in text.replace("_", ",").split(",")]
        if not parts or any(not part.isdigit() for part in parts):
            raise ValueError(f"Invalid source ID value: {value!r}")
        return ",".join(dict.fromkeys(parts))

    def source_id_key(value):
        normalized = normalize_source_ids(value)
        return tuple(sorted(normalized.split(","))) if normalized else tuple()

    available_sources = []
    seen_id_keys = set()
    for _, row in source_df.iterrows():
        normalized_ids = normalize_source_ids(row["ID number(s)"])
        if not normalized_ids:
            continue
        key = source_id_key(normalized_ids)
        if key in seen_id_keys:
            continue
        seen_id_keys.add(key)
        available_sources.append({
            "source_region": str(row["Source region"]).strip(),
            "acronym": str(row["Acronym"]).strip(),
            "ids_text": str(row["ID number(s)"]).strip(),
            "normalized_ids": normalized_ids,
        })

    if not available_sources:
        raise ValueError(f"No valid source-organ IDs were found in {source_csv}")

    for i, source in enumerate(available_sources, start=1):
        print(
            f"[{i}] {source['source_region']} "
            f"(Acronym: {source['acronym']}; IDs: {source['ids_text']})"
        )
    print("[A] All source organs")

    while True:
        source_choice = input("\nSelect source organ group(s): ").strip().upper()
        if source_choice == "A":
            selected_sources = available_sources.copy()
            break
        try:
            indices = [int(item.strip()) - 1 for item in source_choice.split(",")]
            if not indices or any(
                index < 0 or index >= len(available_sources) for index in indices
            ):
                raise ValueError
            selected_sources = [
                available_sources[index] for index in dict.fromkeys(indices)
            ]
            break
        except ValueError:
            print(
                "\nInvalid selection. Enter valid row numbers "
                "separated by commas, or A."
            )

    print("\nSelected source organ group(s):")
    for source in selected_sources:
        print(
            f"  - {source['source_region']} "
            f"({source['acronym']}; IDs: {source['ids_text']})"
        )

    selected_source_keys = {
        source_id_key(source["normalized_ids"]) for source in selected_sources
    }

    # ==========================================================
    # Find Geant4 input files
    # ==========================================================

    input_files = []
    for phantom in selected_phantoms:
        phantom_dir = input_root / phantom
        if not phantom_dir.exists():
            print(f"\nSkipping {phantom}: input directory does not exist.")
            continue

        for infile in sorted(phantom_dir.rglob("*.in")):
            source_id_file = infile.parent / "source_id.txt"
            if not source_id_file.is_file():
                print(f"[WARNING] Missing source_id.txt for {infile.name}")
                continue

            source_id = source_id_file.read_text(encoding="utf-8-sig").strip()
            try:
                source_key = source_id_key(source_id)
            except ValueError:
                print(
                    f"[WARNING] Invalid source_id.txt contents for {infile.name}: "
                    f"{source_id!r}"
                )
                continue

            if source_key in selected_source_keys:
                input_files.append(infile)

    input_files = sorted(set(input_files))

    if not input_files:
        print(
            "\nNo Geant4 input files were found "
            "for the selected phantom/source organs."
        )
        return

    # ======================================================
    # Existing simulation check
    # ======================================================

    print("=" * 90)
    print("Existing Simulation Check")
    print("=" * 90)

    print(
        "[1] Skip jobs that have already been simulated"
    )

    print(
        "[2] Re-run completed jobs (overwrite outputs)"
    )

    while True:

        overwrite_choice = input(
            "Enter your choice (1-2): "
        ).strip()

        if overwrite_choice in {"1", "2"}:
            break

        print(
            "Invalid choice. "
            "Please enter 1 or 2."
        )

    skip_completed = (
        overwrite_choice == "1"
    )

    # ==========================================================
    # Check completion
    # ==========================================================

    def job_completed(infile):

        outfile = infile.with_name(
            infile.stem.replace(
                "geant4_",
                "geant4_deposit_"
            ) + ".out"
        )

        timing_file = infile.with_name(
            infile.stem.replace(
                "geant4_",
                "geant4_timing_"
            ) + ".txt"
        )
        
        required_outputs = [outfile, timing_file,]

        # A simulation is considered completed only when every
        # required output exists AND is non-zero in size.
        return all(output.is_file() and output.stat().st_size > 0
                  for output in required_outputs)

    # ==========================================================
    # Filter completed jobs
    # ==========================================================

    if skip_completed:

        original = len(input_files)

        input_files = [
            infile
            for infile in input_files
            if not job_completed(infile)
        ]

        skipped = (
            original -
            len(input_files)
        )

        print(
            f"\nSkipping {skipped} "
            f"completed job(s)."
        )

        print(
            f"Remaining jobs: "
            f"{len(input_files)}"
        )

    # ==========================================================
    # Check whether anything remains
    # ==========================================================

    if not input_files:

        print(
            "\nNo simulations to run. "
            "All jobs have already been completed."
        )

        return

    # ==========================================================
    # Parallelization
    # ==========================================================

    total_threads = os.cpu_count()

    threads_per_job = params["threads"]

    max_parallel = max(
        1,
        total_threads // threads_per_job
    )

    print(
        f"\nFound {len(input_files)} "
        f"input files."
    )

    print(
        f"Total CPU threads  : "
        f"{total_threads}"
    )

    print(
        f"Threads per Geant4 : "
        f"{threads_per_job}"
    )

    print(
        f"Parallel Geant4 jobs: "
        f"{max_parallel}"
    )

    print(
        f"\nRunning up to "
        f"{max_parallel} jobs simultaneously.\n"
    )

    # ==========================================================
    # Run one Geant4 job
    # ==========================================================

    def run_job(infile):

        relative = infile.relative_to(
            input_root
        )

        phantom = relative.parts[0]

        outfile = infile.with_name(
            infile.stem.replace(
                "geant4_",
                "geant4_deposit_"
            ) + ".out"
        )

        source_id_file = (
            infile.parent /
            "source_id.txt"
        )

        if not source_id_file.exists():

            print(
                f"Missing source_id.txt for "
                f"{infile.name}"
            )

            return infile, 1

        source_id = (
            source_id_file
            .read_text()
            .strip()
        )

        cmd = (
            f'source "{geant4make_path}" && '
            f'"{executable}" '
            f'-i {source_id} '
            f'-m "{infile.resolve()}" '
            f'-o "{outfile.resolve()}"'
        )

        # ------------------------------------------------------
        # AF phantom flag
        # ------------------------------------------------------

        phantom_spec = get_phantom(phantom)

        if phantom_spec.sex == "AF":
            cmd += " -f"

        # Measure total elapsed wall time in Python and write the summary
        # beside this simulation's input/output files.
        timing_file = infile.with_name(
            infile.stem.replace(
                "geant4_",
                "geant4_timing_"
            ) + ".txt"
        )

        start_clock = perf_counter()
        returncode = -1
        launch_error = None

        try:
            process = subprocess.Popen(
                ["bash", "-c", cmd],
                cwd=current_working_dir
            )
            returncode = process.wait()
        except Exception as exc:
            launch_error = repr(exc)
            print(f"Error launching {infile.name}: {exc}")
        finally:
            elapsed_seconds = perf_counter() - start_clock
            try:
                with timing_file.open("w", encoding="utf-8") as f:
                    f.write("Geant4 Timing Summary\n")
                    f.write("=====================\n\n")
                    f.write(f"Total Wall Time (s): {elapsed_seconds:.6f}\n")
                    if launch_error:
                        f.write(f"Launch Error: {launch_error}\n")
            except OSError as exc:
                print(f"WARNING: Could not write timing file {timing_file}: {exc}")

        if returncode == 0:
            print(f"Finished {infile.name}")
        else:
            print(f"Failed {infile.name} (return code {returncode})")

        return infile, returncode

    # ==========================================================
    # Run jobs in parallel
    # ==========================================================

    with ThreadPoolExecutor(
        max_workers=max_parallel
    ) as executor:

        futures = [
            executor.submit(
                run_job,
                infile
            )
            for infile in input_files
        ]

        for future in as_completed(
            futures
        ):

            infile, returncode = (
                future.result()
            )

            if returncode == 0:

                print(
                    f"✓ {infile.name} "
                    f"completed successfully"
                )

            else:

                print(
                    f"✗ {infile.name} "
                    f"failed "
                    f"(return code "
                    f"{returncode})"
                )

    print(
        "\nAll Geant4 jobs finished."
    )