# This script asks the user for various input files and PHITS parameters

# Import necessary libraries
import b_config.a_config as config
from pathlib import Path

def update_config(setting, value):
    config_file = Path(__file__).parent / "b_config" / "a_config.py"

    lines = config_file.read_text().splitlines()

    if isinstance(value, Path):
        new_line = f'{setting} = Path(r"{value}")'
    elif isinstance(value, str):
        new_line = f'{setting} = "{value}"'
    else:
        new_line = f"{setting} = {repr(value)}"

    for i, line in enumerate(lines):
        if line.startswith(f"{setting} ="):
            lines[i] = new_line
            break

    config_file.write_text("\n".join(lines))

def get_user_parameters():
    # User Interface
    print("=" * 60)
    print("Internal Dosimetry Pipeline Configuration")
    print("=" * 60)

    phits_root = input(f"PHITS root [Current = {config.PHITS_ROOT}]: ").strip()
    phits_root = Path(phits_root) if phits_root else config.PHITS_ROOT

    parallelization = (
        input(f"Parallelization (OMP/MPI) [Current = {config.PARALLELIZATION}]: ")
        .strip()
        .upper()
    )
    if parallelization == "":
        parallelization = config.PARALLELIZATION

    threads = input(f"Parallelization Threads [Current = {config.THREADS}]: ").strip()
    threads = int(threads) if threads else config.THREADS

    maxcas = input(f"PHITS maxcas (no. of particle histories per batch) [Current = {config.MAXCAS}]: ").strip()
    maxcas = int(maxcas) if maxcas else config.MAXCAS

    maxbch = input(f"PHITS maxbch (no. of batches) [Current = {config.MAXBCH}]: ").strip()
    maxbch = int(maxbch) if maxbch else config.MAXBCH

    uncertainty_limit = input(
        f"Maximum allowed statistical uncertainty (%) "
        f"[Current = {config.UNCERTAINTY_LIMIT}]: "
    ).strip()

    uncertainty_limit = (
        float(uncertainty_limit)
        if uncertainty_limit
        else config.UNCERTAINTY_LIMIT
    )

    while True:
        source_target_dir = input(
            f"Directory containing source_target_organs.csv "
            f"[Current = {config.SOURCE_TARGET_CSV.parent}]: "
        ).strip()

        source_target_csv = (
            Path(source_target_dir) / "source_target_organs.csv"
            if source_target_dir
            else config.SOURCE_TARGET_CSV
        )

        if source_target_csv.is_file():
            break

        print(f"\nError: '{source_target_csv}' was not found. Please try again.\n")

    while True:
        phantom_input_generation = input(
            f"Phantom Input Generation (AM/AF/BOTH) [Current = {config.PHANTOM_INPUT_GENERATION}]: "
        ).strip().upper()

        if phantom_input_generation == "":
            phantom_input_generation = config.PHANTOM_INPUT_GENERATION

        if phantom_input_generation in ("AM", "AF", "Both"):
            break

        print("\nError: Please enter AM, AF, or Both.\n")
   
    update_config("PHITS_ROOT", phits_root)
    update_config("PARALLELIZATION", parallelization)
    update_config("THREADS", threads)
    update_config("MAXCAS", maxcas)
    update_config("MAXBCH", maxbch)
    update_config("UNCERTAINTY_LIMIT", uncertainty_limit)
    update_config("SOURCE_TARGET_CSV", source_target_csv)
    update_config("PHANTOM_INPUT_GENERATION", phantom_input_generation)

    ############################
    # Create required directories
    ############################

    DIRECTORIES = [
        config.GENERATED_INPUTS_DIR,
        config.RESULTS_DIR,
    ]

    for directory in DIRECTORIES:
        directory.mkdir(parents=True, exist_ok=True)

    return {
            "phits_root": phits_root,
            "parallelization": parallelization,
            "threads": threads,
            "maxcas": maxcas,
            "maxbch": maxbch,
            "uncertainty_limit": uncertainty_limit,
            "source_target_csv": source_target_csv,
            "phantom": phantom_input_generation,
    }
    