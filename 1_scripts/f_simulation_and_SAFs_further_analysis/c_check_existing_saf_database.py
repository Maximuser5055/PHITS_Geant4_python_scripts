# This script checks whether the required publishable SAF databases exist.

# Import necessary libraries
import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom_group

# ============================================================
# EXISTING PUBLISHABLE SAF DATABASE CHECK
# ============================================================

def get_required_saf_database_files(phantom, simulation_code):
    """
    Return all required publishable SAF and STD files for
    the selected phantom.f

    Both photon and electron databases are always required.
    """

    simulation_code = simulation_code.upper()

    publishable_dirs = {
        "PHITS": config.RESULTS_PHITS_PUBLISHABLE_SAF_DATABASE_DIR,
        "GEANT4": config.RESULTS_GEANT4_PUBLISHABLE_SAF_DATABASE_DIR,
    }

    if simulation_code not in publishable_dirs:
        raise ValueError(
            f"Unsupported simulation code: {simulation_code}"
        )

    phantom_codes = get_phantom_group(phantom)

    publishable_dir = publishable_dirs[simulation_code]
    required_files = []

    for phantom_code in phantom_codes:
        phantom_name = phantom_code.lower()

        for source_name in ("photons", "electrons"):
            required_files.extend([
                publishable_dir / f"{phantom_name}_{source_name}_saf.csv",
                publishable_dir / f"{phantom_name}_{source_name}_std.csv",
            ])

    return required_files

def check_existing_saf_database(phantom, simulation_code):
    """
    Check whether the complete publishable SAF database required
    by the selected phantom and simulation code exists.

    Both photon and electron SAF databases are required.

    This function only checks whether the required SAF files exist.
    The user interface in b_input_user_parameters.py decides whether
    to reuse the existing database or redo the SAF calculations.
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
        }

    return {
        "exists": True,
        "complete": True,
        "required_files": required_files,
        "existing_files": existing_files,
        "missing_files": [],
    }