# This stores the configuration files of the internal dosimetry pipeline

# Import necessary libraries
from pathlib import Path
import platform

############################
# Operating System
############################

SYSTEM = None
IS_WINDOWS = False
IS_LINUX = False
IS_WSL = False

def detect_operating_system():

    global SYSTEM, IS_WINDOWS, IS_LINUX, IS_WSL

    SYSTEM = platform.system()  # "Windows", "Linux"
    IS_WINDOWS = SYSTEM == "Windows"
    IS_LINUX = SYSTEM == "Linux"

    # Detect if running inside WSL
    IS_WSL = (
        IS_LINUX and
        ("microsoft" in platform.release().lower() or
         "microsoft" in platform.version().lower())
    )

    return SYSTEM, IS_WINDOWS, IS_LINUX, IS_WSL

def update_config(setting, value):
    config_file = Path(__file__)
    
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
    
############################
# File Paths
############################

# 1_scripts directory
SCRIPT_DIR = Path(__file__).resolve().parent.parent

# Parent Directory
ROOT = SCRIPT_DIR.parent

############################
# General Directories and Files
############################

DATABASE_FILE = ROOT / "1_scripts/c_database/b_organ_database.py"

RESULTS_DIR = ROOT / "4_results"

RERUN_CSV = ROOT / "4_results/5_rerun_required.csv"

OTHER_INPUT_FILES_DIR = ROOT / "5_other_input_files"

ORGAN_ID_CSV = OTHER_INPUT_FILES_DIR / "organ_ID_names.csv"

SOURCE_CSV = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/source_organs.csv")

TARGET_REGION_CSV = OTHER_INPUT_FILES_DIR / "target_regions_Filipino.csv"

SKELETAL_MASSES_CSV = OTHER_INPUT_FILES_DIR / "skeletal_tissue_masses_ICRP116.csv"

############################
# PHITS Directories and Files
############################

PHITS_ROOT = Path(r"C:\phits")

CELL_FILES = {
    "AM": ROOT / "2_phits/phantoms/MRCP-AM.cell",
    "AF": ROOT / "2_phits/phantoms/MRCP-AF.cell",
}

INPUT_TEMPLATE_FILE = ROOT / "2_phits/template_input_files/1_template_MRCP_internal_input.inp"

INCLUDE_FILES_DIR = ROOT / "2_phits/phantoms"

NUCLEAR_DATA_FILE = PHITS_ROOT / "data/xsdir.jnd"

EGS5_DATA_DIR = PHITS_ROOT / "XS/egs"

GENERATED_INPUTS_DIR = ROOT / "2_phits/generated_inputs"

PHITS_BAT = "phits.bat"

PHITS_METADATA_FILE = RESULTS_DIR / "a_phits_all_simulations_log.csv"

############################
# Geant4 Directories and Files
############################

GEANT4_DIR = ROOT / "3_geant4"

INTERNAL_DIR = GEANT4_DIR / "Internal"

GEANT4_BUILD_DIR = INTERNAL_DIR / "build"

GEANT4_EXECUTABLE_FILE = GEANT4_BUILD_DIR / "Internal"

GEANT4_GENERATED_INPUTS_DIR = GEANT4_BUILD_DIR / "generated_inputs"

GEANT4_METADATA_FILE = RESULTS_DIR / "b_geant4_all_simulations_log.csv"

############################
# PHITS or Geant4
############################

SIMULATION_CODE = "GEANT4"

############################
# PHITS Settings
############################

PARALLELIZATION = "OMP"
    
MAXCAS = 500

MAXBCH = 20

PHANTOMS = [
    "AM",
    "AF",
]

SEXINFO = {
    "AM": "MALE",
    "AF": "FEMALE",
}

PHITS_SOURCE_TYPES = [
    "photon",
]

############################
# Geant 4 Settings
############################

REBUILD = False

GEANT4_SOURCE_TYPES = [
    "gamma",
    #"e-", # electron
]

PHANTOM_INPUT_GENERATION = "Both"

PHANTOM_NAMES = {
    "AM": "Adult Male",
    "AF": "Adult Female"
}

NPS = 100000000

SOURCE_TYPE_MAP = {
    "gamma": "photon",
    "e-": "electron",
}

############################
# General Settings
############################

THREADS = 8

SOURCE_ENERGIES = [
    1.0,
]

UNCERTAINTY_LIMIT = 5.0