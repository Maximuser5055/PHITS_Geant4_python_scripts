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

def detect_operating_system():

    global SYSTEM, IS_WINDOWS, IS_LINUX

    SYSTEM = platform.system()

    IS_WINDOWS = SYSTEM == "Windows"
    IS_LINUX = SYSTEM == "Linux"

    return SYSTEM, IS_WINDOWS, IS_LINUX

detect_operating_system()

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

SKELETAL_RESPONSE_FUNCTIONS_CSV = OTHER_INPUT_FILES_DIR / "ICRP116-Table-D-1-Updated.csv"

############################
# PHITS Directories and Files
############################

PHITS_ROOT_WINDOWS = Path(r"C:\phits")

PHITS_ROOT_LINUX = Path(r"/home/clarence/Software/phits/phits")

if IS_WINDOWS:
    PHITS_INSTALLATION_DIR = PHITS_ROOT_WINDOWS
    PHITS_EXECUTABLE = PHITS_INSTALLATION_DIR / "bin/phits.bat"

elif IS_LINUX:
    PHITS_INSTALLATION_DIR = PHITS_ROOT_LINUX
    PHITS_EXECUTABLE = PHITS_INSTALLATION_DIR / "bin/phits.sh"

else:
    raise RuntimeError(
        f"Unsupported operating system: {SYSTEM}"
    )

CELL_FILES = {
    "AM": ROOT / "2_phits/phantoms/MRCP-AM.cell",
    "AF": ROOT / "2_phits/phantoms/MRCP-AF.cell",
}

INPUT_TEMPLATE_FILE = ROOT / "2_phits/template_input_files/1_template_MRCP_internal_input.inp"

INCLUDE_FILES_DIR = ROOT / "2_phits/phantoms"

GENERATED_INPUTS_DIR = ROOT / "2_phits/generated_inputs"

PHITS_METADATA_FILE = RESULTS_DIR / "a_phits_all_simulations_log.csv"

############################
# Geant4 Directories and Files
############################

GEANT4_DIR = ROOT / "3_geant4"

INTERNAL_DIR = GEANT4_DIR / "Internal"

GEANT4_BUILD_DIR = INTERNAL_DIR / "build"

INCLUDE_DIR = INTERNAL_DIR / "include"

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
    
MAXCAS = 100

MAXBCH = 10

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
    "electron"
]

SPONGIOSA_IDS = [
    1400,
    2500,
    2700,
    2900,
    4000,
    4200,
    4400,
    4600,
    4800,
    5000,
    5200,
    5400,
    5600,
]

############################
# Geant 4 Settings
############################

REBUILD = False

GEANT4_SOURCE_TYPES = [
    "gamma",
    "e-",
]

PHANTOM_INPUT_GENERATION = "Both"

PHANTOM_NAMES = {
    "AM": "Adult Male",
    "AF": "Adult Female"
}

NPS = 1000

GEANT4_SOURCE_TYPE_MAP = {
    "gamma": "photon",
    "e-": "electron",
}

############################
# General Settings
############################

MEV_TO_J = 1.6021766339999e-13

THREADS = 12

SOURCE_ENERGIES = [10.0]

SELECTED_SOURCE_TYPE = "electron"

# Mass fraction from ICRP 116 Table 3.2
MARROW_MASS_KG = {
    1400: .0269,
    2500: .0093,
    2700: .0889,
    2900: .0784,
    4000: .0094,
    4200: .2052,
    4400: .1888,
    4600: .0328,
    4800: .0456,
    5000: .1888,
    5200: .1439,
    5400: .1159,
    5600: .0363,
}

FLUENCE_SOURCE_TYPES = [
    "photon",
    "gamma",
]

# Logarithmic bins for fluence-to-dose response functions for RBM and endosteum dosimetry
ENERGY_BINS = 100

ENERGY_MIN = 0.01

ENERGY_MAX = 10

UNCERTAINTY_LIMIT = 5.0