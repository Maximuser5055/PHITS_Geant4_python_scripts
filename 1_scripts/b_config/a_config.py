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

############################
# File Paths
############################

# 1_scripts directory
SCRIPT_DIR = Path(__file__).resolve().parent.parent

# Parent Directory
ROOT = SCRIPT_DIR.parent

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

############################
# Geant4 Directories and Files
############################

GEANT4_DIR = ROOT / "3_geant4"

INTERNAL_DIR = ROOT / "3_geant4" / "Internal"

GEANT4_BUILD_DIR = INTERNAL_DIR / "build"

GEANT4_GENERATED_INPUTS_DIR = GEANT4_BUILD_DIR / "generated_inputs"

############################
# General Directories and Files
############################

ORGAN_ID_CSV = ROOT / "5_other_input_files/organ_ID_names.csv"

SOURCE_CSV = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/source_organs.csv")

DATABASE_FILE = ROOT / "1_scripts/c_database/b_organ_database.py"

RESULTS_DIR = ROOT / "4_results"

RERUN_CSV = ROOT / "4_results/5_rerun_required.csv"

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
]

PHANTOM_INPUT_GENERATION = "Both"

NPS = 10000

############################
# General Settings
############################

THREADS = 8

SOURCE_ENERGIES = [
    1.0,
]

UNCERTAINTY_LIMIT = 5.0