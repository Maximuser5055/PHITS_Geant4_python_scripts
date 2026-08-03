# This stores the configuration files of the internal dosimetry pipeline

# Import necessary libraries
from pathlib import Path

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

GEANT4_FOLDER = ROOT / "3_geant4"

INTERNAL_FOLDER = ROOT / "3_geant4" / "Internal"

############################
# General Directories and Files
############################

ORGAN_ID_CSV = ROOT / "5_other_input_files/organ_ID_names.csv"

SOURCE_CSV = Path(r"/home/clarence/Geant4_SAF_Calculations/python_scripts/5_other_input_files/source_organs.csv")

DATABASE_FILE = ROOT / "1_scripts/d_organ_database.py"

RESULTS_DIR = ROOT / "4_results"

RERUN_CSV = ROOT / "4_results/5_rerun_required.csv"

############################
# PHITS or Geant4
############################

SIMULATION_CODE = "PHITS"

############################
# PHITS Settings
############################

PARALLELIZATION = "OMP"

THREADS = 8
    
MAXCAS = 500

MAXBCH = 20

PHANTOMS = [
    "AM",
    "AF",
]

PHANTOM_INPUT_GENERATION = "Both"

SEXINFO = {
    "AM": "MALE",
    "AF": "FEMALE",
}

SOURCE_TYPES = [
    "photon",
]

SOURCE_ENERGIES = [
    1.0,
]

############################
# Geant 4 Settings
############################

############################
# General Settings
############################

UNCERTAINTY_LIMIT = 5.0