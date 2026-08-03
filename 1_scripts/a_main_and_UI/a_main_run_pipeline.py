# This script acts as the main python file for running all the Python files below

# Import sub-Python files
from logging import config
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from a_main_and_UI.b_input_user_parameters import get_user_parameters
from c_database.a_parse_cell_and_csv import parse_cell_csv_inputs
from d_phits.a_phits_generating_inputs import phits_generate_inputs
from d_phits.b_phits_running_inputs import run_phits
from d_phits.c_phits_extracting_metadata import extract_metadata_stats
from d_phits.d_calculating_extra_metadata import calculate_metadata
from d_phits.e_phits_extracting_dose_and_SAFs import calculate_dose_and_safs
from d_phits.f_check_uncertainty import check_uncertainty
from e_geant4.a_geant4_setup_and_build_executable import build_geant4
from e_geant4.b_geant4_generating_inputs import geant4_generate_inputs
from e_geant4.c_geant4_running_inputs import run_geant4

def main():
    params = get_user_parameters()

    if params["simulation_code"] == "PHITS":

        parse_cell_csv_inputs()
        phits_generate_inputs()
        run_phits()
        extract_metadata_stats()
        calculate_metadata()
        calculate_dose_and_safs()
        check_uncertainty()

    elif params["simulation_code"] == "GEANT4":

        parse_cell_csv_inputs()
        build_geant4()
        geant4_generate_inputs()
        run_geant4()

if __name__ == "__main__":
    main()