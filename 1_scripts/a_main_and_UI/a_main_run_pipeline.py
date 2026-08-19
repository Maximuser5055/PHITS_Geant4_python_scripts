# This script acts as the main python file for running all the Python files below

# Import sub-Python files
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from a_main_and_UI.b_input_user_parameters import get_user_parameters

from c_database.a_parse_cell_and_csv import parse_cell_csv_inputs

from d_phits.a_phits_generating_inputs import phits_generate_inputs
from d_phits.b_phits_running_inputs import run_phits
from d_phits.c_phits_extracting_metadata import phits_extract_metadata_stats
from d_phits.d_phits_extracting_dose_and_SAFs import phits_calculate_dose_and_safs
from d_phits.e_phits_calculate_marrow_endosteum_SAFs import phits_calculate_marrow_endosteum_SAFs

from e_geant4.a_geant4_setup_and_build_executable import build_geant4
from e_geant4.b_geant4_generating_inputs import geant4_generate_inputs
from e_geant4.c_geant4_running_inputs import run_geant4
from e_geant4.d_geant4_extracting_metadata import geant4_extract_metadata_stats
from e_geant4.e_geant4_extracting_dose_and_SAFs import geant4_calculate_dose_and_SAFs
from e_geant4.f_geant4_calculate_marrow_endosteum_SAFs import geant4_calculate_marrow_endosteum_SAFs

from f_simulation_and_SAFs_further_analysis.a_calculating_extra_metadata import calculate_extra_metadata
from f_simulation_and_SAFs_further_analysis.b_configure_target_organs_and_calculate_true_dose_and_SAFs import combine_target_organs_and_calculate_true_dose_and_SAFs
# Not updated
from f_simulation_and_SAFs_further_analysis.c_check_uncertainty import check_uncertainty
# ===========
from f_simulation_and_SAFs_further_analysis.d_create_and_update_master_saf_database import update_master_saf_database
from f_simulation_and_SAFs_further_analysis.e_create_and_update_publishable_saf_database import create_publishable_saf_database
def main():

    params = get_user_parameters()
    
    if params["simulation_code"] == "PHITS":

        parse_cell_csv_inputs()
        phits_generate_inputs(params)
        run_phits()
        phits_extract_metadata_stats()
        calculate_extra_metadata(params)
        phits_calculate_dose_and_safs()
        phits_calculate_marrow_endosteum_SAFs(params)

    elif params["simulation_code"] == "GEANT4":

        parse_cell_csv_inputs()
        build_geant4(params)
        geant4_generate_inputs(params)
        run_geant4(params)
        geant4_extract_metadata_stats()
        calculate_extra_metadata(params)
        geant4_calculate_dose_and_SAFs()
        geant4_calculate_marrow_endosteum_SAFs(params)

    combine_target_organs_and_calculate_true_dose_and_SAFs(params)
    check_uncertainty(params)
    update_master_saf_database(params)
    create_publishable_saf_database(params)

if __name__ == "__main__":
    main()