# This script acts as the main python file for running all the Python files below

# Import sub-Python files
from a_main_and_UI.b_input_user_parameters import get_user_parameters, update_config
from c_database.a_parse_cell_and_csv import parse_cell_csv_inputs
from d_phits.a_phits_generating_inputs import generate_inputs
from d_phits.b_phits_running_inputs import run_phits
from d_phits.c_phits_extracting_metadata import extract_metadata_stats
from d_phits.d_calculating_extra_metadata import calculate_metadata
from d_phits.e_phits_extracting_dose_and_SAFs import calculate_dose_and_safs
from d_phits.f_check_uncertainty import check_uncertainty

def main():
    get_user_parameters()
    parse_cell_csv_inputs()
    generate_inputs()
    run_phits()
    extract_metadata_stats()
    calculate_metadata()
    calculate_dose_and_safs()
    check_uncertainty()

if __name__ == "__main__":
    main()