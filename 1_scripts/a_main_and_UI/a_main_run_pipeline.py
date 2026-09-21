# This script acts as the main python file for running all the Python files below

# Import sub-Python files
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

def main():

    from a_main_and_UI.b_input_user_parameters import get_user_parameters
    params = get_user_parameters() # Get user inputs

    if params["redo_saf_calculations"]:
        from c_database.a_parse_cell_and_csv import parse_cell_csv_inputs
        parse_cell_csv_inputs(params) # Parse immediately

        from d_phits import (
            phits_generate_inputs,
            run_phits,
            phits_extract_metadata_stats,
            phits_calculate_dose_and_safs,
            phits_calculate_marrow_endosteum_SAFs,
        )

        from e_geant4 import (
            build_geant4,
            geant4_generate_inputs,
            run_geant4,
            geant4_extract_metadata_stats,
            geant4_calculate_dose_and_SAFs,
            geant4_calculate_marrow_endosteum_SAFs,
        )

        from f_simulation_and_SAFs_further_analysis import (
            calculate_extra_metadata,
            combine_target_organs_and_calculate_true_dose_and_SAFs,
            update_master_saf_database,
            create_publishable_saf_database,
            # Not implemented yet
            #limiting_SAF_approach_calculation
        )

        # --------------------------------------------------------
        # Transport simulation
        # --------------------------------------------------------

        if params["simulation_code"] == "PHITS":

            phits_generate_inputs(params)
            run_phits(params)
            phits_extract_metadata_stats()
            calculate_extra_metadata(params)
            phits_calculate_dose_and_safs(params)
            phits_calculate_marrow_endosteum_SAFs(params)

        elif params["simulation_code"] == "GEANT4":

            build_geant4(params)
            geant4_generate_inputs(params)
            run_geant4(params)
            geant4_extract_metadata_stats()
            calculate_extra_metadata(params)
            geant4_calculate_dose_and_SAFs(params)
            geant4_calculate_marrow_endosteum_SAFs(params)

        else:
            raise ValueError(f"Unsupported simulation code: {params['simulation_code']}")
        
        combine_target_organs_and_calculate_true_dose_and_SAFs(params)
        update_master_saf_database(params)
        create_publishable_saf_database(params)
        # Not implemented yet
        # limiting_SAF_approach_calculation(params)

    # else:
        # Not implemented yet
        # from g_radionuclides.extract_ICRP_107_radionuclide_information import extract_radionuclide_info
        # extract_radionuclide_info(params)
        # from h_S_values_and_whole_body_doses.calculate_S_values import calculate_S_values
        # calculate_S_values(params)

if __name__ == "__main__":
    main()