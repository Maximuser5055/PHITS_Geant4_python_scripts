# This script automates input and macro file generation for Geant4.
# Reads the example.in and source.mac as templates then generates new input and macro files for each source organ and energy combination.

#from pathlib import Path
#import shutil

import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS

def geant4_generate_inputs(params):

    total_files = 0

    # Read templates
    input_template = (config.GEANT4_BUILD_DIR / "example.in").read_text()
    source_template = (config.GEANT4_BUILD_DIR / "source.mac").read_text()

    base_output_dir = config.GEANT4_GENERATED_INPUTS_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Adult male or female phantom or both
    if params["phantom"] == "AM":
        phantoms = ["AM"]
    elif params["phantom"] == "AF":
        phantoms = ["AF"]
    else:
        phantoms = ["AM", "AF"]

    for phantom in phantoms:

        threads = params["threads"]
        nps = params["nps"]      
        source_type = params["source_type"]

        for energy in params["source_energies"]:

            for region, organ_name in SOURCE_ORGANS[phantom].items():

                safe_name = (
                    organ_name
                    .replace(",", "")
                    .replace(" ", "_")
                )

                basename = (
                    f"Geant4_MRCP_"
                    f"{phantom}_"
                    f"source_{safe_name}_"
                    f"{source_type}_"
                    f"energy_{energy}"
                )

                #############################
                # Job directory
                #############################

                job_dir = base_output_dir / phantom / basename
                job_dir.mkdir(parents=True, exist_ok=True)

                #############################
                # source.mac
                #############################

                mac_text = source_template

                mac_text = mac_text.replace(
                    "{{PARTICLE}}",
                    source_type
                )

                mac_text = mac_text.replace(
                    "{{ENERGY}}",
                    str(energy)
                )

                mac_filename = f"{basename}.mac"

                (job_dir / mac_filename).write_text(mac_text)

                #############################
                # example.in
                #############################

                in_text = input_template

                in_text = in_text.replace(
                    "{{THREADS}}",
                    str(threads)
                )

                in_text = in_text.replace(
                    "{{NPS}}",
                    str(nps)
                )

                in_text = in_text.replace(
                    "{{SOURCE_MAC}}",
                    mac_filename
                )

                in_filename = f"{basename}.in"

                (job_dir / in_filename).write_text(in_text)

                (job_dir / "source_id.txt").write_text(str(region))

                total_files += 1

    print(f"\nGenerated Geant4 inputs in {base_output_dir}")

    num_phantoms = len(phantoms)
    num_organs_per_phantom = len(SOURCE_ORGANS[phantoms[0]])
    num_source_types = int(params["source_type"] in config.GEANT4_SOURCE_TYPES)
    num_energy_bins = len(params["source_energies"])

    print(
        f"Generated {total_files} Geant4 input file(s) [{int(total_files/2)} .in and {int(total_files/2)} .mac files]"
        f"({num_phantoms} phantom type(s) × "
        f"{num_organs_per_phantom} source organ(s) × "
        f"{num_source_types} source type(s) × "
        f"{num_energy_bins} source energy bin(s))."
    )