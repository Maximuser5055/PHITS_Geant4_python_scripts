# This script automates input and macro file generation for Geant4.
# Reads the example.in and source.mac as templates then generates new input and macro files for each source organ and energy combination.

# temporary
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pathlib import Path
import shutil

import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS

def geant4_generate_inputs():

    # Read templates
    input_template = (config.GEANT4_BUILD_DIR / "example.in").read_text()
    source_template = (config.GEANT4_BUILD_DIR / "source.mac").read_text()

    base_output_dir = config.GEANT4_GENERATED_INPUTS_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Adult male or female phantom or both
    if config.PHANTOM_INPUT_GENERATION == "AM":
        phantoms = ["AM"]
    elif config.PHANTOM_INPUT_GENERATION == "AF":
        phantoms = ["AF"]
    else:
        phantoms = ["AM", "AF"]

    for phantom in phantoms:

        threads = config.THREADS
        nps = config.NPS       
        particle = config.GEANT4_SOURCE_TYPES[0]

        for energy in config.SOURCE_ENERGIES:

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
                    f"{particle}_"
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
                    particle
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

    print(f"Generated Geant4 inputs in {base_output_dir}")

geant4_generate_inputs()