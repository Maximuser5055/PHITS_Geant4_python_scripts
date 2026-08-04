# This script generates PHITS input files based on an input template and the b_organ_database.py. 
# It replaces placeholders in the template with actual values for phantom type, source energy, source region, and other parameters, and 
# saves the generated input files in their own directories under the base output directory. 
# It also copies the necessary phantom files to the mother directory (AM and AF directory).

# Import necessary libraries
from pathlib import Path
import shutil
import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS

def phits_generate_inputs():

    total_files = 0

    # Read input template and define paths for phantom and interaction model files
    template = config.INPUT_TEMPLATE_FILE.read_text()
    infl_file_directory = config.INCLUDE_FILES_DIR
    nuclear_data_file = config.NUCLEAR_DATA_FILE
    EGS5_data_directory = config.EGS5_DATA_DIR
    base_output_dir = config.GENERATED_INPUTS_DIR

    # Adult male or female phantom or both
    if config.PHANTOM_INPUT_GENERATION == "AM":
        phantoms = ["AM"]
    elif config.PHANTOM_INPUT_GENERATION == "AF":
        phantoms = ["AF"]
    else:
        phantoms = ["AM", "AF"]

    for phantom in phantoms:

        # Replace placeholders with actual values
        parallelization = config.PARALLELIZATION
        threads = config.THREADS
        maxcas = config.MAXCAS
        maxbch = config.MAXBCH
        source_type = config.PHITS_SOURCE_TYPES
        source_energies = config.SOURCE_ENERGIES
        target_regions = "all"
        sexinfo = config.SEXINFO.get(phantom)

        # Loop through source energies and regions to generate input files
        for energy in source_energies:

            for region, organ_name in SOURCE_ORGANS[phantom].items():

                safe_name = organ_name.replace(",", "").replace(" ", "_")

                phits_output_file = (
                    f"phits_MRCP_{phantom}_source_{safe_name}_{source_type[0]}_energy_{energy}.out"
                )

                deposit_output_file = (
                    f"phits_deposit_MRCP_{phantom}_source_{safe_name}_{source_type[0]}_energy_{energy}.out"
                )

                text = template

                text = text.replace("{{SEXINFO}}", sexinfo)
                text = text.replace("{{PARALLELIZATION}}", parallelization)
                text = text.replace("{{THREADS}}", str(threads))
                text = text.replace("{{MAXCAS}}", f"{maxcas}")
                text = text.replace("{{MAXBCH}}", f"{maxbch}")
                text = text.replace("{{NUCLEARDATAFILE}}", str(nuclear_data_file))
                text = text.replace("{{EGS5DIRECTORY}}", str(EGS5_data_directory))
                text = text.replace("{{SOURCETYPE}}", source_type[0])
                text = text.replace("{{SOURCEREGION}}", str(region))
                text = text.replace("{{SOURCEENERGY}}", f"{energy}")
                text = text.replace("{{SEX}}", str(phantom))
                text = text.replace("{{TARGETREGIONS}}", target_regions)
                text = text.replace("{{PHITSOUTPUTFILE}}", phits_output_file)
                text = text.replace("{{DEPOSITOUTPUTFILE}}", deposit_output_file)

                filename = f"PHITS_MRCP_{phantom}_source_{safe_name}_{source_type[0]}_energy_{energy}.inp"

                phantom_source = Path(infl_file_directory)

                filename = f"PHITS_MRCP_{phantom}_source_{safe_name}_{source_type[0]}_energy_{energy}.inp"
                job_name = filename.removesuffix(".inp")

                # Job directory
                job_dir = base_output_dir / phantom / job_name
                job_dir.mkdir(parents=True, exist_ok=True)

                # Copy phantom files
                phantom_dir = base_output_dir / phantom / "phantoms"
                phantom_dir.mkdir(exist_ok=True)

                for ext in ("cell", "material", "node", "ele"):
                    shutil.copy2(
                        phantom_source / f"MRCP-{phantom}.{ext}",
                        phantom_dir / f"MRCP-{phantom}.{ext}"
                    )

                # Write the input file
                (job_dir / filename).write_text(text)

                total_files += 1

    print(f"\nGenerated PHITS inputs in {base_output_dir}")

    num_phantoms = len(phantoms)
    num_organs_per_phantom = len(SOURCE_ORGANS[phantoms[0]])
    num_source_types = len(config.PHITS_SOURCE_TYPES)
    num_energy_bins = len(config.SOURCE_ENERGIES)

    print(
        f"Generated {total_files} PHITS input file(s) "
        f"({num_phantoms} phantom type(s) × "
        f"{num_organs_per_phantom} source organ(s) × "
        f"{num_source_types} source type(s) × "
        f"{num_energy_bins} source energy bin(s))."
    )