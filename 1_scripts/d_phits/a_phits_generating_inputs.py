# This script generates PHITS input files based on an input template and the b_organ_database.py. 
# It replaces placeholders in the template with actual values for phantom type, source energy, source region, and other parameters, and 
# saves the generated input files in their own directories under the base output directory. 
# It also copies the necessary phantom files to the mother directory (AM and AF directory).

# Import necessary libraries
from pathlib import Path
import shutil
import b_config.a_config as config
from c_database.b_organ_database import SOURCE_ORGANS

def phits_generate_inputs(params):

    total_files = 0

    # Read input template and define paths for phantom and interaction model files
    template = config.INPUT_TEMPLATE_FILE.read_text()
    infl_file_directory = config.INCLUDE_FILES_DIR
    phits_installation_dir = config.PHITS_INSTALLATION_DIR
    base_output_dir = config.GENERATED_INPUTS_DIR

    # T-track configuration
    spongiosa_regions = " ".join(str(region)for region in config.SPONGIOSA_IDS)
    fluence_source_type = config.FLUENCE_SOURCE_TYPES
    energy_bins = config.ENERGY_BINS
    energy_min = config.ENERGY_MIN
    energy_max = config.ENERGY_MAX
    if params["source_type"].lower() in config.FLUENCE_SOURCE_TYPES:
        ttrack_status = ""
    else:
        ttrack_status = "off"

    # Adult male or female phantom or both
    if params["phantom"] == "AM":
        phantoms = ["AM"]
    elif params["phantom"] == "AF":
        phantoms = ["AF"]
    else:
        phantoms = ["AM", "AF"]

    for phantom in phantoms:

        # Replace placeholders with actual values
        parallelization =  params["parallelization"]
        threads =  params["threads"]
        maxcas =  params["maxcas"]
        maxbch =  params["maxbch"]
        source_type = params["source_type"]
        source_energies = params["source_energies"]
        target_regions = "all"
        sexinfo = config.SEXINFO.get(phantom)

        # Loop through source energies and regions to generate input files
        for energy in source_energies:

            for region, organ_name in SOURCE_ORGANS[phantom].items():

                safe_name = organ_name.replace(",", "").replace(" ", "_")

                phits_output_file = (
                    f"phits_MRCP_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out"
                )

                deposit_output_file = (
                    f"phits_deposit_MRCP_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out"
                )

                fluence_output_file = (
                    f"phits_fluence_MRCP_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out"
                )

                text = template

                text = text.replace("{{SEXINFO}}", sexinfo)
                text = text.replace("{{PARALLELIZATION}}", parallelization)
                text = text.replace("{{THREADS}}", str(threads))
                text = text.replace("{{MAXCAS}}", f"{maxcas}")
                text = text.replace("{{MAXBCH}}", f"{maxbch}")
                text = text.replace("{{PHITSINSTALLATION}}", str(phits_installation_dir))
                text = text.replace("{{SOURCETYPE}}", source_type)
                text = text.replace("{{SOURCEREGION}}", str(region))
                text = text.replace("{{SOURCEENERGY}}", f"{energy}")
                text = text.replace("{{SEX}}", str(phantom))
                text = text.replace("{{TARGETREGIONS}}", target_regions)
                text = text.replace("{{PHITSOUTPUTFILE}}", phits_output_file)
                text = text.replace("{{DEPOSITOUTPUTFILE}}", deposit_output_file)
                text = text.replace("{{TTRACKSTATUS}}", ttrack_status)
                text = text.replace("{{FLUENCEOUTPUTFILE}}", fluence_output_file)
                text = text.replace("{{SPONGIOSAREGIONS}}", spongiosa_regions)
                text = text.replace("{{FLUENCESOURCETYPE}}", fluence_source_type[0])
                text = text.replace("{{ENERGYBINS}}", str(energy_bins))
                text = text.replace("{{ENERGYMIN}}", str(energy_min))
                text = text.replace("{{ENERGYMAX}}", str(energy_max))

                filename = f"PHITS_MRCP_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.inp"

                phantom_source = Path(infl_file_directory)

                filename = f"PHITS_MRCP_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.inp"
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
    num_source_types = int(params["source_type"] in config.PHITS_SOURCE_TYPES)
    num_energy_bins = len(params["source_energies"])

    print(
        f"Generated {total_files} PHITS input file(s) "
        f"({num_phantoms} phantom type(s) × "
        f"{num_organs_per_phantom} source organ(s) × "
        f"{num_source_types} source type(s) × "
        f"{num_energy_bins} source energy bin(s))."
    )