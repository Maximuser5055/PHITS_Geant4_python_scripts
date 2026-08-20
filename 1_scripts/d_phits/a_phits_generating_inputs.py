# This script generates PHITS input files based on an input template and the b_organ_database.py. 
# It replaces placeholders in the template with actual values for phantom type, source energy, source region, and other parameters, and 
# saves the generated input files in their own directories under the base output directory. 
# It also copies the necessary phantom files to the mother directory (AM and AF directory).

# Import necessary libraries
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

    # Configs
    sex_info = config.SEXINFO
    phits_source_types = config.PHITS_SOURCE_TYPES
    
    # T-track configuration
    skeletal_regions = " ".join(str(region)for region in config.SKELETAL_IDS)
    fluence_source_type = config.FLUENCE_SOURCE_TYPES
    energy_bins = config.ENERGY_BINS
    energy_min = config.ENERGY_MIN
    energy_max = config.ENERGY_MAX
    if params["source_type"].lower() in config.FLUENCE_SOURCE_TYPES:
        ttrack_status = ""
    else:
        ttrack_status = "off"

    # ============================================================
    # PHANTOM SELECTION
    # ============================================================

    phantom_selection = params["phantom"]

    if phantom_selection == "MRCP_AM":
        phantoms = ["MRCP_AM"]

    elif phantom_selection == "MRCP_AF":
        phantoms = ["MRCP_AF"]

    elif phantom_selection == "MRCP_AF_AM":
        phantoms = ["MRCP_AM", "MRCP_AF"]

    elif phantom_selection == "MFCP_AM":
        phantoms = ["MFCP_AM"]

    elif phantom_selection == "MFCP_AF":
        phantoms = ["MFCP_AF"]

    elif phantom_selection == "MFCP_AF_AM":
        phantoms = ["MFCP_AM", "MFCP_AF"]

    else:
        raise ValueError(
            f"Unknown phantom selection: {phantom_selection}"
        )

    for phantom in phantoms:

        phantom_prefix, AM_or_AF = phantom.split("_")
        
        # Replace placeholders with actual values
        parallelization =  params["parallelization"]
        threads =  params["threads"]
        maxcas =  params["maxcas"]
        maxbch =  params["maxbch"]
        source_type = params["source_type"]
        source_energies = params["source_energies"]
        target_regions = "all"
        sexinfo = sex_info.get(phantom)

        # Loop through source energies and regions to generate input files
        for energy in source_energies:

            for region, organ_name in SOURCE_ORGANS[phantom].items():

                safe_name = organ_name.replace(",", "").replace(" ", "_")

                phits_output_file = (
                    f"phits_{phantom_prefix}_{AM_or_AF}_source_{safe_name}_{source_type}_energy_{energy}.out"
                )

                deposit_output_file = (
                    f"phits_deposit_{phantom_prefix}_{AM_or_AF}_source_{safe_name}_{source_type}_energy_{energy}.out"
                )

                fluence_output_file = (
                    f"phits_fluence_{phantom_prefix}_{AM_or_AF}_source_{safe_name}_{source_type}_energy_{energy}.out"
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
                text = text.replace("{{SEX}}", str(phantom_prefix) + "-" + str(AM_or_AF))
                text = text.replace("{{TARGETREGIONS}}", target_regions)
                text = text.replace("{{PHITSOUTPUTFILE}}", phits_output_file)
                text = text.replace("{{DEPOSITOUTPUTFILE}}", deposit_output_file)
                text = text.replace("{{TTRACKSTATUS}}", ttrack_status)
                text = text.replace("{{FLUENCEOUTPUTFILE}}", fluence_output_file)
                text = text.replace("{{SKELETALREGIONS}}", skeletal_regions)
                text = text.replace("{{FLUENCESOURCETYPE}}", fluence_source_type[0])
                text = text.replace("{{ENERGYBINS}}", str(energy_bins))
                text = text.replace("{{ENERGYMIN}}", str(energy_min))
                text = text.replace("{{ENERGYMAX}}", str(energy_max))

                filename = f"phits_{phantom_prefix}_{AM_or_AF}_source_{safe_name}_{source_type}_energy_{energy}.inp"
                job_name = filename.removesuffix(".inp")

                # Job directory
                job_dir = base_output_dir / phantom / job_name
                job_dir.mkdir(parents=True, exist_ok=True)

                # Copy phantom files
                phantom_dir = base_output_dir / phantom / "phantoms"
                phantom_dir.mkdir(exist_ok=True)

                for ext in ("cell", "material", "node", "ele"):
                    shutil.copy2(
                        infl_file_directory / f"{phantom_prefix}-{AM_or_AF}.{ext}",
                        phantom_dir / f"{phantom_prefix}-{AM_or_AF}.{ext}"
                    )

                # Write the input file
                (job_dir / filename).write_text(text)

                total_files += 1

    print(f"\nGenerated PHITS inputs in {base_output_dir}")

    num_phantoms = len(phantoms)
    num_organs_per_phantom = len(SOURCE_ORGANS[phantoms[0]])
    num_source_types = int(params["source_type"].lower() in phits_source_types)
    num_energy_bins = len(params["source_energies"])

    print(
        f"Generated {total_files} PHITS input file(s) "
        f"({num_phantoms} phantom type(s) × "
        f"{num_organs_per_phantom} source organ(s) × "
        f"{num_source_types} source type(s) × "
        f"{num_energy_bins} source energy bin(s))."
    )