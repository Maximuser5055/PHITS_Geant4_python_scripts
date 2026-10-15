# This script generates PHITS input files based on an input template and the b_organ_database.py. 
# It replaces placeholders in the template with actual values for phantom type, source energy, source region, and other parameters, and 
# saves the generated input files in their own directories under the base output directory. 
# It also copies the necessary phantom files to the mother directory (AM and AF directory).

# Import necessary libraries
import re
import shutil

import b_config.a_config as config
from b_config.b_phantom_registry import (
    PHANTOMS,
    get_phantom_group,
    get_skeletal_ids
)
from c_database.b_organ_database import SOURCE_ORGANS

def parse_cell_organ_masses(cell_file_path):
    """Return organ masses (g) from a phantom cell file.

    Cell-file volume is in cm^3 and density is in g/cm^3, so mass = |density| * volume.
    The cell-file path is obtained from the phantom registry.
    """
    cell_pattern = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(-?\d+(?:\.\d+)?)\s+"
        r"-?\d+\s+u=\d+\s+VOL=(\d+(?:\.\d+)?)"
    )
    masses = {}
    with open(cell_file_path, encoding="utf-8") as cell_file:
        for line in cell_file:
            match = cell_pattern.match(line)
            if match:
                organ_id = int(match.group(2))
                density = abs(float(match.group(3)))
                volume = float(match.group(4))
                masses[organ_id] = masses.get(organ_id, 0.0) + density * volume
    return masses


def group_source_organs(source_organs):
    """Group organ IDs by acronym, preserving their first-seen order."""
    grouped = {}
    for organ_id, acronym in source_organs.items():
        grouped.setdefault(acronym, []).append(int(organ_id))
    return grouped


def build_multi_source_block(source_ids, masses, source_type, energy):
    """Build a PHITS multi-source definition weighted by source-region mass."""
    missing = [organ_id for organ_id in source_ids if organ_id not in masses]
    if missing:
        raise ValueError(
            f"Source organ ID(s) {missing} were not found in the registered cell file."
        )

    total_mass = sum(masses[organ_id] for organ_id in source_ids)
    if total_mass <= 0:
        raise ValueError(
            f"Total source mass must be positive for organ IDs {source_ids}; "
            f"got {total_mass}."
        )

    lines = [
        "[ Source ]",
        "  totfact = 1.00         # global factor",
    ]
    for organ_id in source_ids:
        weight = masses[organ_id] / total_mass
        lines.extend([
            f"  <source> = {weight:.12g}  # mass fraction for organ ID {organ_id}",
            "    s-type = 24",
            f"      proj = {source_type}",
            f"    tetreg = {organ_id}",
            f"        e0 = {energy}",
            "       dir = all",
        ])
    return "\n".join(lines) + "\n"


def phits_generate_inputs(params):

    total_files = 0

    # Read input template and define paths for phantom and interaction model files
    template = config.INPUT_TEMPLATE_FILE.read_text()
    phits_installation_dir = config.PHITS_INSTALLATION_DIR
    base_output_dir = config.GENERATED_INPUTS_DIR

    # Configs
    phits_source_types = config.PHITS_SOURCE_TYPES
    
    # T-track configuration
    skeletal_ids = get_skeletal_ids(params["phantom"])
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

    phantoms = get_phantom_group(phantom_selection)

    skeletal_regions = " ".join(str(region) for region in skeletal_ids)
    
    for phantom in phantoms:

        phantom_spec = PHANTOMS[phantom]
        source_masses = parse_cell_organ_masses(phantom_spec.cell_file)
        grouped_sources = group_source_organs(SOURCE_ORGANS[phantom])
        
        # Replace placeholders with actual values
        parallelization =  params["parallelization"]
        threads =  params["threads"]
        maxcas =  params["maxcas"]
        maxbch =  params["maxbch"]
        source_type = params["source_type"]
        source_energies = params["source_energies"]
        target_regions = "all"
        namesexinfo = phantom_spec.display_name

        # Loop through source energies and regions to generate input files
        for energy in source_energies:

            for organ_name, source_ids in grouped_sources.items():

                safe_name = organ_name.replace(",", "").replace(" ", "_")
                phits_output_file = (f"phits_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out")
                deposit_output_file = (f"phits_deposit_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out")
                fluence_output_file = (f"phits_fluence_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.out")

                text = template

                text = text.replace("{{NAMESEXINFO}}", namesexinfo)
                text = text.replace("{{PARALLELIZATION}}", parallelization)
                text = text.replace("{{THREADS}}", str(threads))
                text = text.replace("{{MAXCAS}}", f"{maxcas}")
                text = text.replace("{{MAXBCH}}", f"{maxbch}")
                text = text.replace("{{PHITSINSTALLATION}}", str(phits_installation_dir))
                text = text.replace("{{SOURCETYPE}}", source_type)
                text = text.replace("{{SOURCEREGION}}", str(source_ids[0]))
                text = text.replace("{{SOURCEENERGY}}", f"{energy}")
                text = text.replace("{{MATERIAL_FILE}}", phantom_spec.material_file.name)
                text = text.replace("{{CELL_FILE}}", phantom_spec.cell_file.name)
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

                # Keep the existing [Source] template for a single-ID source.
                # For grouped organs, replace it with PHITS multi-source blocks;
                # each source's relative weight is its mass / total grouped mass.
                if len(source_ids) > 1:
                    source_block = build_multi_source_block(
                        source_ids, source_masses, source_type, energy
                    )
                    text = re.sub(
                        r"(?ms)^\[ Source \].*?(?=^\[ Material \])",
                        source_block + "\n",
                        text,
                        count=1,
                    )

                filename = (f"phits_{phantom}_source_{safe_name}_{source_type}_energy_{energy}.inp")
                job_name = filename.removesuffix(".inp")

                # Job directory
                job_dir = base_output_dir / phantom / job_name
                job_dir.mkdir(parents=True, exist_ok=True)

                # Copy phantom files
                phantom_dir = base_output_dir / phantom / "phantoms"
                phantom_dir.mkdir(exist_ok=True)

                for phantom_file in (
                    phantom_spec.cell_file,
                    phantom_spec.material_file,
                    phantom_spec.node_file,
                    phantom_spec.element_file,
                ):
                    shutil.copy2(phantom_file, phantom_dir / phantom_file.name)

                # Write the input file
                (job_dir / filename).write_text(text)

                total_files += 1

    print(f"\nGenerated PHITS inputs in {base_output_dir}")

    num_phantoms = len(phantoms)
    num_organs_per_phantom = len(group_source_organs(SOURCE_ORGANS[phantoms[0]]))
    num_source_types = int(params["source_type"].lower() in phits_source_types)
    num_energy_bins = len(params["source_energies"])

    print(
        f"Generated {total_files} PHITS input file(s) "
        f"({num_phantoms} phantom type(s) × "
        f"{num_organs_per_phantom} source organ(s) × "
        f"{num_source_types} source type(s) × "
        f"{num_energy_bins} source energy bin(s))."
    )