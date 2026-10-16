# This script automates input and macro file generation for Geant4.
# It creates one job per source-organ name and source energy. Multiple organ IDs
# assigned to the same source-organ name are written together to source_id.txt.

import b_config.a_config as config
from b_config.b_phantom_registry import get_phantom_group
from c_database.b_organ_database import SOURCE_ORGANS


def _group_source_regions(source_organs):
    """Group {organ_id: organ_name} entries into {organ_name: [organ_ids]}."""
    grouped = {}
    for organ_id, organ_name in source_organs.items():
        grouped.setdefault(organ_name, []).append(int(organ_id))

    # Stable ordering makes source_id.txt reproducible.
    return {
        organ_name: sorted(organ_ids)
        for organ_name, organ_ids in grouped.items()
    }


def geant4_generate_inputs(params):
    total_jobs = 0

    geant4_source_types = config.GEANT4_SOURCE_TYPES

    # Read templates.
    input_template = (config.GEANT4_BUILD_DIR / "example.in").read_text()
    source_template = (config.GEANT4_BUILD_DIR / "source.mac").read_text()

    base_output_dir = config.GEANT4_GENERATED_INPUTS_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)

    phantom_selection = params["phantom"]
    phantoms = get_phantom_group(phantom_selection)

    threads = params["threads"]
    nps = params["nps"]
    source_type = params["source_type"]
    source_energies = params["source_energies"]

    if source_type not in geant4_source_types:
        raise ValueError(
            f"Unsupported Geant4 source type {source_type!r}. "
            f"Supported source types: {geant4_source_types}"
        )

    for phantom in phantoms:
        if phantom not in SOURCE_ORGANS:
            raise KeyError(
                f"No SOURCE_ORGANS entry found for phantom {phantom!r}."
            )

        grouped_sources = _group_source_regions(SOURCE_ORGANS[phantom])

        for energy in source_energies:
            for organ_name, organ_ids in grouped_sources.items():
                safe_name = organ_name.replace(",", "").replace(" ", "_")

                basename = (
                    f"geant4_{phantom}_"
                    f"source_{safe_name}_"
                    f"{source_type}_"
                    f"energy_{energy}"
                )

                job_dir = base_output_dir / phantom / basename
                job_dir.mkdir(parents=True, exist_ok=True)

                # Generate source.mac.
                mac_text = source_template
                mac_text = mac_text.replace("{{PARTICLE}}", source_type)
                mac_text = mac_text.replace("{{ENERGY}}", str(energy))

                mac_filename = f"{basename}.mac"
                (job_dir / mac_filename).write_text(mac_text)

                # Generate example.in.
                in_text = input_template
                in_text = in_text.replace("{{THREADS}}", str(threads))
                in_text = in_text.replace("{{NPS}}", str(nps))
                in_text = in_text.replace("{{SOURCE_MAC}}", mac_filename)

                in_filename = f"{basename}.in"
                (job_dir / in_filename).write_text(in_text)

                # Multiple IDs belonging to one named source organ are
                # comma-separated, e.g. "9700,9900" for both lungs.
                source_ids_text = ",".join(str(organ_id) for organ_id in organ_ids)
                (job_dir / "source_id.txt").write_text(source_ids_text + "\n")

                total_jobs += 1

    print(f"\\nGenerated Geant4 inputs in {base_output_dir}")

    num_phantoms = len(phantoms)
    num_source_organs = {
        phantom: len(_group_source_regions(SOURCE_ORGANS[phantom]))
        for phantom in phantoms
    }
    num_energies = len(source_energies)
    expected_jobs = sum(num_source_organs.values()) * num_energies

    print(
        f"Generated {total_jobs} Geant4 job(s) "
        f"({total_jobs} .in files, {total_jobs} .mac files, "
        f"and {total_jobs} source_id.txt files)."
    )
    print(
        f"Configuration: {num_phantoms} phantom type(s), "
        f"{sum(num_source_organs.values())} source-organ groups across selected "
        f"phantoms, {num_energies} energy bin(s)."
    )
