# This script parses the .cell and .material files to extract organ
# information and generates a Python database file containing organ IDs,
# names, densities, volumes, and masses for the selected phantoms.

# import necessary libraries
import re
import pandas as pd
import b_config.a_config as config
from b_config.b_phantom_registry import PHANTOMS, get_phantom_group


# Organ database file path
database_file_path = config.DATABASE_FILE


# Define the regex pattern to match the cell data in the .cell file
cell_pattern = re.compile(
    r"^\s*(\d+)\s+"          # Cell number
    r"(\d+)\s+"              # Organ/material ID
    r"(-?\d+\.\d+)\s+"       # Density
    r"-?\d+\s+"              # Surface (ignore)
    r"u=\d+\s+"              # Universe (ignore)
    r"VOL=(\d+\.\d+)"        # Volume
)


# Material-file header:
# $ Organ_name density g/cm3
material_header_pattern = re.compile(
    r"^\s*\$\s*(.*?)\s+"
    r"[-+]?\d+(?:\.\d+)?\s+g/cm3\s*$"
)

# Material definition:
# MAT[organ_id]
material_id_pattern = re.compile(
    r"^\s*MAT\[(\d+)\]\s*$"
)


def parse_material(material_file_path):
    """
    Parse the PHITS .material file and return an organ-ID-to-name mapping.

    The .material files contain blocks of the form:

        $ Liver 1.060 g/cm3
        MAT[9500]
        ...

    Therefore, the organ name is taken directly from the material-file
    comment rather than from ORGAN_ID_CSV.
    """
    material_names = {}
    pending_name = None

    with open(material_file_path, encoding="utf-8") as f:

        for line in f:

            header_match = material_header_pattern.match(line)

            if header_match:
                pending_name = header_match.group(1).strip()
                continue

            material_match = material_id_pattern.match(line)

            if material_match:

                if pending_name is None:
                    raise ValueError(
                        f"Found {line.strip()} without a preceding "
                        f"material-name header in {material_file_path}."
                    )

                organ_id = int(material_match.group(1))
                material_names[organ_id] = pending_name
                pending_name = None

    return material_names


def parse_cell_csv_inputs(params):

    source_organs_csv_file_path = params["source_csv"]

    phantom_selection = params["phantom"]

    phantom_codes = get_phantom_group(phantom_selection)

    def parse_cell(cell_file_path):

        # Initialize an empty dictionary to store organ data
        organs = {}

        with open(cell_file_path, encoding="utf-8") as f:

            for line in f:

                match = cell_pattern.match(line)

                if match:

                    organ_id = int(match.group(2))
                    density = abs(float(match.group(3)))
                    volume = float(match.group(4))

                    mass = density * volume

                    organs[organ_id] = {
                        "organ_id": organ_id,
                        "density": density,
                        "volume": volume,
                        "mass": mass,
                    }

        return organs

    phantoms = []

    for phantom_code in phantom_codes:

        phantom = PHANTOMS[phantom_code]

        # Parse organ names directly from this phantom's material file.
        material_names = parse_material(phantom.material_file)

        # Parse organ IDs, density, volume, and mass from the cell file.
        organs = parse_cell(phantom.cell_file)

        for organ_id, organ in organs.items():

            organ_name = material_names.get(organ_id)

            if organ_name is None:
                raise ValueError(
                    f"Organ ID {organ_id} exists in "
                    f"{phantom.cell_file.name} but has no matching "
                    f"material name in {phantom.material_file.name}."
                )

            organ["name"] = organ_name

        phantoms.append((phantom_code, organs))

    # Reading the source organs
    organ_groups = pd.read_csv(source_organs_csv_file_path)

    main_source_organs = (
        organ_groups["source_organ_ID"]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    # -------------------------------------------------------------------------
    # Write the organ database
    # -------------------------------------------------------------------------

    with open(database_file_path, "w", encoding="utf-8") as f:

        # =====================================================================
        # ORGANS
        # =====================================================================

        f.write("ORGANS = {\n")

        for phantom_code, organs in phantoms:

            f.write(f'    "{phantom_code}": {{\n')

            for organ_id, organ in organs.items():

                f.write(f"        {organ_id}: {{\n")
                f.write(f'            "organ_id": {organ["organ_id"]},\n')
                f.write(f'            "name": "{organ["name"]}",\n')
                f.write(f'            "density": {organ["density"]},\n')
                f.write(f'            "volume": {organ["volume"]},\n')
                f.write(f'            "mass": {organ["mass"]},\n')
                f.write("        },\n\n")

            f.write("    },\n\n")

        f.write("}\n\n")

        # =====================================================================
        # SOURCE ORGANS
        # =====================================================================

        f.write("SOURCE_ORGANS = {\n")

        for phantom_code, organs in phantoms:

            f.write(f'    "{phantom_code}": {{\n')

            for organ_id in main_source_organs:

                organ = organs.get(organ_id)

                if organ is None:

                    print(
                        f"Warning: Source organ ID {organ_id} does not "
                        f"exist in phantom {phantom_code}. Skipping."
                    )

                    continue

                f.write(
                    f'        {organ_id}: "{organ["name"]}",\n'
                )

            f.write("    },\n\n")

        f.write("}\n\n")
