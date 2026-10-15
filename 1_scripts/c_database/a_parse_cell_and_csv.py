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
    r"(-?\d+(?:\.\d+)?)\s+" # Density (integer or decimal)
    r"-?\d+\s+"              # Surface (ignore)
    r"u=\d+\s+"              # Universe (ignore)
    r"VOL=(\d+(?:\.\d+)?)"  # Volume (integer or decimal)
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

    # Read source organs directly from source_organs.csv.
    # Each row can contain multiple IDs separated by underscores, and the
    # CSV acronym is used as the source-organ name.
    organ_groups = pd.read_csv(source_organs_csv_file_path)

    required_columns = {"Acronym", "ID number(s)"}
    missing_columns = required_columns.difference(organ_groups.columns)
    if missing_columns:
        raise ValueError(
            "The source-organ CSV is missing required column(s): "
            + ", ".join(sorted(missing_columns))
        )

    # Preserve CSV order and group IDs that share the same acronym.
    grouped_source_organs = {}
    for _, row in organ_groups.iterrows():
        if pd.isna(row["Acronym"]) or pd.isna(row["ID number(s)"]):
            continue

        acronym = str(row["Acronym"]).strip()
        organ_ids = [
            int(value.strip())
            for value in str(row["ID number(s)"]).split("_")
            if value.strip()
        ]

        if not acronym:
            continue

        grouped_source_organs.setdefault(acronym, [])
        for organ_id in organ_ids:
            if organ_id not in grouped_source_organs[acronym]:
                grouped_source_organs[acronym].append(organ_id)

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

            for acronym, organ_ids in grouped_source_organs.items():
                # Keep only IDs that exist in this phantom.
                valid_ids = []
                for organ_id in organ_ids:
                    if organ_id not in organs:
                        print(
                            f"Warning: Source organ ID {organ_id} does not "
                            f"exist in phantom {phantom_code}. Skipping."
                        )
                        continue
                    valid_ids.append(organ_id)

                if not valid_ids:
                    continue

                # Use dict.fromkeys so multiple IDs sharing an acronym can
                # be written compactly while remaining valid Python.
                ids_repr = repr(valid_ids)
                acronym_repr = repr(acronym)
                f.write(
                    f"        **dict.fromkeys({ids_repr}, {acronym_repr}),\n"
                )

            f.write("    },\n\n")

        f.write("}\n\n")
