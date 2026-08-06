# This script parses the .cell files to extract organ information and 
# generates a Python database file containing organ IDs, names, densities, volumes, and masses for both male and female phantoms.

# import necessary libraries
import re
import pandas as pd
import b_config.a_config as config

# Define path parameters
male_cell_file_path = config.CELL_FILES["AM"]
female_cell_file_path = config.CELL_FILES["AF"]
csv_file_path = config.ORGAN_ID_CSV
source_organs_csv_file_path = config.SOURCE_CSV

# Organ database file path
database_file_path = config.DATABASE_FILE

# Define the regex pattern to match the cell data in the .cell file
pattern = re.compile(
    r"^\s*(\d+)\s+"          # Cell number
    r"(\d+)\s+"              # Organ/material ID
    r"(-?\d+\.\d+)\s+"       # Density
    r"-?\d+\s+"              # Surface (ignore)
    r"u=\d+\s+"              # Universe (ignore)
    r"VOL=(\d+\.\d+)"        # Volume
)
def parse_cell_csv_inputs():
    def parse_cell(cell_file_path):
        # Initialize an empty dictionary to store organ data
        organs = {}

        with open(cell_file_path) as f:

            for line in f:

                match = pattern.match(line)

                if match:

                    organ_id = int(match.group(1))
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

    # Use parse function to extract organ data from male and female cell files, and then map organ IDs to names using the CSV file.
    male_organs = parse_cell(male_cell_file_path)
    female_organs = parse_cell(female_cell_file_path)

    names = pd.read_csv(csv_file_path)

    name_dict = dict(zip(names.organ_id, names.name))

    for organs in (male_organs, female_organs):

        for organ_id, organ in organs.items():
            organ["name"] = (
                name_dict.get(organ_id, "Unknown")
                .replace(", ", "_")
            )

    # Reading the source organs
    organ_groups = pd.read_csv(source_organs_csv_file_path)

    main_source_organs = (
        organ_groups["source_organ_ID"]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    # Write the organ data to a Python file in the specified format
    with open(database_file_path, "w") as f:

        f.write("ORGANS = {\n")

        for phantom_name, organs in [
            ("AM", male_organs),
            ("AF", female_organs)
        ]:

            f.write(f'    "{phantom_name}": {{\n')

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
        
        # -------------------------------------------------------------------------
        # Source organs
        # -------------------------------------------------------------------------

        f.write("SOURCE_ORGANS = {\n")

        for phantom_name, organs in [
            ("AM", male_organs),
            ("AF", female_organs),
        ]:

            f.write(f'    "{phantom_name}": {{\n')

            for organ_id in main_source_organs:

                organ = organs.get(organ_id)

                if organ is None:
                    organ_name = name_dict.get(organ_id, "Unknown")

                    print(
                        f"Warning: Source organ ID {organ_id} ({organ_name})"
                        f"does not exist in phantom {phantom_name}. Skipping."
                    )
                    continue

                f.write(
                    f'        {organ_id}: "{organ["name"]}",\n'
                )

            f.write("    },\n\n")

        f.write("}\n\n")