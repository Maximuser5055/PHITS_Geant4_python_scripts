from pathlib import Path
import re
import numpy as np
import pandas as pd

# ============================================================
# USER CONFIGURATION
# ============================================================

# ------------------------------------------------------------
# Reference phantom group
# ------------------------------------------------------------
REFERENCE_GROUP = "ICRP145_RESIZED"

# ------------------------------------------------------------
# Coordinate unit in NODE files
# ------------------------------------------------------------
NODE_UNIT = "cm"


# ------------------------------------------------------------
# Input CSV files
# ------------------------------------------------------------

ORGAN_ID_NAMES_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/organ_ID_names.csv")
TARGET_REGIONS_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/target_regions_Filipino.csv")

# ------------------------------------------------------------
# Phantom groups
# ------------------------------------------------------------

PHANTOM_GROUPS = {

    "ICRP145": {

        "Male": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AM.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AM.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AM.material"),
        },

        "Female": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AF.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AF.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/MRCP-AF.material"),
        },
    },


    "ICRP145_RESIZED": {

        "Male": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/M_H165W65.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/M_H165W65.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/M_H165W65.material"),
        },

        "Female": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/F_H150W55.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/F_H150W55.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/F_H150W55.material"),
        },
    },


    "FILIPINO": {

        "Male": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AM.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AM.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AM.material"),
        },

        "Female": {
            "node": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AF.node"),
            "ele": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AF.ele"),
            "material": Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/2_phits/phantoms/Filipino-MRCP-AF.material"),
        },
    },
}


# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

OUTPUT_ORGAN_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/phantom_comparison.csv")
OUTPUT_COMPOUND_FILE = Path(r"/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/compound_organ_weights.csv")


# ============================================================
# GENERAL UTILITIES
# ============================================================

def convert_to_cm(value, unit):

    if unit.lower() == "cm":
        return value

    if unit.lower() == "mm":
        return value / 10.0

    if unit.lower() == "m":
        return value * 100.0

    raise ValueError(
        f"Unknown NODE unit: {unit}"
    )

# ============================================================
# ORGAN ID NAMES
# ============================================================

def load_organ_names(filename):

    df = pd.read_csv(filename)

    df.columns = [
        str(c).strip().lstrip("\ufeff")
        for c in df.columns
    ]

    id_column = "organ_id"
    name_column = "name"

    if id_column not in df.columns:
        raise ValueError(
            f"Could not find '{id_column}' in {filename}"
        )

    if name_column not in df.columns:
        raise ValueError(
            f"Could not find '{name_column}' in {filename}"
        )

    return {
        int(row[id_column]): str(row[name_column])
        for _, row in df.iterrows()
    }


# ============================================================
# NODE FILE
# ============================================================

def read_node_file(filename):

    print(f"Reading NODE: {filename}")

    with open(filename, "r") as f:

        # Find header
        n_elements = 0
        n_elements = None
        nodes_per_element = None
        n_nodes = None
        dimension = None

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            n_nodes = int(parts[0])
            dimension = int(parts[1])

            break

        if n_nodes is None or dimension != 3:
            raise ValueError(
                f"{filename} is not a valid 3D NODE file."
            )

        coordinates = np.empty(
            (n_nodes, 3),
            dtype=np.float64
        )

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            node_id = int(parts[0])

            coordinates[node_id] = [
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            ]

    # Convert to cm
    if NODE_UNIT == "mm":
        coordinates /= 10.0

    elif NODE_UNIT == "m":
        coordinates *= 100.0

    elif NODE_UNIT != "cm":
        raise ValueError(
            "NODE_UNIT must be mm, cm, or m."
        )

    return coordinates


# ============================================================
# MATERIAL FILE
# ============================================================

def read_material_densities(filename):

    print(f"Reading MATERIAL: {filename}")

    densities = {}

    current_density = None

    density_pattern = re.compile(
        r"^\s*\$\s+.*?([0-9]+(?:\.[0-9]+)?)\s+g/cm3",
        re.IGNORECASE
    )

    mat_pattern = re.compile(
        r"^\s*MAT\[(\d+)\]"
    )

    with open(filename, "r") as f:

        for line in f:

            density_match = density_pattern.match(line)

            if density_match:

                current_density = float(
                    density_match.group(1)
                )

                continue

            mat_match = mat_pattern.match(line)

            if mat_match:

                if current_density is None:

                    raise ValueError(
                        f"MAT found without density:\n{line}"
                    )

                material_id = int(
                    mat_match.group(1)
                )

                densities[material_id] = current_density

                current_density = None

    return densities


# ============================================================
# CALCULATE HEIGHT
# ============================================================

def calculate_height(coordinates):

    x_min = np.min(coordinates[:, 0])
    x_max = np.max(coordinates[:, 0])

    y_min = np.min(coordinates[:, 1])
    y_max = np.max(coordinates[:, 1])

    z_min = np.min(coordinates[:, 2])
    z_max = np.max(coordinates[:, 2])

    return {
        "height_cm": z_max - z_min,
        "width_cm": x_max - x_min,
        "depth_cm": y_max - y_min,
    }


# ============================================================
# CALCULATE ORGAN VOLUMES AND MASSES
# ============================================================

def calculate_organ_masses(
    ele_file,
    coordinates,
    densities,
    chunk_size=200_000,
):

    print(f"Reading ELE: {ele_file}")

    organ_volumes = {}
    organ_masses = {}

    with open(ele_file, "r") as f:

        # ----------------------------------------------------
        # Read header
        # ----------------------------------------------------

        n_elements = None
        nodes_per_element = None

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            header = line.split()

            n_elements = int(header[0])
            nodes_per_element = int(header[1])

            break

        if (
            n_elements is None
            or nodes_per_element is None
            or nodes_per_element != 4
        ):

            raise ValueError(
                "Invalid ELE header; this script expects "
                "tetrahedral elements."
            )

        # ----------------------------------------------------
        # Process chunks
        # ----------------------------------------------------

        chunk = []

        def process_chunk(lines):

            if not lines:
                return

            data = np.fromstring(
                "\n".join(lines),
                sep=" ",
                dtype=np.float64,
            )

            data = data.reshape(-1, 6)

            node_ids = data[:, 1:5].astype(
                np.int64
            )

            material_ids = data[:, 5].astype(
                np.int64
            )

            # Coordinates of tetrahedron vertices

            a = coordinates[node_ids[:, 0]]
            b = coordinates[node_ids[:, 1]]
            c = coordinates[node_ids[:, 2]]
            d = coordinates[node_ids[:, 3]]

            ab = b - a
            ac = c - a
            ad = d - a

            volumes = np.abs(
                np.einsum(
                    "ij,ij->i",
                    ab,
                    np.cross(ac, ad),
                )
            ) / 6.0

            # ----------------------------------------------
            # Accumulate each organ separately
            # ----------------------------------------------

            unique_materials = np.unique(
                material_ids
            )

            for material_id in unique_materials:

                mask = (
                    material_ids == material_id
                )

                material_volume = np.sum(
                    volumes[mask]
                )

                material_id = int(
                    material_id
                )

                organ_volumes[material_id] = (
                    organ_volumes.get(
                        material_id,
                        0.0
                    )
                    + material_volume
                )

                if material_id not in densities:

                    raise ValueError(
                        f"Material ID {material_id} "
                        f"from ELE was not found "
                        f"in MATERIAL file."
                    )

                organ_masses[material_id] = (
                    organ_masses.get(
                        material_id,
                        0.0
                    )
                    + material_volume
                    * densities[material_id]
                )

        # ----------------------------------------------------
        # Read elements
        # ----------------------------------------------------

        count = 0

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            chunk.append(line)

            if len(chunk) >= chunk_size:

                process_chunk(chunk)

                count += len(chunk)

                chunk.clear()

                if count % 1_000_000 < chunk_size:

                    print(
                        f"  Processed "
                        f"{count:,} / "
                        f"{n_elements:,} elements"
                    )

        # Remaining elements
        process_chunk(chunk)

    return organ_volumes, organ_masses


# ============================================================
# PROCESS ONE PHANTOM
# ============================================================

def process_phantom(
    group_name,
    sex,
    config,
):

    print()
    print("=" * 60)
    print(f"{group_name} - {sex}")
    print("=" * 60)

    coordinates = read_node_file(
        config["node"]
    )

    dimensions = calculate_height(
        coordinates
    )

    densities = read_material_densities(
        config["material"]
    )

    organ_volumes, organ_masses = (
        calculate_organ_masses(
            config["ele"],
            coordinates,
            densities,
        )
    )

    total_mass = sum(
        organ_masses.values()
    )

    return {
        "group": group_name,
        "sex": sex,
        "height_cm": dimensions["height_cm"],
        "width_cm": dimensions["width_cm"],
        "depth_cm": dimensions["depth_cm"],
        "total_mass_g": total_mass,
        "organ_masses": organ_masses,
        "organ_volumes": organ_volumes,
    }


# ============================================================
# PERCENT ERROR
# ============================================================

def percent_error(value, reference):

    if reference == 0:

        return np.nan

    return (
        (value - reference)
        / reference
        * 100.0
    )


# ============================================================
# ORGAN COMPARISON TABLE
# ============================================================

def create_organ_comparison(
    results,
    organ_names,
):

    # --------------------------------------------------------
    # Group results by sex
    # --------------------------------------------------------

    reference_results = {
        r["sex"]: r
        for r in results
        if r["group"] == REFERENCE_GROUP
    }

    rows = []

    for result in results:

        sex = result["sex"]

        reference = reference_results.get(
            sex
        )

        if reference is None:

            raise ValueError(
                f"No reference {sex} phantom "
                f"found in {REFERENCE_GROUP}."
            )

        # ----------------------------------------------------
        # Height row
        # ----------------------------------------------------

        height_error = percent_error(
            result["height_cm"],
            reference["height_cm"],
        )

        rows.append({
            "Phantom Group": result["group"],
            "Sex": sex,
            "Quantity": "Height",
            "Organ ID": "",
            "Organ/tissue": "",
            "Reference (cm or g)": reference["height_cm"],
            "Value (cm or g)": result["height_cm"],
            "% Error": height_error,
        })

        # ----------------------------------------------------
        # Total mass
        # ----------------------------------------------------

        total_mass_error = percent_error(
            result["total_mass_g"],
            reference["total_mass_g"],
        )

        rows.append({
            "Phantom Group": result["group"],
            "Sex": sex,
            "Quantity": "Total Body Mass",
            "Organ ID": "",
            "Organ/tissue": "",
            "Reference (cm or g)": reference["total_mass_g"],
            "Value (cm or g)": result["total_mass_g"],
            "% Error":  total_mass_error,
        })

        # ----------------------------------------------------
        # Individual organs
        # ----------------------------------------------------

        all_ids = set(
            reference["organ_masses"]
        ) | set(
            result["organ_masses"]
        )

        for organ_id in sorted(all_ids):

            reference_mass = (
                reference["organ_masses"]
                .get(organ_id, 0.0)
            )

            phantom_mass = (
                result["organ_masses"]
                .get(organ_id, 0.0)
            )

            error = percent_error(
                phantom_mass,
                reference_mass,
            )

            rows.append({
                "Phantom Group": result["group"],
                "Sex": sex,
                "Quantity": "Organ",
                "Organ ID": organ_id,
                "Organ/tissue": organ_names.get(
                    organ_id,
                    f"Unknown ({organ_id})"
                ),
                "Reference (cm or g)": reference_mass,
                "Value (cm or g)": phantom_mass,
                "% Error": error,
            })

    return pd.DataFrame(rows)


# ============================================================
# TARGET REGION CSV
# ============================================================

def load_compound_regions(filename):

    df = pd.read_csv(filename)

    df.columns = [
        str(c).strip().lstrip("\ufeff")
        for c in df.columns
    ]

    required = [
        "Target region",
        "Acronym",
        "ID number(s)",
    ]

    for column in required:

        if column not in df.columns:

            raise ValueError(
                f"Column '{column}' not found "
                f"in {filename}"
            )

    compound_regions = []

    for _, row in df.iterrows():

        ids = [
            int(x)
            for x in str(
                row["ID number(s)"]
            ).split("_")
            if x.strip()
        ]

        compound_regions.append({
            "name": row["Target region"],
            "acronym": row["Acronym"],
            "ids": ids,
        })

    return compound_regions


# ============================================================
# COMPOUND ORGAN MASS TABLE
# ============================================================

def create_compound_comparison(
    results,
    compound_regions,
):

    reference_results = {
        r["sex"]: r
        for r in results
        if r["group"] == REFERENCE_GROUP
    }

    rows = []

    for result in results:

        sex = result["sex"]

        reference = reference_results.get(
            sex
        )

        if reference is None:

            continue

        for compound in compound_regions:

            ids = compound["ids"]

            # ----------------------------------------------
            # Reference compound mass
            # ----------------------------------------------

            reference_mass = sum(
                reference["organ_masses"].get(
                    organ_id,
                    0.0
                )
                for organ_id in ids
            )

            # ----------------------------------------------
            # Phantom compound mass
            # ----------------------------------------------

            phantom_mass = sum(
                result["organ_masses"].get(
                    organ_id,
                    0.0
                )
                for organ_id in ids
            )

            error = percent_error(
                phantom_mass,
                reference_mass,
            )

            rows.append({
                "Phantom Group": result["group"],
                "Sex": sex,
                "Compound Organ": compound["name"],
                "Acronym": compound["acronym"],
                "Constituent Organ IDs": "_".join(
                    str(x)
                    for x in ids
                ),
                "Reference Mass (g)": (
                    reference_mass
                ),
                "Phantom Mass (g)": (
                    phantom_mass
                ),
                "% Error": error,
            })

    return pd.DataFrame(rows)


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 70)
    print("PHANTOM GROUP COMPARISON")
    print("=" * 70)

    if REFERENCE_GROUP not in PHANTOM_GROUPS:

        raise ValueError(
            f"Reference group '{REFERENCE_GROUP}' "
            f"is not defined."
        )

    # --------------------------------------------------------
    # Load organ names
    # --------------------------------------------------------

    organ_names = load_organ_names(
        ORGAN_ID_NAMES_FILE
    )

    # --------------------------------------------------------
    # Process all phantoms
    # --------------------------------------------------------

    results = []

    for group_name, group_config in (
        PHANTOM_GROUPS.items()
    ):

        for sex, config in group_config.items():

            result = process_phantom(
                group_name,
                sex,
                config,
            )

            results.append(result)

    # --------------------------------------------------------
    # Individual organ comparison
    # --------------------------------------------------------

    organ_df = create_organ_comparison(
        results,
        organ_names,
    )

    organ_df.to_csv(
        OUTPUT_ORGAN_FILE,
        index=False,
    )

    # --------------------------------------------------------
    # Compound organs
    # --------------------------------------------------------

    compound_regions = (
        load_compound_regions(
            TARGET_REGIONS_FILE
        )
    )

    compound_df = (
        create_compound_comparison(
            results,
            compound_regions,
        )
    )

    compound_df.to_csv(
        OUTPUT_COMPOUND_FILE,
        index=False,
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("OUTPUT")
    print("=" * 70)

    print()
    print(
        f"Organ comparison:"
        f"\n  {OUTPUT_ORGAN_FILE}"
    )

    print(
        f"\nCompound-organ comparison:"
        f"\n  {OUTPUT_COMPOUND_FILE}"
    )

    print()
    print(
        f"Reference group: {REFERENCE_GROUP}"
    )

    print()
    print("Done.")


if __name__ == "__main__":
    main()