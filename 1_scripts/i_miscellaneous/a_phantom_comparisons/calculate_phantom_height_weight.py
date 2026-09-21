'''
Use: 

python 1_scripts/i_miscellaneous/a_phantom_comparisons/calculate_phantom_height_weight.py \
    --node 2_phits/phantoms/M_H165W65.node \
    --ele 2_phits/phantoms/M_H165W65.ele \
    --material 2_phits/phantoms/M_H165W65.material \
    --coordinate-unit cm
'''

import argparse
import re
import numpy as np


# ============================================================
# Read NODE file
# ============================================================

def read_node_file(node_file, coordinate_unit="cm"):
    """
    Read a TetGen-style .node file.

    Expected format:
        <number of nodes> <dimension> <attributes> <boundary markers>
        node_id x y z

    Returns
    -------
    coordinates : numpy.ndarray
        Shape (n_nodes, 3)
        Coordinates in cm.
    """

    print(f"Reading NODE file: {node_file}")

    with open(node_file, "r") as f:
        # Skip blank/comment lines
        dimension = None
        n_nodes = None
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            header = line.split()
            n_nodes = int(header[0])
            dimension = int(header[1])
            break

        if dimension is None or n_nodes is None:
            raise ValueError("NODE file does not contain a valid header.")

        if dimension != 3:
            raise ValueError(
                f"Expected 3D coordinates, but NODE file has dimension {dimension}."
            )

        coordinates = np.empty((n_nodes, 3), dtype=np.float64)

        count = 0

        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            node_id = int(parts[0])

            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])

            coordinates[node_id] = [x, y, z]

            count += 1

    if count != n_nodes:
        raise ValueError(
            f"Expected {n_nodes} nodes, but read {count} nodes."
        )

    # Convert coordinates to cm if necessary
    if coordinate_unit.lower() == "mm":
        coordinates /= 10.0

    elif coordinate_unit.lower() == "m":
        coordinates *= 100.0

    elif coordinate_unit.lower() != "cm":
        raise ValueError(
            "coordinate_unit must be 'mm', 'cm', or 'm'."
        )

    print(f"  Nodes read: {n_nodes:,}")

    return coordinates


# ============================================================
# Read MATERIAL file
# ============================================================

def read_material_densities(material_file):
    """
    Reads material densities from the .material file.

    Example:

        $ Adrenal_left 1 g/cm3
        MAT[100]
            1000 -0.104
            ...

    Returns
    -------
    densities : dict
        {material_id: density_g_per_cm3}
    """

    print(f"Reading MATERIAL file: {material_file}")

    densities = {}

    pending_density = None

    density_pattern = re.compile(
        r"^\s*\$\s+.*?([0-9]+(?:\.[0-9]+)?)\s+g/cm3\s*$",
        re.IGNORECASE
    )

    mat_pattern = re.compile(
        r"^\s*MAT\[(\d+)\]"
    )

    with open(material_file, "r") as f:

        for line in f:

            # Look for:
            # $ Tissue_name 1.06 g/cm3
            density_match = density_pattern.match(line)

            if density_match:
                pending_density = float(density_match.group(1))
                continue

            # Look for:
            # MAT[100]
            mat_match = mat_pattern.match(line)

            if mat_match and pending_density is not None:

                material_id = int(mat_match.group(1))

                densities[material_id] = pending_density

                pending_density = None

    print(f"  Materials found: {len(densities):,}")

    return densities


# ============================================================
# Calculate phantom dimensions
# ============================================================

def calculate_height(coordinates):

    x_min = np.min(coordinates[:, 0])
    x_max = np.max(coordinates[:, 0])

    y_min = np.min(coordinates[:, 1])
    y_max = np.max(coordinates[:, 1])

    z_min = np.min(coordinates[:, 2])
    z_max = np.max(coordinates[:, 2])

    width = x_max - x_min
    depth = y_max - y_min
    height = z_max - z_min

    return {
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "z_min": z_min,
        "z_max": z_max,
        "width": width,
        "depth": depth,
        "height": height,
    }


# ============================================================
# Calculate volume and mass
# ============================================================

def calculate_volume_and_mass(
    ele_file,
    coordinates,
    densities,
    chunk_size=200_000
):
    """
    Calculate total tetrahedral volume and mass.

    ELE format:

        element_id node1 node2 node3 node4 material_id

    Tetrahedron volume:

        V = |(b-a) . ((c-a) x (d-a))| / 6

    Density is taken from MAT[material_id].

    Volume is assumed to be in cm^3 because NODE coordinates
    were converted to cm.
    """

    print(f"Reading ELE file: {ele_file}")
    print("Calculating tetrahedral volumes and mass...")

    total_volume = 0.0
    total_mass = 0.0

    n_elements = 0
    n_total_elements = None
    missing_materials = set()

    lines = []

    with open(ele_file, "r") as f:

        # ----------------------------------------------------
        # Read header
        # ----------------------------------------------------

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            header = line.split()

            n_total_elements = int(header[0])
            nodes_per_element = int(header[1])

            if nodes_per_element != 4:
                raise ValueError(
                    f"Expected tetrahedral elements with 4 nodes, "
                    f"but found {nodes_per_element}."
                )

            break

        if n_total_elements is None:
            raise ValueError("ELE file does not contain a valid header.")

        # ----------------------------------------------------
        # Process in chunks
        # ----------------------------------------------------

        def process_chunk(lines):

            nonlocal total_volume
            nonlocal total_mass
            nonlocal n_elements

            if not lines:
                return

            data = np.fromstring(
                "\n".join(lines),
                sep=" ",
                dtype=np.float64
            )

            data = data.reshape(-1, 6)

            # ------------------------------------------------
            # Columns:
            #
            # 0 = element ID
            # 1 = node 1
            # 2 = node 2
            # 3 = node 3
            # 4 = node 4
            # 5 = material ID
            # ------------------------------------------------

            node_ids = data[:, 1:5].astype(np.int64)
            material_ids = data[:, 5].astype(np.int64)

            # Get coordinates
            a = coordinates[node_ids[:, 0]]
            b = coordinates[node_ids[:, 1]]
            c = coordinates[node_ids[:, 2]]
            d = coordinates[node_ids[:, 3]]

            # Vectors
            ab = b - a
            ac = c - a
            ad = d - a

            # Tetrahedron volume
            volumes = np.abs(
                np.einsum(
                    "ij,ij->i",
                    ab,
                    np.cross(ac, ad)
                )
            ) / 6.0

            total_volume += np.sum(volumes)

            # ------------------------------------------------
            # Calculate mass
            # ------------------------------------------------

            chunk_densities = np.empty(
                len(material_ids),
                dtype=np.float64
            )

            for i, material_id in enumerate(material_ids):

                if material_id not in densities:

                    missing_materials.add(int(material_id))

                    # Don't silently use zero density
                    chunk_densities[i] = np.nan

                else:
                    chunk_densities[i] = densities[material_id]

            if np.any(np.isnan(chunk_densities)):

                raise ValueError(
                    "The following material IDs from the ELE file "
                    "were not found in the MATERIAL file:\n"
                    + ", ".join(
                        str(x)
                        for x in sorted(missing_materials)
                    )
                )

            total_mass += np.sum(
                volumes * chunk_densities
            )

            n_elements += len(data)

        # ----------------------------------------------------
        # Main loop
        # ----------------------------------------------------

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            lines.append(line)

            if len(lines) >= chunk_size:

                process_chunk(lines)

                lines.clear()

                if n_elements % 1_000_000 < chunk_size:
                    print(
                        f"  Processed approximately "
                        f"{n_elements:,} / "
                        f"{n_total_elements:,} elements..."
                    )

        # Process remaining elements
        process_chunk(lines)

    print(f"  Elements processed: {n_elements:,}")

    if n_elements != n_total_elements:
        print(
            f"WARNING: Header says {n_total_elements:,} elements, "
            f"but {n_elements:,} were read."
        )

    return total_volume, total_mass


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate phantom height, volume, and mass/weight "
            "from NODE, ELE, and MATERIAL files."
        )
    )

    parser.add_argument(
        "--node",
        required=True,
        help="Path to .node file"
    )

    parser.add_argument(
        "--ele",
        required=True,
        help="Path to .ele file"
    )

    parser.add_argument(
        "--material",
        required=True,
        help="Path to .material file"
    )

    parser.add_argument(
        "--coordinate-unit",
        choices=["mm", "cm", "m"],
        default="cm",
        help="Coordinate unit in NODE file (default: cm)"
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=200_000,
        help="Number of ELE rows processed at once"
    )

    args = parser.parse_args()

    print()
    print("=" * 60)
    print("PHANTOM HEIGHT AND WEIGHT CALCULATOR")
    print("=" * 60)
    print()

    # --------------------------------------------------------
    # Read nodes
    # --------------------------------------------------------

    coordinates = read_node_file(
        args.node,
        coordinate_unit=args.coordinate_unit
    )

    # --------------------------------------------------------
    # Calculate dimensions
    # --------------------------------------------------------

    dimensions = calculate_height(coordinates)

    # --------------------------------------------------------
    # Read material densities
    # --------------------------------------------------------

    densities = read_material_densities(
        args.material
    )

    # --------------------------------------------------------
    # Calculate volume and mass
    # --------------------------------------------------------

    volume_cm3, mass_g = calculate_volume_and_mass(
        args.ele,
        coordinates,
        densities,
        chunk_size=args.chunk_size
    )

    mass_kg = mass_g / 1000.0
    volume_liters = volume_cm3 / 1000.0

    # Weight force
    weight_newtons = mass_kg * 9.80665

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)

    print()
    print("DIMENSIONS")
    print("-" * 60)

    print(
        f"X range       : "
        f"{dimensions['x_min']:.3f} to "
        f"{dimensions['x_max']:.3f} cm"
    )

    print(
        f"Y range       : "
        f"{dimensions['y_min']:.3f} to "
        f"{dimensions['y_max']:.3f} cm"
    )

    print(
        f"Z range       : "
        f"{dimensions['z_min']:.3f} to "
        f"{dimensions['z_max']:.3f} cm"
    )

    print()
    print(f"Height        : {dimensions['height']:.3f} cm")
    print(f"Width         : {dimensions['width']:.3f} cm")
    print(f"Depth         : {dimensions['depth']:.3f} cm")

    print()
    print("MASS / WEIGHT")
    print("-" * 60)

    print(f"Total volume  : {volume_cm3:,.2f} cm³")
    print(f"Total volume  : {volume_liters:,.3f} L")

    print(f"Mass          : {mass_g:,.2f} g")
    print(f"Mass          : {mass_kg:,.3f} kg")

    print(
        f"Weight force  : {weight_newtons:,.2f} N"
    )

    print()
    print("=" * 60)
    print("Done.")
    print("=" * 60)


if __name__ == "__main__":
    main()