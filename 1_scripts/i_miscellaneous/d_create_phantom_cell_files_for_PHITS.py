"""
Generate a PHITS .cell file from a .vol and .material file.

The first five cells are fixed. The phantom name in cell 00001 is taken
from the .vol filename:

    F_H150W55.vol
        -> tfile=../phantoms/F_H150W55

The remaining cells are generated from the IDs/volumes in the .vol file
and densities in the .material file.

Usage:
    python 1_scripts/i_miscellaneous/d_create_phantom_cell_files_for_PHITS.py M_H165W65.vol M_H165W65.material

Optional:
    python 1_scripts/i_miscellaneous/d_create_phantom_cell_files_for_PHITS.py M_H165W65.vol M_H165W65.material -o output.cell
"""

import argparse
import re
from pathlib import Path


MATERIAL_HEADER_RE = re.compile(
    r"^\s*C\s+(.+?)\s+([0-9]+(?:\.[0-9]+)?)\s+g/cm3\s*$",
    re.IGNORECASE,
)

MATERIAL_ID_RE = re.compile(r"^\s*m(\d+)\s+", re.IGNORECASE)


def read_vol_file(vol_file):
    """Read material/organ IDs and volumes from the .vol file."""
    volumes = {}

    with open(vol_file, "r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()

            if not line or line.lower().startswith("id"):
                continue

            parts = line.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Invalid .vol line {line_number}: {line!r}"
                )

            try:
                material_id = int(parts[0])
                volume = float(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid .vol line {line_number}: {line!r}"
                ) from exc

            volumes[material_id] = volume

    return volumes


def read_material_file(material_file):
    """Read material IDs and their densities from the .material file."""
    densities = {}
    current_density = None

    with open(material_file, "r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):

            # Example:
            # C Adrenal_left 1.035 g/cm3
            header_match = MATERIAL_HEADER_RE.match(line)

            if header_match:
                current_density = float(header_match.group(2))
                continue

            # Example:
            # m100
            material_match = MATERIAL_ID_RE.match(line)

            if material_match:
                material_id = int(material_match.group(1))

                if current_density is None:
                    raise ValueError(
                        f"No density found for m{material_id} "
                        f"before line {line_number}."
                    )

                densities[material_id] = current_density

    return densities


def make_fixed_header(phantom_name):
    """Create the fixed first five lines of the .cell file."""
    return [
        "$ CELLS FOR PHANTOM",
        (
            f" 00001         0        -20"
            f"            U=15000 LAT=3 tfile=../phantoms/{phantom_name}"
        ),
        " 00002         0        -10            FILL=15000",
        " 00003         0        -20     10",
        " 00004         0        -90     20",
        " 00005        -1         90",
    ]


def make_cell_file(vol_file, material_file, output_file):
    volumes = read_vol_file(vol_file)
    densities = read_material_file(material_file)

    # The phantom name is the .vol filename without its extension.
    phantom_name = Path(vol_file).stem

    # Check that every volume ID has a corresponding material.
    missing_materials = sorted(set(volumes) - set(densities))

    if missing_materials:
        raise ValueError(
            "These IDs are present in the .vol file but missing from "
            f"the .material file: {missing_materials}"
        )

    # Create output.
    with open(output_file, "w", encoding="utf-8", newline="\n") as f:

        # ---------------------------------------------------------
        # Fixed phantom cells
        # ---------------------------------------------------------
        for line in make_fixed_header(phantom_name):
            f.write(line + "\n")

        # ---------------------------------------------------------
        # Organ/tissue cells
        # ---------------------------------------------------------
        for material_id in sorted(volumes):

            volume = volumes[material_id]
            density = densities[material_id]

            # Cell number is the material/organ ID padded to 5 digits.
            cell_number = f"{material_id:05d}"

            # PHITS uses negative mass density in the cell definition.
            f.write(
                f" {cell_number:<5}"
                f" {material_id:>10}"
                f" {-density:>11.3f}"
                f"    -90"
                f"     u={material_id:<8}"
                f"VOL={volume:.10f}"
                f"\n"
            )

    print(f"Created: {output_file}")
    print(f"Phantom: {phantom_name}")
    print(f"Generated {len(volumes)} organ/tissue cells.")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a PHITS .cell file from a .vol and .material file."
        )
    )

    parser.add_argument(
        "vol_file",
        type=Path,
        help="Input .vol file",
    )

    parser.add_argument(
        "material_file",
        type=Path,
        help="Input .material file",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .cell file (default: same name as .vol file)",
    )

    args = parser.parse_args()

    # Example:
    # F_H150W55.vol -> F_H150W55.cell
    output_file = (
        args.output
        if args.output is not None
        else args.vol_file.with_suffix(".cell")
    )

    make_cell_file(
        args.vol_file,
        args.material_file,
        output_file,
    )


if __name__ == "__main__":
    main()
