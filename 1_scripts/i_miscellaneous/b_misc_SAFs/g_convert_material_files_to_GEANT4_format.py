"""
Convert PHITS MAT[] .material files to the Geant4-style .material format.

Example input:

    $ Adrenal_left 1.02 g/cm3
    MAT[100]
         1000      -0.104
         6000      -0.222
         7000      -0.028
    mt100

becomes:

    C Adrenal_left 1.02 g/cm3
    m100     1000      -0.104
         6000      -0.222
         7000      -0.028
    C

The material composition data is preserved. Only the material-definition
format is changed.

Usage:
    python 1_scripts/i_miscellaneous/b_misc_SAFs/g_convert_material_files_to_GEANT4_format.py Filipino-MRCP-AF.material

Output:
    Filipino-MRCP-AF_Geant4.material

Or specify an output file:
    python 1_scripts/i_miscellaneous/b_misc_SAFs/g_convert_material_files_to_GEANT4_format.py Filipino-MRCP-AF.material -o Filipino-MRCP-AF_Geant4.material
"""

from pathlib import Path
import argparse
import re


# Matches:
#   $ Adrenal_left 1.02 g/cm3
PHITS_HEADER_RE = re.compile(
    r"^\s*\$\s+(.+?)\s+([+-]?\d+(?:\.\d+)?)\s+g/cm3\s*$",
    re.IGNORECASE,
)

# Matches:
#   MAT[100]
PHITS_MAT_RE = re.compile(
    r"^\s*MAT\[(\d+)\]\s*$",
    re.IGNORECASE,
)

# Matches:
#   mt100
PHITS_MT_RE = re.compile(
    r"^\s*mt\d+\s*$",
    re.IGNORECASE,
)

# Matches a material composition line:
#   1000      -0.104
#   6000      -0.222
COMPOSITION_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*$"
)


def convert_material_file(input_file, output_file):
    input_path = Path(input_file)
    output_path = Path(output_file)

    lines = input_path.read_text(encoding="utf-8").splitlines()

    output = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Look for the PHITS material header.
        header_match = PHITS_HEADER_RE.match(line)

        if not header_match:
            # Preserve anything that is not a material definition.
            output.append(line)
            i += 1
            continue

        material_name = header_match.group(1)
        density = header_match.group(2)

        # The next meaningful line should be MAT[ID].
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1

        if j >= len(lines):
            raise ValueError(
                f"Could not find MAT[] line after header: {line!r}"
            )

        mat_match = PHITS_MAT_RE.match(lines[j])

        if not mat_match:
            raise ValueError(
                f"Expected MAT[] after header {line!r}, "
                f"but found: {lines[j]!r}"
            )

        material_id = mat_match.group(1)

        # Collect the material composition until mt<ID>.
        composition = []
        k = j + 1

        while k < len(lines):
            current = lines[k]

            # End of this material.
            if PHITS_MT_RE.match(current):
                break

            # A new material header before mt<ID> means the input is malformed.
            if PHITS_HEADER_RE.match(current):
                raise ValueError(
                    f"Found a new material before mt{material_id}: "
                    f"{current!r}"
                )

            # Ignore blank lines inside a material block.
            if not current.strip():
                k += 1
                continue

            # Preserve only valid composition lines.
            composition_match = COMPOSITION_RE.match(current)

            if composition_match:
                element_id = composition_match.group(1)
                fraction = composition_match.group(2)
                composition.append((element_id, fraction))
            else:
                raise ValueError(
                    f"Unrecognized line inside material {material_id}: "
                    f"{current!r}"
                )

            k += 1

        if k >= len(lines) or not PHITS_MT_RE.match(lines[k]):
            raise ValueError(
                f"Could not find mt{material_id} for material "
                f"{material_id} ({material_name})."
            )

        if not composition:
            raise ValueError(
                f"Material {material_id} ({material_name}) has no "
                f"composition entries."
            )

        # Geant4-style material header.
        output.append(
            f"C {material_name} {density} g/cm3"
        )

        # Put the first composition entry on the m<ID> line,
        # matching the structure of M_H165W65(3).material.
        first_element, first_fraction = composition[0]

        output.append(
            f"m{material_id:<6}{first_element:<10}{first_fraction}"
        )

        # Remaining composition entries stay on separate lines.
        for element_id, fraction in composition[1:]:
            output.append(
                f"     {element_id:<10}{fraction}"
            )

        # Geant4 material blocks use C as the separator.
        output.append("C")

        # Move to the line after mt<ID>.
        i = k + 1

    output_path.write_text(
        "\n".join(output) + "\n",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert PHITS MAT[] .material files to "
            "Geant4-style material format."
        )
    )

    parser.add_argument(
        "input",
        help="Input PHITS .material file",
    )

    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output Geant4-style .material file. "
            "Defaults to <input>_Geant4.material"
        ),
    )

    args = parser.parse_args()

    input_path = Path(args.input)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input material file not found: {input_path}"
        )

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(
            f"{input_path.stem}_Geant4{input_path.suffix}"
        )

    convert_material_file(input_path, output_path)

    print(f"Converted: {input_path}")
    print(f"Output:    {output_path}")


if __name__ == "__main__":
    main()
