"""
Convert Geant4 .material files to the PHITS MAT[] format.

    C Adrenal_left 1.035 g/cm3
    m100
         1000      -0.104
         6000      -0.233
         7000      -0.028
    C
    mt100

or, when the first composition entry is on the m-line:

    C Adrenal_left 1.035 g/cm3
    m100     1000      -0.104
         6000      -0.233
         7000      -0.028
    C
    mt100

Both become:

    $ Adrenal_left 1.035 g/cm3
    MAT[100]
         1000      -0.104
         6000      -0.233
         7000      -0.028
    mt100

The material composition data is preserved. Only the material-definition
format is changed.

Usage:
    python 1_scripts/i_miscellaneous/b_misc_SAFs/f_convert_material_files_to_PHITS_format.py M_H165W65.material

Output:
    input_converted.material

Or specify an output file:
    python 1_scripts/i_miscellaneous/b_misc_SAFs/f_convert_material_files_to_PHITS_format.py input.material -o output.material
"""

from pathlib import Path
import argparse
import re


MATERIAL_HEADER_RE = re.compile(
    r"^\s*C\s+(.+?)\s+([+-]?\d+(?:\.\d+)?)\s+g/cm3\s*$",
    re.IGNORECASE,
)

# Matches the old material ID and optionally captures a composition entry
# appearing on the same line:
#     m100
#     m100  1000  -0.104
MATERIAL_ID_RE = re.compile(
    r"^\s*m(\d+)"
    r"(?:\s+(\d+)\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?))?"
    r"\s*$",
    re.IGNORECASE,
)

MT_RE = re.compile(r"^\s*mt\d+\s*$", re.IGNORECASE)
SEPARATOR_RE = re.compile(r"^\s*C\s*$", re.IGNORECASE)


def convert_material_file(input_file, output_file):
    input_path = Path(input_file)
    output_path = Path(output_file)

    lines = input_path.read_text(encoding="utf-8").splitlines()

    output = []
    i = 0

    while i < len(lines):
        line = lines[i]
        header_match = MATERIAL_HEADER_RE.match(line)

        if not header_match:
            output.append(line)
            i += 1
            continue

        material_name = header_match.group(1)
        density = header_match.group(2)

        # Find the corresponding m<ID> line.
        j = i + 1
        material_id = None
        first_composition = None

        while j < len(lines):
            id_match = MATERIAL_ID_RE.match(lines[j])

            if id_match:
                material_id = id_match.group(1)

                # If composition was written on the m<ID> line,
                # move it below MAT[ID].
                if id_match.group(2) is not None:
                    first_composition = (
                        id_match.group(2),
                        id_match.group(3),
                    )
                break

            if MATERIAL_HEADER_RE.match(lines[j]):
                break

            j += 1

        if material_id is None:
            raise ValueError(
                f"Could not find material ID after header: {line!r}"
            )

        # Converted header and material ID.
        output.append(f"$ {material_name} {density} g/cm3")
        output.append(f"MAT[{material_id}]")

        # Move any composition entry that was on m<ID> below MAT[ID].
        if first_composition is not None:
            element_id, fraction = first_composition
            output.append(f"     {element_id:<10}{fraction}")

        # Read the remaining material composition lines.
        k = j + 1

        while k < len(lines):
            current = lines[k]

            # Next material starts here.
            if MATERIAL_HEADER_RE.match(current):
                break

            # Old m<ID> marker: already converted.
            if MATERIAL_ID_RE.match(current):
                k += 1
                continue

            # Standalone C was only an old-format separator.
            if SEPARATOR_RE.match(current):
                k += 1
                continue

            # Old mt<ID> marker is rebuilt below.
            if MT_RE.match(current):
                k += 1
                continue

            # Preserve every composition entry exactly as supplied.
            output.append(current)
            k += 1

        # Rebuild mt<ID> after all composition entries.
        output.append(f"mt{material_id}")

        i = k

    output_path.write_text("\n".join(output) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Convert old-format PHITS .material files to MAT[] format."
    )
    parser.add_argument("input", help="Input .material file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output .material file. Defaults to <input>_converted.material",
    )

    args = parser.parse_args()

    input_path = Path(args.input)

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(
            f"{input_path.stem}_converted{input_path.suffix}"
        )

    convert_material_file(input_path, output_path)

    print(f"Converted: {input_path}")
    print(f"Output:    {output_path}")


if __name__ == "__main__":
    main()
