"""
Check ICRP 145 organ IDs against TWO phantom types simultaneously.

The ICRP 145 phantom can be male OR female, and the other phantom can
also be male OR female.

Each phantom is checked using its:
    - .cell file
    - .material file

The script produces:
    1. A full comparison table showing YES/NO for every ICRP 145 organ ID.
    2. A second output containing ONLY organs where at least one of the
       four checks is NO.

Usage
-----
python 1_scripts/i_miscellaneous/b_misc_SAFs/e_check_organ_IDs.py \
    organ_ID_names.csv \
    MRCP-AM.cell MRCP-AM.material \
    M_H165W65.cell M_H165W65.material

Example
-------
python 1_scripts/i_miscellaneous/b_misc_SAFs/e_check_organ_IDs.py \
    organ_ID_names.csv \
    MRCP-AM.cell MRCP-AM.material \
    MFCP-AM.cell MFCP-AM.material
"""

import argparse
import csv
import re
from pathlib import Path


CELL_LINE_RE = re.compile(r"^\s*\d{5}\s+(-?\d+)\s+")
MATERIAL_ID_RE = re.compile(r"^\s*m(\d+)\b", re.IGNORECASE)


def read_organ_csv(csv_file):
    """Read ICRP 145 organ IDs and names."""
    organs = {}

    with open(csv_file, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames:
            raise ValueError(f"CSV file is empty: {csv_file}")

        field_map = {
            field.strip().lower(): field
            for field in reader.fieldnames
        }

        id_field = None
        name_field = None

        for candidate in ("organ_id", "id", "organ id", "organid"):
            if candidate in field_map:
                id_field = field_map[candidate]
                break

        for candidate in ("name", "organ_name", "organ name", "organ"):
            if candidate in field_map:
                name_field = field_map[candidate]
                break

        if id_field is None or name_field is None:
            raise ValueError(
                "Could not find the organ ID/name columns in "
                f"{csv_file}. Found columns: {reader.fieldnames}"
            )

        for row in reader:
            if not str(row[id_field]).strip():
                continue

            organ_id = int(str(row[id_field]).strip())
            name = str(row[name_field]).strip()
            organs[organ_id] = name

    return organs


def read_cell_ids(cell_file):
    """Extract positive material/organ IDs from a PHITS .cell file."""
    cell_ids = set()

    with open(cell_file, "r", encoding="utf-8") as f:
        for line in f:
            match = CELL_LINE_RE.match(line)

            if match:
                material_id = int(match.group(1))

                if material_id > 0:
                    cell_ids.add(material_id)

    return cell_ids


def read_material_ids(material_file):
    """Extract m<ID> definitions from a PHITS .material file."""
    material_ids = set()

    with open(material_file, "r", encoding="utf-8") as f:
        for line in f:
            match = MATERIAL_ID_RE.match(line)

            if match:
                material_ids.add(int(match.group(1)))

    return material_ids


def check_phantom(cell_file, material_file, expected_ids):
    """Check cell and material IDs for one phantom."""
    cell_ids = read_cell_ids(cell_file)
    material_ids = read_material_ids(material_file)

    return {
        "cell": {
            organ_id: organ_id in cell_ids
            for organ_id in expected_ids
        },
        "material": {
            organ_id: organ_id in material_ids
            for organ_id in expected_ids
        },
    }


def yes_no(value):
    """Convert boolean to YES/NO."""
    return "YES" if value else "NO"


def write_full_csv(
    output_file,
    organs,
    phantom1,
    phantom2,
    result1,
    result2,
):
    """Write CSV containing every ICRP 145 organ ID."""
    fieldnames = [
        "organ_id",
        "organ_name",
        f"{phantom1['name']}_cell",
        f"{phantom1['name']}_material",
        f"{phantom2['name']}_cell",
        f"{phantom2['name']}_material",
    ]

    with open(
        output_file,
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for organ_id in sorted(organs):
            writer.writerow({
                "organ_id": organ_id,
                "organ_name": organs[organ_id],
                f"{phantom1['name']}_cell":
                    yes_no(result1["cell"][organ_id]),
                f"{phantom1['name']}_material":
                    yes_no(result1["material"][organ_id]),
                f"{phantom2['name']}_cell":
                    yes_no(result2["cell"][organ_id]),
                f"{phantom2['name']}_material":
                    yes_no(result2["material"][organ_id]),
            })


def write_missing_csv(
    output_file,
    organs,
    phantom1,
    phantom2,
    result1,
    result2,
):
    """
    Write ONLY organs for which at least one check is NO.
    """
    fieldnames = [
        "organ_id",
        "organ_name",
        f"{phantom1['name']}_cell",
        f"{phantom1['name']}_material",
        f"{phantom2['name']}_cell",
        f"{phantom2['name']}_material",
    ]

    with open(
        output_file,
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for organ_id in sorted(organs):

            checks = (
                result1["cell"][organ_id],
                result1["material"][organ_id],
                result2["cell"][organ_id],
                result2["material"][organ_id],
            )

            # Skip if everything is YES.
            if all(checks):
                continue

            writer.writerow({
                "organ_id": organ_id,
                "organ_name": organs[organ_id],
                f"{phantom1['name']}_cell":
                    yes_no(result1["cell"][organ_id]),
                f"{phantom1['name']}_material":
                    yes_no(result1["material"][organ_id]),
                f"{phantom2['name']}_cell":
                    yes_no(result2["cell"][organ_id]),
                f"{phantom2['name']}_material":
                    yes_no(result2["material"][organ_id]),
            })


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Check ICRP 145 organ IDs against two phantom "
            "cell/material pairs and save two CSV reports."
        )
    )

    parser.add_argument(
        "organ_csv",
        type=Path,
        help="CSV containing ICRP 145 organ IDs and names",
    )

    parser.add_argument(
        "phantom1_cell",
        type=Path,
        help="Phantom 1 .cell file",
    )

    parser.add_argument(
        "phantom1_material",
        type=Path,
        help="Phantom 1 .material file",
    )

    parser.add_argument(
        "phantom2_cell",
        type=Path,
        help="Phantom 2 .cell file",
    )

    parser.add_argument(
        "phantom2_material",
        type=Path,
        help="Phantom 2 .material file",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        default=None,
        help=(
            "Prefix for output CSV files. If omitted, the two phantom "
            "names are combined automatically."
        ),
    )

    args = parser.parse_args()

    # Read expected ICRP 145 organ IDs.
    organs = read_organ_csv(args.organ_csv)

    # Phantom names are taken from the .cell filenames.
    phantom1 = {
        "name": args.phantom1_cell.stem,
        "cell": args.phantom1_cell,
        "material": args.phantom1_material,
    }

    phantom2 = {
        "name": args.phantom2_cell.stem,
        "cell": args.phantom2_cell,
        "material": args.phantom2_material,
    }

    # Check both phantoms.
    result1 = check_phantom(
        args.phantom1_cell,
        args.phantom1_material,
        set(organs),
    )

    result2 = check_phantom(
        args.phantom2_cell,
        args.phantom2_material,
        set(organs),
    )

    # Output prefix.
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        prefix = f"{phantom1['name']}_{phantom2['name']}"

    full_output = Path(f"{prefix}_organ_id_check.csv")
    missing_output = Path(f"{prefix}_missing_organ_ids.csv")

    # Write both CSVs.
    write_full_csv(
        full_output,
        organs,
        phantom1,
        phantom2,
        result1,
        result2,
    )

    write_missing_csv(
        missing_output,
        organs,
        phantom1,
        phantom2,
        result1,
        result2,
    )

    # Summary for terminal.
    missing_count = sum(
        not all((
            result1["cell"][organ_id],
            result1["material"][organ_id],
            result2["cell"][organ_id],
            result2["material"][organ_id],
        ))
        for organ_id in organs
    )

    print()
    print("ICRP 145 ORGAN ID CHECK")
    print("=" * 60)
    print(f"Phantom 1: {phantom1['name']}")
    print(f"Phantom 2: {phantom2['name']}")
    print(f"Total ICRP 145 IDs: {len(organs)}")
    print(f"IDs with at least one NO: {missing_count}")
    print()
    print(f"Full output:    {full_output}")
    print(f"Missing output: {missing_output}")
    print()

    if missing_count == 0:
        print("All checks are YES.")
        print(
            "The missing-organ CSV was created with headers only."
        )
    else:
        print(
            f"{missing_count} organ(s) have at least one NO. "
            "See the missing-organ CSV."
        )


if __name__ == "__main__":
    main()
