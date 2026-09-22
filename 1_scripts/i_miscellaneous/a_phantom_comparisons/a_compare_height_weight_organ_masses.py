from pathlib import Path
import re
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.properties import PageSetupProperties


# ============================================================
# USER CONFIGURATION
# ============================================================

# Pick the TWO phantom groups you want to compare.
#
# Available in the configuration below:
#     "ICRP145"
#     "ICRP145_Resized"
#     "Filipino"
#
COMPARISON_GROUP_A = "ICRP145"
COMPARISON_GROUP_B = "Filipino"


# ------------------------------------------------------------
# Coordinate unit in NODE files
# ------------------------------------------------------------

NODE_UNIT = "cm"


# ------------------------------------------------------------
# Input CSV
# ------------------------------------------------------------

ORGAN_ID_NAMES_FILE = Path(
    r"/home/clarence/Geant4_SAF_Calculations/"
    r"PHITS_Geant4_python_scripts/5_other_input_files/"
    r"organ_ID_names.csv"
)


# ------------------------------------------------------------
# Phantom groups
# ------------------------------------------------------------

PHANTOM_GROUPS = {

    "ICRP145": {

        "Male": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AM.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AM.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AM.material"
            ),
        },

        "Female": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AF.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AF.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"MRCP-AF.material"
            ),
        },
    },


    "ICRP145_Resized": {

        "Male": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"M_H165W65.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"M_H165W65.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"M_H165W65.material"
            ),
        },

        "Female": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"F_H150W55.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"F_H150W55.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"F_H150W55.material"
            ),
        },
    },


    "Filipino": {

        "Male": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AM.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AM.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AM.material"
            ),
        },

        "Female": {
            "node": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AF.node"
            ),
            "ele": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AF.ele"
            ),
            "material": Path(
                r"/home/clarence/Geant4_SAF_Calculations/"
                r"PHITS_Geant4_python_scripts/2_phits/phantoms/"
                r"Filipino-MRCP-AF.material"
            ),
        },
    },
}


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUTPUT_EXCEL_FILE = Path(
    r"/home/clarence/Geant4_SAF_Calculations/"
    r"PHITS_Geant4_python_scripts/5_other_input_files/"
    r"phantom_height_weight_organ_comparison.xlsx"
)


# ============================================================
# DATA READING
# ============================================================

def load_organ_names(filename):
    df = pd.read_csv(filename)
    df.columns = [
        str(c).strip().lstrip("\ufeff")
        for c in df.columns
    ]

    if "organ_id" not in df.columns:
        raise ValueError(
            f"Could not find 'organ_id' in {filename}"
        )

    if "name" not in df.columns:
        raise ValueError(
            f"Could not find 'name' in {filename}"
        )

    return {
        int(row["organ_id"]): str(row["name"])
        for _, row in df.iterrows()
    }


def read_node_file(filename):
    print(f"Reading NODE: {filename}")

    with open(filename, "r") as f:
        dimension = None
        n_nodes = None

        nodes_per_element = None

        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            n_nodes = int(parts[0])
            dimension = int(parts[1])
            break

        if dimension is None or n_nodes is None:
            raise ValueError(
                f"Could not find a NODE header in {filename}."
            )

        if dimension != 3:
            raise ValueError(
                f"{filename} is not a 3D NODE file."
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

    if NODE_UNIT == "mm":
        coordinates /= 10.0
    elif NODE_UNIT == "m":
        coordinates *= 100.0
    elif NODE_UNIT != "cm":
        raise ValueError(
            "NODE_UNIT must be 'mm', 'cm', or 'm'."
        )

    return coordinates


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
# CALCULATIONS
# ============================================================

def calculate_dimensions(coordinates):
    return {
        "height_cm": (
            np.max(coordinates[:, 2])
            - np.min(coordinates[:, 2])
        ),
        "width_cm": (
            np.max(coordinates[:, 0])
            - np.min(coordinates[:, 0])
        ),
        "depth_cm": (
            np.max(coordinates[:, 1])
            - np.min(coordinates[:, 1])
        ),
    }


def calculate_organ_masses(
    ele_file,
    coordinates,
    densities,
    chunk_size=200_000,
):
    print(f"Reading ELE: {ele_file}")

    organ_masses = {}

    with open(ele_file, "r") as f:
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

        if n_elements is None or nodes_per_element is None:
            raise ValueError("The ELE file does not contain an element header.")

        if nodes_per_element != 4:
            raise ValueError(
                "This script expects tetrahedral elements."
            )

        chunk = []

        def process_chunk(lines):
            if not lines:
                return

            data = np.fromstring(
                "\n".join(lines),
                sep=" ",
                dtype=np.float64,
            ).reshape(-1, 6)

            node_ids = data[:, 1:5].astype(np.int64)
            material_ids = data[:, 5].astype(np.int64)

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

            for material_id in np.unique(material_ids):

                material_id = int(material_id)
                mask = material_ids == material_id

                if material_id not in densities:
                    raise ValueError(
                        f"Material ID {material_id} "
                        f"from {ele_file} was not found "
                        f"in the MATERIAL file."
                    )

                mass = np.sum(
                    volumes[mask]
                ) * densities[material_id]

                organ_masses[material_id] = (
                    organ_masses.get(material_id, 0.0)
                    + mass
                )

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

        process_chunk(chunk)

    return organ_masses


def process_phantom(group_name, sex, config):

    print()
    print("=" * 70)
    print(f"{group_name} - {sex}")
    print("=" * 70)

    coordinates = read_node_file(config["node"])
    dimensions = calculate_dimensions(coordinates)

    densities = read_material_densities(
        config["material"]
    )

    organ_masses = calculate_organ_masses(
        config["ele"],
        coordinates,
        densities,
    )

    total_mass_g = sum(
        organ_masses.values()
    )

    return {
        "group": group_name,
        "sex": sex,
        "height_cm": dimensions["height_cm"],
        "width_cm": dimensions["width_cm"],
        "depth_cm": dimensions["depth_cm"],
        "total_mass_g": total_mass_g,
        "organ_masses": organ_masses,
        "organ_densities": densities,
    }


def percent_difference(group_a_value, group_b_value):
    """
    Percentage difference relative to Group A:

        (Group B - Group A) / Group A * 100
    """
    if group_a_value == 0:
        return np.nan

    return (
        (group_b_value - group_a_value)
        / group_a_value
        * 100.0
    )


# ============================================================
# EXCEL FORMATTING
# ============================================================

THIN = Side(
    style="thin",
    color="000000"
)

MEDIUM = Side(
    style="medium",
    color="000000"
)

BORDER = Border(
    left=THIN,
    right=THIN,
    top=THIN,
    bottom=THIN,
)

MEDIUM_BORDER = Border(
    left=MEDIUM,
    right=MEDIUM,
    top=MEDIUM,
    bottom=MEDIUM,
)

BLUE_FILL = PatternFill(
    "solid",
    fgColor="D9EAF7"
)

PINK_FILL = PatternFill(
    "solid",
    fgColor="F4DFE5"
)

GRAY_FILL = PatternFill(
    "solid",
    fgColor="E7E7E7"
)

WHITE_FILL = PatternFill(
    "solid",
    fgColor="FFFFFF"
)


def set_cell(
    ws,
    row,
    column,
    value,
    bold=False,
    fill=None,
    align="center",
    border=BORDER,
    number_format=None,
):
    cell = ws.cell(
        row=row,
        column=column,
        value=value,
    )

    cell.font = Font(
        name="Arial",
        size=10,
        bold=bold,
    )

    cell.alignment = Alignment(
        horizontal=align,
        vertical="center",
        wrap_text=True,
    )

    cell.border = border

    if fill is not None:
        cell.fill = fill

    if number_format is not None:
        cell.number_format = number_format

    return cell


def write_sex_table(
    ws,
    start_row,
    start_col,
    sex,
    group_a,
    group_b,
    result_a,
    result_b,
    organ_names,
):
    """
    Creates one table like the figure:

    Male
    ------------------------------------------------
    Organ       Group A     Group B     % Difference
    Height      ...
    Weight      ...
    ------------------------------------------------
    Organ mass (g)
    Brain       ...
    ...
    """

    group_a_display = group_a.replace("_", " ")
    group_b_display = group_b.replace("_", " ")

    # --------------------------------------------------------
    # Title
    # --------------------------------------------------------

    ws.merge_cells(
        start_row=start_row,
        start_column=start_col,
        end_row=start_row,
        end_column=start_col + 5,
    )

    title = ws.cell(
        start_row,
        start_col,
        sex
    )

    title.font = Font(
        name="Arial",
        size=11,
        bold=True,
    )

    title.alignment = Alignment(
        horizontal="center",
        vertical="center",
    )

    title_fill = (
        BLUE_FILL
        if sex == "Male"
        else PINK_FILL
    )

    title.fill = title_fill

    for col in range(start_col, start_col + 6):
        ws.cell(
            start_row,
            col
        ).border = MEDIUM_BORDER
        ws.cell(
            start_row,
            col
        ).fill = title_fill

    # --------------------------------------------------------
    # Group names
    # --------------------------------------------------------

    header_row = start_row + 1

    headers = [
        "",
        f"{group_a_display} density",
        f"{group_a_display} mass",
        f"{group_b_display} density",
        f"{group_b_display} mass",
        "% Difference",
    ]

    for i, value in enumerate(headers):
        set_cell(
            ws,
            header_row,
            start_col + i,
            value,
            bold=True,
            fill=GRAY_FILL,
        )

    # --------------------------------------------------------
    # Height
    # --------------------------------------------------------

    row = start_row + 2

    set_cell(
        ws,
        row,
        start_col,
        "Height (cm)",
        bold=True,
        align="left",
    )

    set_cell(ws, row, start_col + 1, "")
    set_cell(
        ws,
        row,
        start_col + 2,
        result_a["height_cm"],
        number_format="0.00",
    )
    set_cell(ws, row, start_col + 3, "")
    set_cell(
        ws,
        row,
        start_col + 4,
        result_b["height_cm"],
        number_format="0.00",
    )
    set_cell(
        ws,
        row,
        start_col + 5,
        percent_difference(
            result_a["height_cm"],
            result_b["height_cm"],
        ),
        number_format="0.00",
    )

    # --------------------------------------------------------
    # Weight
    # --------------------------------------------------------

    row += 1

    set_cell(
        ws,
        row,
        start_col,
        "Weight (kg)",
        bold=True,
        align="left",
    )

    set_cell(ws, row, start_col + 1, "")
    set_cell(
        ws,
        row,
        start_col + 2,
        result_a["total_mass_g"] / 1000.0,
        number_format="0.00",
    )
    set_cell(ws, row, start_col + 3, "")
    set_cell(
        ws,
        row,
        start_col + 4,
        result_b["total_mass_g"] / 1000.0,
        number_format="0.00",
    )
    set_cell(
        ws,
        row,
        start_col + 5,
        percent_difference(
            result_a["total_mass_g"],
            result_b["total_mass_g"],
        ),
        number_format="0.00",
    )

    # --------------------------------------------------------
    # Organ mass header
    # --------------------------------------------------------

    row += 2

    ws.merge_cells(
        start_row=row,
        start_column=start_col,
        end_row=row,
        end_column=start_col + 5,
    )

    cell = ws.cell(
        row,
        start_col,
        "Density (g cm⁻³) / Organ mass (g)"
    )

    cell.font = Font(
        name="Arial",
        size=10,
        bold=True,
    )

    cell.alignment = Alignment(
        horizontal="center",
        vertical="center",
    )

    cell.fill = GRAY_FILL

    for col in range(start_col, start_col + 6):
        ws.cell(row, col).border = BORDER
        ws.cell(row, col).fill = GRAY_FILL

    # --------------------------------------------------------
    # Organ masses
    # --------------------------------------------------------

    all_ids = sorted(
        set(result_a["organ_masses"])
        | set(result_b["organ_masses"])
    )

    row += 1

    for organ_id in all_ids:

        mass_a = result_a["organ_masses"].get(
            organ_id,
            0.0
        )

        mass_b = result_b["organ_masses"].get(
            organ_id,
            0.0
        )

        density_a = result_a["organ_densities"].get(
            organ_id,
            np.nan
        )

        density_b = result_b["organ_densities"].get(
            organ_id,
            np.nan
        )

        difference = percent_difference(
            mass_a,
            mass_b,
        )

        organ_name = organ_names.get(
            organ_id,
            f"Unknown ({organ_id})"
        )

        set_cell(
            ws,
            row,
            start_col,
            organ_name,
            bold=True,
            align="left",
        )

        # Density immediately to the left of Group A mass.
        set_cell(
            ws,
            row,
            start_col + 1,
            density_a,
            number_format="0.000",
        )

        set_cell(
            ws,
            row,
            start_col + 2,
            mass_a,
            number_format="0.00",
        )

        # Density immediately to the left of Group B mass.
        set_cell(
            ws,
            row,
            start_col + 3,
            density_b,
            number_format="0.000",
        )

        set_cell(
            ws,
            row,
            start_col + 4,
            mass_b,
            number_format="0.00",
        )

        set_cell(
            ws,
            row,
            start_col + 5,
            difference,
            number_format="0.00",
        )

        row += 1

    return row


# ============================================================
# CREATE EXCEL
# ============================================================

def create_excel(
    results,
    organ_names,
    output_file,
):

    result_lookup = {
        (r["group"], r["sex"]): r
        for r in results
    }

    wb = Workbook()
    ws = wb.active
    assert ws is not None
    ws.title = "Comparison"

    group_a_display = (
        COMPARISON_GROUP_A.replace("_", " ")
    )

    group_b_display = (
        COMPARISON_GROUP_B.replace("_", " ")
    )

    # --------------------------------------------------------
    # Main title
    # --------------------------------------------------------

    ws.merge_cells(
        "A1:M1"
    )

    title = ws["A1"]
    title.value = (
        "Comparison of Height, Weight, and Organ Masses"
    )
    title.font = Font(
        name="Arial",
        size=14,
        bold=True,
    )
    title.alignment = Alignment(
        horizontal="center"
    )

    # --------------------------------------------------------
    # Definition of percentage difference
    # --------------------------------------------------------

    ws.merge_cells("A2:M2")

    ws["A2"] = (
        f"% Difference = "
        f"({group_b_display} − {group_a_display}) "
        f"/ {group_a_display} × 100"
    )

    ws["A2"].font = Font(
        name="Arial",
        size=9,
        italic=True,
    )

    ws["A2"].alignment = Alignment(
        horizontal="center"
    )

    # --------------------------------------------------------
    # Male table
    # --------------------------------------------------------

    male_a = result_lookup[
        (COMPARISON_GROUP_A, "Male")
    ]

    male_b = result_lookup[
        (COMPARISON_GROUP_B, "Male")
    ]

    write_sex_table(
        ws,
        start_row=4,
        start_col=1,
        sex="Male",
        group_a=COMPARISON_GROUP_A,
        group_b=COMPARISON_GROUP_B,
        result_a=male_a,
        result_b=male_b,
        organ_names=organ_names,
    )

    # --------------------------------------------------------
    # Female table
    # --------------------------------------------------------

    female_a = result_lookup[
        (COMPARISON_GROUP_A, "Female")
    ]

    female_b = result_lookup[
        (COMPARISON_GROUP_B, "Female")
    ]

    write_sex_table(
        ws,
        start_row=4,
        start_col=8,
        sex="Female",
        group_a=COMPARISON_GROUP_A,
        group_b=COMPARISON_GROUP_B,
        result_a=female_a,
        result_b=female_b,
        organ_names=organ_names,
    )

    # --------------------------------------------------------
    # Column widths
    # --------------------------------------------------------

    widths = {
        "A": 27,
        "B": 13,
        "C": 13,
        "D": 13,
        "E": 13,
        "F": 13,
        "G": 3,
        "H": 27,
        "I": 13,
        "J": 13,
        "K": 13,
        "L": 13,
        "M": 13,
    }

    for column, width in widths.items():
        ws.column_dimensions[column].width = width

    # --------------------------------------------------------
    # Freeze panes
    # --------------------------------------------------------

    ws.freeze_panes = "A6"

    # --------------------------------------------------------
    # Page setup
    # --------------------------------------------------------

    ws.sheet_view.showGridLines = False

    ws.page_setup.orientation = "landscape"
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0

    ws.sheet_properties.pageSetUpPr = PageSetupProperties(fitToPage=True)

    ws.print_title_rows = "1:5"

    # --------------------------------------------------------
    # Workbook metadata
    # --------------------------------------------------------

    wb.properties.title = (
        "Phantom Height Weight Organ Mass Comparison"
    )

    wb.properties.subject = (
        f"{group_a_display} vs {group_b_display}"
    )

    wb.save(output_file)


# ============================================================
# MAIN
# ============================================================

def main():

    if COMPARISON_GROUP_A not in PHANTOM_GROUPS:
        raise ValueError(
            f"Unknown phantom group: "
            f"{COMPARISON_GROUP_A}"
        )

    if COMPARISON_GROUP_B not in PHANTOM_GROUPS:
        raise ValueError(
            f"Unknown phantom group: "
            f"{COMPARISON_GROUP_B}"
        )

    if (
        COMPARISON_GROUP_A
        == COMPARISON_GROUP_B
    ):
        raise ValueError(
            "COMPARISON_GROUP_A and "
            "COMPARISON_GROUP_B must be different."
        )

    print()
    print("=" * 80)
    print("PHANTOM COMPARISON")
    print("=" * 80)

    print(
        f"\nComparing:"
        f"\n  A = {COMPARISON_GROUP_A}"
        f"\n  B = {COMPARISON_GROUP_B}"
    )

    organ_names = load_organ_names(
        ORGAN_ID_NAMES_FILE
    )

    results = []

    for group_name in (
        COMPARISON_GROUP_A,
        COMPARISON_GROUP_B,
    ):

        for sex in ("Male", "Female"):

            config = PHANTOM_GROUPS[
                group_name
            ][sex]

            result = process_phantom(
                group_name,
                sex,
                config,
            )

            results.append(result)

    create_excel(
        results,
        organ_names,
        OUTPUT_EXCEL_FILE,
    )

    print()
    print("=" * 80)
    print("DONE")
    print("=" * 80)

    print(
        f"\nExcel file:"
        f"\n{OUTPUT_EXCEL_FILE}"
    )

    print(
        "\n% Difference is calculated relative to "
        f"{COMPARISON_GROUP_A}."
    )


if __name__ == "__main__":
    main()
