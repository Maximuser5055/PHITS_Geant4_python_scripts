"""
Compare chord-length distributions (CLDs) for mesh-type computational phantoms.

Default comparison:
    - Filipino Adult MRCP
    - ICRP 145 Adult MRCP resized to Filipino height/weight

The script uses:
    SOURCE_CSV
    target_region_file
    NODE/ELE/MATERIAL files from the phantom registry

For every source-region/target-region pair:
    1. Tetrahedra belonging to each region are collected.
    2. If a region contains multiple organ IDs, each ID is selected
       with probability proportional to its mass.
    3. A tetrahedron is selected proportional to its volume within
       that organ.
    4. A point is sampled uniformly inside the tetrahedron.
    5. Straight-line source-to-target distances are accumulated into
       1-mm bins.

The output is:
    - one CLD CSV for every phantom/source/target combination
    - one summary CSV containing mean and standard deviation
    - one user-selected multi-panel vector PDF, with all selected phantom
      groups combined in every plot

Existing CLD CSVs are reused. If all expected results already exist, the
user can skip recalculation or redo all calculations; if only some exist,
only the missing CLDs are calculated.

The NODE/ELE format follows the ICRP mesh phantom convention:
NODE files contain node coordinates and ELE files contain tetrahedra
plus their organ-ID attribute.
"""

from __future__ import annotations

from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import hashlib
import math
import os
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd

# Figures are written as vector PDF, so use a non-interactive backend
# suitable for headless systems.
# This also avoids Qt-specific rendering problems on headless systems.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Preferred journal font. Arial is used whenever it is available.
JOURNAL_FONT = "Arial"

# Linux/WSL fallbacks used only when Arial is unavailable. For final journal
# submission, install/use Arial and regenerate the figures.
FONT_FALLBACKS = (
    "Liberation Sans",
    "DejaVu Sans",
)

from matplotlib import font_manager


def configure_journal_font() -> str:
    """Select an installed sans-serif font, preferring Arial."""
    selected_font = None

    # Arial already registered with Matplotlib.
    try:
        font_manager.findfont(
            JOURNAL_FONT,
            fallback_to_default=False,
        )
        selected_font = JOURNAL_FONT
    except ValueError:
        pass

    # In WSL, Windows fonts are commonly accessible here.
    if selected_font is None:
        windows_font_dir = Path("/mnt/c/Windows/Fonts")
        arial_candidates = [
            windows_font_dir / "arial.ttf",
            windows_font_dir / "Arial.ttf",
            windows_font_dir / "ARIAL.TTF",
        ]

        for font_file in arial_candidates:
            if font_file.is_file():
                font_manager.fontManager.addfont(str(font_file))
                selected_font = JOURNAL_FONT
                break

    # Portable Linux fallbacks.
    if selected_font is None:
        for fallback in FONT_FALLBACKS:
            try:
                font_manager.findfont(
                    fallback,
                    fallback_to_default=False,
                )
                selected_font = fallback
                break
            except ValueError:
                continue

    if selected_font is None:
        raise RuntimeError(
            "No suitable journal sans-serif font was found. "
            "Install Arial, Liberation Sans, or DejaVu Sans."
        )

    if selected_font != JOURNAL_FONT:
        print(
            f"WARNING: {JOURNAL_FONT!r} was not found. "
            f"Using {selected_font!r} for this run. "
            f"For final journal submission, install Arial and "
            f"regenerate the PDF/TIFF figures."
        )

    plt.rcParams.update({
        "font.family": selected_font,
        "font.size": 8,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 5.5,

        # Embed TrueType fonts in vector PDF output.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    return selected_font

import b_config.a_config as config
from b_config import b_phantom_registry as registry
from c_database import b_organ_database as organ_database

# ======================================================================
# USER SETTINGS
# ======================================================================

# Registry group codes to compare.
# Phantom file paths are NOT hard-coded here.
PHANTOM_GROUPS_TO_COMPARE = [
    "MRCP_AF_AM",
    "MRCP_AF_AM_Filipino_Resized",
]

# Output directory name under config.RESULTS_DIR.
OUTPUT_DIR_NAME = "chord_length_distributions"

# Number of source-target point pairs per CLD.
N_SAMPLES = 10_000_000

# 1-mm bins.
BIN_WIDTH_MM = 1.0

# Number of samples processed by each worker at one time.
CHUNK_SIZE = 250_000

# Reproducible random sampling. Set to None for a random run.
RANDOM_SEED = 20261004

# Units of coordinates in the TETGEN NODE files.
NODE_UNIT = "cm"

# Plot settings.
# FIGURE_WIDTH_MM is the width of one subplot/panel. The final figure
# size scales automatically with the user-selected rows x columns layout.
FIGURE_WIDTH_MM = 70.0
FIGURE_WIDTH = FIGURE_WIDTH_MM / 25.4

# Round each x-axis maximum upward to the nearest 50 mm.
# Examples: 1723 mm -> 1750 mm, 1895 mm -> 1900 mm.
PLOT_X_ROUNDING_MM = 50.0

# Calculate CLDs for all source organs. Plot selection is handled separately
# after the CSV results have been generated.
ASK_SOURCE_SELECTION = False


# ======================================================================
# DATA CLASSES
# ======================================================================

@dataclass
class MeshData:
    """Tetrahedral mesh data."""

    # Coordinates remain in cm to match the reference mass calculation.
    nodes_cm: np.ndarray
    tetrahedra: np.ndarray
    organ_ids: np.ndarray


@dataclass
class RegionSampler:
    """Sampling data for one individual or compound organ region."""

    region_name: str
    organ_ids: tuple[int, ...]
    tetrahedra: np.ndarray
    cumulative_weights: np.ndarray
    total_mass_g: float


@dataclass
class PhantomInfo:
    """Resolved phantom information."""

    code: str
    display_name: str
    sex: str
    node_file: Path
    element_file: Path
    material_file: Path


# ======================================================================
# GENERAL UTILITIES
# ======================================================================

def available_threads() -> int:
    """Return the number of logical CPUs available to Python."""
    return os.cpu_count() or 1


def ask_threads() -> int:
    """Ask for the number of worker threads and cap at available CPUs."""

    maximum = available_threads()

    print()
    print("=" * 70)
    print("MULTITHREADING")
    print("=" * 70)
    print(f"Available logical CPU threads: {maximum}")

    while True:
        value = input(
            f"Number of threads to use [1-{maximum}]: "
        ).strip()

        try:
            threads = int(value)
        except ValueError:
            print("Please enter an integer.")
            continue

        if threads < 1:
            print("Threads must be at least 1.")
            continue

        if threads > maximum:
            print(
                f"Requested {threads} threads, but only {maximum} "
                f"are available. Using {maximum}."
            )
            return maximum

        return threads


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Strip whitespace/BOM from CSV column names."""

    df = df.copy()
    df.columns = [
        str(column).strip().lstrip("\ufeff")
        for column in df.columns
    ]
    return df


# ======================================================================
# SOURCE / TARGET CSV
# ======================================================================
# Source regions use the same compound-region representation as targets:
# one region can contain one or more organ IDs, stored in `ids`.


def parse_region_ids(
    ids_text: str,
    region_name: str,
    row_number: int,
    source_or_target: str,
) -> tuple[int, ...]:
    """Parse underscore-separated organ IDs from a region definition."""
    ids = []

    for value in str(ids_text).split("_"):
        value = value.strip()

        if not value:
            continue

        try:
            ids.append(int(float(value)))
        except ValueError as exc:
            raise ValueError(
                f"Invalid organ ID '{value}' for '{region_name}' "
                f"at row {row_number} in the {source_or_target}-region file."
            ) from exc

    if not ids:
        raise ValueError(
            f"No valid organ IDs for '{region_name}' "
            f"at row {row_number} in the {source_or_target}-region file."
        )

    # Remove accidental duplicate IDs while preserving the CSV order.
    return tuple(dict.fromkeys(ids))


def load_source_regions(source_csv: Path) -> list[dict]:
    """
    Load source-region definitions.

    Preferred source CSV columns:
        Source region
        Acronym
        ID number(s)

    IDs are underscore-separated, e.g.:
        12000_12100
        10000_10100_10200_10300_10400_10500

    For backward compatibility, a CSV containing only:
        source_organ_ID
    is also accepted. Each ID becomes a one-ID source region whose key is
    the numeric ID as a string.
    """

    if not source_csv.is_file():
        raise FileNotFoundError(
            f"SOURCE_CSV was not found:\n{source_csv}"
        )

    df = clean_columns(pd.read_csv(source_csv))

    compound_columns = [
        "Source region",
        "Acronym",
        "ID number(s)",
    ]

    if all(column in df.columns for column in compound_columns):
        regions = []

        for row_number, (_, row) in enumerate(df.iterrows(), start=2):
            if pd.isna(row["Source region"]):
                continue

            name = str(row["Source region"]).strip()
            acronym = (
                str(row["Acronym"]).strip()
                if pd.notna(row["Acronym"])
                else ""
            )
            ids_text = (
                str(row["ID number(s)"]).strip()
                if pd.notna(row["ID number(s)"])
                else ""
            )

            if not name:
                raise ValueError(
                    f"Empty 'Source region' at row {row_number} "
                    f"in {source_csv}."
                )

            if not ids_text:
                raise ValueError(
                    f"Empty 'ID number(s)' for '{name}' "
                    f"at row {row_number} in {source_csv}."
                )

            ids = parse_region_ids(
                ids_text,
                name,
                row_number,
                "source",
            )

            # The acronym is the stable internal key when available;
            # otherwise the source-region name is used.
            key = acronym or name

            if any(existing["key"] == key for existing in regions):
                raise ValueError(
                    f"Duplicate source-region key '{key}' in {source_csv}. "
                    "Use unique values in the 'Acronym' column."
                )

            regions.append(
                {
                    "key": key,
                    "name": name,
                    "acronym": acronym,
                    "ids": ids,
                }
            )

        if not regions:
            raise ValueError(
                f"No source regions were found in {source_csv}."
            )

        return regions

    # Backward-compatible support for the previous one-ID format.
    if "source_organ_ID" in df.columns:
        regions = []
        seen_ids = set()

        for value in df["source_organ_ID"].dropna():
            try:
                source_key = int(float(str(value).strip()))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid source organ ID '{value}' in {source_csv}."
                ) from exc

            if source_key in seen_ids:
                continue

            seen_ids.add(source_key)
            regions.append(
                {
                    "key": str(source_key),
                    "name": "",
                    "acronym": str(source_key),
                    "ids": (source_key,),
                }
            )

        if not regions:
            raise ValueError(
                f"No source-organ IDs were found in {source_csv}."
            )

        print(
            "WARNING: SOURCE_CSV is using the legacy 'source_organ_ID' "
            "format. Compound source regions require the columns "
            "'Source region', 'Acronym', and 'ID number(s)'."
        )

        return regions

    raise ValueError(
        f"SOURCE_CSV must contain either the compound source-region columns "
        f"{compound_columns} or the legacy column 'source_organ_ID'.\n"
        f"Available columns: {list(df.columns)}"
    )


def load_target_regions(target_region_file: Path) -> list[dict]:
    """
    Load target-region definitions.

    Expected columns:
        Target region
        Acronym
        ID number(s)

    IDs are underscore-separated, e.g.:
        9700_9900
    """

    if not target_region_file.is_file():
        raise FileNotFoundError(
            f"Target-region file was not found:\n"
            f"{target_region_file}"
        )

    df = clean_columns(
        pd.read_csv(
            target_region_file,
            dtype={
                "Target region": "string",
                "Acronym": "string",
                "ID number(s)": "string",
            },
        )
    )

    required = [
        "Target region",
        "Acronym",
        "ID number(s)",
    ]

    missing = [
        column
        for column in required
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing target-region column(s): {missing}\n"
            f"Expected: {required}"
        )

    regions = []

    for row_number, (_, row) in enumerate(df.iterrows(), start=2):

        if pd.isna(row["Target region"]):
            continue

        name = str(row["Target region"]).strip()
        acronym = (
            str(row["Acronym"]).strip()
            if pd.notna(row["Acronym"])
            else ""
        )
        ids_text = (
            str(row["ID number(s)"]).strip()
            if pd.notna(row["ID number(s)"])
            else ""
        )

        if not ids_text:
            raise ValueError(
                f"Empty 'ID number(s)' for '{name}' "
                f"at row {row_number}."
            )

        ids = []

        for value in ids_text.split("_"):
            value = value.strip()

            if not value:
                continue

            try:
                ids.append(int(float(value)))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid organ ID '{value}' for '{name}' "
                    f"at row {row_number}."
                ) from exc

        if not ids:
            raise ValueError(
                f"No valid organ IDs for '{name}' "
                f"at row {row_number}."
            )

        regions.append(
            {
                "name": name,
                "acronym": acronym,
                "ids": tuple(ids),
            }
        )

    if not regions:
        raise ValueError(
            f"No target regions were found in {target_region_file}."
        )

    return regions


def merge_gonads_target_regions(target_regions: list[dict]) -> list[dict]:
    """
    Keep separate Testes and Ovaries target regions and also add a combined
    Gonads region.

    The source target-region files are allowed to use either:
        - Testes / LOvary / ROvary
        - Testes / Ovaries

    If only left/right ovary entries are present, they are combined into the
    canonical ``Ovaries`` target region.  The original Testes and Ovaries
    regions are retained, and ``Gonads`` is added as their union.
    """
    testes_ids = []
    ovary_ids = []
    merged = []
    has_ovaries_region = False

    for region in target_regions:
        name = str(region["name"]).strip()
        name_lower = name.lower()

        if "test" in name_lower:
            for organ_id in region["ids"]:
                if organ_id not in testes_ids:
                    testes_ids.append(organ_id)
            # Keep the separate Testes region.
            merged.append(region)

        elif "ovar" in name_lower:
            for organ_id in region["ids"]:
                if organ_id not in ovary_ids:
                    ovary_ids.append(organ_id)

            # Keep an existing compound Ovaries region, but do not also keep
            # separate LOvary/ROvary entries in the plot-selection list.
            if name_lower == "ovaries":
                has_ovaries_region = True
                merged.append(region)

        else:
            merged.append(region)

    # If the target file only contains LOvary/ROvary, create one canonical
    # Ovaries region so the plot selection is consistent across target files.
    if ovary_ids and not has_ovaries_region:
        merged.append({
            "name": "Ovaries",
            "acronym": "Ovaries",
            "ids": tuple(ovary_ids),
        })

    gonad_ids = []
    for organ_id in (*testes_ids, *ovary_ids):
        if organ_id not in gonad_ids:
            gonad_ids.append(organ_id)

    if gonad_ids:
        merged.append({
            "name": "Gonads",
            "acronym": "Gonads",
            "ids": tuple(gonad_ids),
        })

    return merged


# ======================================================================
# SEX-SPECIFIC TARGET REGIONS
# ======================================================================

# Organ IDs that do not exist in the opposite-sex phantom.
# These are skipped rather than treated as errors when building
# target-region samplers.
SEX_EXCLUDED_ORGAN_IDS = {
    "AF": {
        11500,       # Prostate
        12900, 13000 # Testes
    },
    "AM": {
        11100, 11200, # Ovaries
        13900         # Uterus/cervix
    },
}


def filter_target_regions_for_sex(
    target_regions: list[dict],
    sex: str,
) -> list[dict]:
    """
    Remove sex-specific organ IDs that are not present in this phantom.

    Female (AF):
        skip prostate and testes IDs.

    Male (AM):
        skip ovaries and uterus/cervix IDs.

    If removing the sex-specific IDs leaves a target region with no IDs,
    that target region is skipped for that phantom.
    """

    excluded_ids = SEX_EXCLUDED_ORGAN_IDS.get(sex, set())

    filtered_regions = []

    for region in target_regions:

        remaining_ids = tuple(
            organ_id
            for organ_id in region["ids"]
            if organ_id not in excluded_ids
        )

        if not remaining_ids:
            print(
                f"    [SKIP] {region['name']} "
                f"({sex} phantom): sex-specific organ IDs "
                f"{region['ids']}"
            )
            continue

        if len(remaining_ids) != len(region["ids"]):

            skipped_ids = tuple(
                organ_id
                for organ_id in region["ids"]
                if organ_id in excluded_ids
            )

            print(
                f"    [SKIP IDs] {region['name']} "
                f"({sex} phantom): {skipped_ids}"
            )

        filtered_region = dict(region)
        filtered_region["ids"] = remaining_ids
        filtered_regions.append(filtered_region)

    return filtered_regions


# ======================================================================
# NODE / ELE / MATERIAL FILE READING
# ======================================================================

def read_node_file(filename: Path) -> np.ndarray:
    """
    Read a TETGEN .node file using the same convention as the
    compare_height_weight_organ_masses reference.

    Coordinates are returned in the original NODE_UNIT, which is cm.
    """

    if not filename.is_file():
        raise FileNotFoundError(
            f"NODE file not found:\n{filename}"
        )

    with filename.open(
        "r",
        encoding="utf-8",
        errors="ignore",
    ) as file:

        dimension = None
        n_nodes = None

        for line in file:
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
            dtype=np.float64,
        )

        for line in file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()
            node_id = int(parts[0])

            if not (0 <= node_id < n_nodes):
                raise ValueError(
                    f"Node ID {node_id} in {filename} is outside "
                    f"the expected range 0-{n_nodes - 1}."
                )

            coordinates[node_id] = [
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            ]

    if NODE_UNIT == "cm":
        return coordinates

    raise ValueError("NODE_UNIT must be 'cm'.")


def read_ele_file(
    filename: Path,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Read a TETGEN .ele file.

    The fifth field after the element ID contains the organ/material ID,
    following the same convention used by the reference script.
    """

    if not filename.is_file():
        raise FileNotFoundError(
            f"ELE file not found:\n{filename}"
        )

    with filename.open(
        "r",
        encoding="utf-8",
        errors="ignore",
    ) as file:

        number_of_elements = None
        nodes_per_element = None

        for line in file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            header = line.split()

            number_of_elements = int(header[0])
            nodes_per_element = int(header[1])
            break

        if (
            number_of_elements is None
            or nodes_per_element is None
        ):
            raise ValueError(
                f"The ELE file does not contain an element header:\n"
                f"{filename}"
            )

        if nodes_per_element != 4:
            raise ValueError(
                f"{filename} does not contain tetrahedral elements."
            )

        tetrahedra = np.empty(
            (number_of_elements, 4),
            dtype=np.int64,
        )

        organ_ids = np.empty(
            number_of_elements,
            dtype=np.int64,
        )

        count = 0

        for line in file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 6:
                raise ValueError(
                    f"Invalid ELE line in {filename}:\n{line}"
                )

            tetrahedra[count] = [
                int(parts[1]),
                int(parts[2]),
                int(parts[3]),
                int(parts[4]),
            ]

            organ_ids[count] = int(parts[5])
            count += 1

            if count == number_of_elements:
                break

    if count != number_of_elements:
        raise ValueError(
            f"Expected {number_of_elements:,} elements in {filename}, "
            f"but read {count:,}."
        )

    # Match the reference file convention.
    minimum_node_id = tetrahedra.min()

    if minimum_node_id == 1:
        tetrahedra -= 1
    elif minimum_node_id != 0:
        tetrahedra -= minimum_node_id

    return tetrahedra, organ_ids


def read_material_densities(
    filename: Path,
) -> dict[int, float]:
    """
    Read organ/material densities using the material-file convention
    used by compare_height_weight_organ_masses.py:

        $ <material name> <density> g/cm3
        MAT[<material ID>]

    The density is associated with the following MAT[ID].
    """

    if not filename.is_file():
        raise FileNotFoundError(
            f"MATERIAL file not found:\n{filename}"
        )

    densities: dict[int, float] = {}
    current_density: float | None = None

    density_pattern = re.compile(
        r"^\s*\$\s+.*?([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+g/cm3",
        re.IGNORECASE,
    )

    material_pattern = re.compile(
        r"^\s*MAT\[(\d+)\]",
        re.IGNORECASE,
    )

    with filename.open(
        "r",
        encoding="utf-8",
        errors="ignore",
    ) as file:

        for line in file:
            line = line.strip()

            if not line:
                continue

            density_match = density_pattern.match(line)

            if density_match:
                current_density = float(
                    density_match.group(1)
                )
                continue

            material_match = material_pattern.match(line)

            if material_match:
                if current_density is None:
                    raise ValueError(
                        f"MAT found without a preceding density in "
                        f"{filename}: {line}"
                    )

                material_id = int(
                    material_match.group(1)
                )

                densities[material_id] = current_density
                current_density = None

    if not densities:
        raise ValueError(
            f"No organ densities could be read from:\n{filename}"
        )

    return densities


# ======================================================================
# MESH PREPARATION
# ======================================================================

def tetrahedron_volumes(
    nodes_cm: np.ndarray,
    tetrahedra: np.ndarray,
) -> np.ndarray:
    """
    Calculate tetrahedron volumes in cm^3.

    This follows the reference calculation:

        V = |(b-a) dot ((c-a) x (d-a))| / 6
    """

    a = nodes_cm[tetrahedra[:, 0]]
    b = nodes_cm[tetrahedra[:, 1]]
    c = nodes_cm[tetrahedra[:, 2]]
    d = nodes_cm[tetrahedra[:, 3]]

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

    return volumes


def load_mesh(
    node_file: Path,
    element_file: Path,
) -> MeshData:
    """Load NODE and ELE files using the reference convention."""

    nodes_cm = read_node_file(node_file)

    tetrahedra, organ_ids = read_ele_file(
        element_file
    )

    if (
        tetrahedra.min() < 0
        or tetrahedra.max() >= len(nodes_cm)
    ):
        raise ValueError(
            f"Tetrahedron node indices in {element_file} "
            "are inconsistent with the NODE file."
        )

    return MeshData(
        nodes_cm=nodes_cm,
        tetrahedra=tetrahedra,
        organ_ids=organ_ids,
    )


def calculate_tetrahedron_masses(
    mesh: MeshData,
    densities: dict[int, float],
) -> np.ndarray:
    """
    Calculate tetrahedron masses in grams.

    NODE coordinates are kept in cm, so:

        volume = cm^3
        density = g/cm^3
        mass = volume * density
    """

    volumes_cm3 = tetrahedron_volumes(
        mesh.nodes_cm,
        mesh.tetrahedra,
    )

    density_array = np.empty(
        len(mesh.organ_ids),
        dtype=np.float64,
    )

    missing = []

    for index, organ_id in enumerate(
        mesh.organ_ids
    ):

        organ_id = int(organ_id)

        if organ_id not in densities:
            missing.append(organ_id)
        else:
            density_array[index] = densities[organ_id]

    if missing:
        missing = sorted(set(missing))

        raise ValueError(
            "No density was found for organ ID(s): "
            + ", ".join(map(str, missing))
        )

    return volumes_cm3 * density_array


def build_region_samplers(
    mesh: MeshData,
    tetra_masses_g: np.ndarray,
    target_regions: list[dict],
) -> dict[str, RegionSampler]:
    """
    Build mass-weighted samplers.

    For a compound region:
        P(organ ID) = organ mass / total region mass

    Then:
        P(tetrahedron | organ) = tetrahedron mass / organ mass

    Since each organ has a fixed density, this is equivalent to
    volume-weighted sampling within an individual organ while ensuring
    that compound-organ components are selected according to their
    mass fractions.
    """

    samplers = {}

    unique_ids = np.unique(
        mesh.organ_ids
    )

    for region in target_regions:

        name = region["name"]
        organ_ids = tuple(
            region["ids"]
        )

        missing_ids = sorted(
            set(organ_ids)
            - set(unique_ids.tolist())
        )

        if missing_ids:
            raise ValueError(
                f"Region '{name}' contains organ ID(s) not found "
                f"in the mesh: {missing_ids}"
            )

        region_tet_indices = np.flatnonzero(
            np.isin(
                mesh.organ_ids,
                np.asarray(organ_ids),
            )
        )

        if len(region_tet_indices) == 0:
            raise ValueError(
                f"No tetrahedra were found for region '{name}'."
            )

        weights = tetra_masses_g[
            region_tet_indices
        ]

        total_mass = float(
            weights.sum()
        )

        if total_mass <= 0:
            raise ValueError(
                f"Region '{name}' has zero mass."
            )

        cumulative = np.cumsum(
            weights
        )

        cumulative /= cumulative[-1]

        samplers[name] = RegionSampler(
            region_name=name,
            organ_ids=organ_ids,
            tetrahedra=region_tet_indices,
            cumulative_weights=cumulative,
            total_mass_g=total_mass,
        )

    return samplers


# ======================================================================
# RANDOM POINT SAMPLING
# ======================================================================

def sample_points_in_tetrahedra(
    vertices: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Uniformly sample one point in each tetrahedron.

    Normalized exponential variates generate uniform barycentric
    coordinates over a tetrahedron.
    """

    number_of_points = len(vertices)

    exponential = rng.exponential(
        size=(number_of_points, 4)
    )

    barycentric = (
        exponential
        / exponential.sum(axis=1, keepdims=True)
    )

    # NODE coordinates are stored in cm, while CLD distances are required
    # in mm.
    points_cm = np.einsum(
        "ij,ijk->ik",
        barycentric,
        vertices,
    )

    return points_cm * 10.0


def sample_region_points(
    mesh: MeshData,
    sampler: RegionSampler,
    number_of_points: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample points from a region according to the sampler weights."""

    random_values = rng.random(
        number_of_points
    )

    tetra_positions = np.searchsorted(
        sampler.cumulative_weights,
        random_values,
        side="right",
    )

    tetra_indices = sampler.tetrahedra[
        tetra_positions
    ]

    tetra_nodes = mesh.tetrahedra[
        tetra_indices
    ]

    vertices = mesh.nodes_cm[
        tetra_nodes
    ]

    return sample_points_in_tetrahedra(
        vertices,
        rng,
    )


# ======================================================================
# CLD CALCULATION
# ======================================================================

def calculate_distance_chunk(
    mesh: MeshData,
    source_sampler: RegionSampler,
    target_sampler: RegionSampler,
    number_of_points: int,
    seed: int,
) -> tuple[np.ndarray, float, float]:
    """
    Calculate one chunk of distances.

    Returns:
        histogram
        sum(distance)
        sum(distance^2)
    """

    rng = np.random.default_rng(seed)

    source_points = sample_region_points(
        mesh,
        source_sampler,
        number_of_points,
        rng,
    )

    target_points = sample_region_points(
        mesh,
        target_sampler,
        number_of_points,
        rng,
    )

    differences = (
        source_points
        - target_points
    )

    distances = np.sqrt(
        np.einsum(
            "ij,ij->i",
            differences,
            differences,
        )
    )

    # 1-mm histogram.
    maximum_distance = (
        math.ceil(
            float(distances.max())
            / BIN_WIDTH_MM
        )
        * BIN_WIDTH_MM
    )

    bins = np.arange(
        0.0,
        maximum_distance + BIN_WIDTH_MM,
        BIN_WIDTH_MM,
    )

    histogram, _ = np.histogram(
        distances,
        bins=bins,
    )

    return (
        histogram,
        float(distances.sum()),
        float(np.dot(distances, distances)),
    )


def calculate_cld(
    mesh: MeshData,
    source_sampler: RegionSampler,
    target_sampler: RegionSampler,
    n_samples: int,
    n_threads: int,
    seed: int | None,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Calculate a CLD using multiple worker threads.

    The returned x values are the centers of 1-mm bins.
    The y values are relative numbers:
        count in bin / total number of samples.
    """

    if seed is None:
        seed_sequence = np.random.SeedSequence()
    else:
        seed_sequence = np.random.SeedSequence(seed)

    # Split samples between workers.
    counts = [
        n_samples // n_threads
        for _ in range(n_threads)
    ]

    for index in range(
        n_samples % n_threads
    ):
        counts[index] += 1

    child_sequences = seed_sequence.spawn(
        n_threads
    )

    tasks = [
        (
            mesh,
            source_sampler,
            target_sampler,
            count,
            int(child.generate_state(1)[0]),
        )
        for child, count in zip(
            child_sequences,
            counts,
        )
        if count > 0
    ]

    histograms = []
    sums = []
    sums_squared = []

    with ThreadPoolExecutor(
        max_workers=n_threads
    ) as executor:

        futures = [
            executor.submit(
                calculate_distance_chunk,
                *task,
            )
            for task in tasks
        ]

        for future in futures:
            histogram, distance_sum, distance_sum_sq = (
                future.result()
            )

            histograms.append(histogram)
            sums.append(distance_sum)
            sums_squared.append(distance_sum_sq)

    # Workers can have different maximum distances.
    maximum_bins = max(
        len(histogram)
        for histogram in histograms
    )

    combined_histogram = np.zeros(
        maximum_bins,
        dtype=np.int64,
    )

    for histogram in histograms:
        combined_histogram[
            :len(histogram)
        ] += histogram

    total = combined_histogram.sum()

    if total != n_samples:
        raise RuntimeError(
            f"Histogram contains {total} samples, "
            f"expected {n_samples}."
        )

    relative_number = (
        combined_histogram
        / n_samples
    )

    bin_centers = (
        np.arange(
            maximum_bins,
            dtype=float,
        )
        + 0.5
    ) * BIN_WIDTH_MM

    distance_sum = sum(sums)
    distance_sum_squared = sum(
        sums_squared
    )

    mean = distance_sum / n_samples

    # Population standard deviation of the sampled distances.
    variance = (
        distance_sum_squared / n_samples
        - mean**2
    )

    # Protect against tiny negative round-off.
    variance = max(
        0.0,
        variance,
    )

    standard_deviation = math.sqrt(
        variance
    )

    return (
        bin_centers,
        relative_number,
        mean,
        standard_deviation,
    )


# ======================================================================
# CSV OUTPUT
# ======================================================================

def safe_filename(text: str) -> str:
    """Convert a region/phantom name into a filesystem-safe name."""

    text = str(text).strip()
    text = re.sub(
        r"[^\w\-.]+",
        "_",
        text,
    )
    return text.strip("_") or "unnamed"


def save_cld_csv(
    output_file: Path,
    phantom: PhantomInfo,
    source_region: dict,
    target_region: str,
    distances_mm: np.ndarray,
    relative_number: np.ndarray,
    mean_mm: float,
    std_mm: float,
):
    """Save one CLD to CSV."""

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    source_ids_text = "_".join(
        str(organ_id)
        for organ_id in source_region["ids"]
    )

    df = pd.DataFrame(
        {
            "Phantom": phantom.display_name,
            "Phantom Code": phantom.code,
            "Sex": (
                "Male"
                if phantom.sex == "AM"
                else "Female"
            ),
            "Source Region": source_region["name"],
            "Source Acronym": source_region["acronym"],
            "Source Organ ID(s)": source_ids_text,
            "Target Region": target_region,
            "Distance (mm)": distances_mm,
            "Relative Number": relative_number,
            "Mean Chord Length (mm)": mean_mm,
            "Standard Deviation (mm)": std_mm,
        }
    )

    # The phantom/source/target metadata and the summary statistics are
    # constant for the entire CLD. Keep them only on the first row so
    # the CSV is easier to read without changing the distance distribution.
    repeated_columns = [
        "Phantom",
        "Phantom Code",
        "Sex",
        "Source Region",
        "Source Acronym",
        "Source Organ ID(s)",
        "Target Region",
        "Mean Chord Length (mm)",
        "Standard Deviation (mm)",
    ]

    if len(df) > 1:
        # Cast these columns to object first so blank strings can be written
        # without changing the numeric columns or raising dtype warnings.
        df[repeated_columns] = df[repeated_columns].astype(object)
        df.loc[1:, repeated_columns] = ""

    df.to_csv(
        output_file,
        index=False,
    )


def save_summary_csv(
    rows: list[dict],
    output_file: Path,
):
    """Save mean/std summary for all calculated CLDs."""

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.sort_values(
            [
                "Source Region",
                "Target Region",
                "Phantom",
            ]
        ).reset_index(drop=True)

    df.to_csv(
        output_file,
        index=False,
    )


def cld_csv_path(
    output_dir: Path,
    phantom_code: str,
    source_region: dict,
    target_region: str,
) -> Path:
    """Return the expected CSV path for one phantom/source/target CLD."""

    source_label = (
        source_region["acronym"]
        or source_region["name"]
        or "_".join(map(str, source_region["ids"]))
    )

    csv_name = (
        f"{safe_filename(phantom_code)}_"
        f"source_{safe_filename(source_label)}_"
        f"target_{safe_filename(target_region)}_"
        f"cld.csv"
    )

    return output_dir / "csv" / csv_name


def existing_cld_csvs(
    output_dir: Path,
    resolved_groups: dict[str, list[PhantomInfo]],
    source_regions: list[dict],
    target_regions_by_group: dict[str, list[dict]],
) -> dict[tuple[str, str, str, str], Path]:
    """
    Find existing CLD CSVs for the requested source/target combinations.

    The key is:
        (group_name, phantom_code, source_key, target_region)

    A source region can contain one or many organ IDs.
    """

    existing = {}

    for group_name, phantoms in resolved_groups.items():
        group_target_regions = target_regions_by_group[group_name]

        for phantom in phantoms:
            phantom_target_regions = filter_target_regions_for_sex(
                group_target_regions,
                phantom.sex,
            )

            valid_target_names = {
                region["name"]
                for region in phantom_target_regions
            }

            for source_region in source_regions:
                source_key = source_region["key"]

                for target_region in group_target_regions:
                    target_name = target_region["name"]

                    if target_name not in valid_target_names:
                        continue

                    csv_file = cld_csv_path(
                        output_dir,
                        phantom.code,
                        source_region,
                        target_name,
                    )

                    if csv_file.is_file() and csv_file.stat().st_size > 0:
                        existing[
                            (
                                group_name,
                                phantom.code,
                                source_key,
                                target_name,
                            )
                        ] = csv_file

    return existing


def inspect_cld_csv_status(
    output_dir: Path,
    resolved_groups: dict[str, list[PhantomInfo]],
    source_regions: list[dict],
    target_regions_by_group: dict[str, list[dict]],
) -> tuple[
    set[tuple[str, str, str, str]],
    list[tuple[str, str, str, str]],
]:
    """
    Check which expected CLD CSVs already exist.

    Returns:
        existing_keys
        missing_keys
    """

    existing = existing_cld_csvs(
        output_dir,
        resolved_groups,
        source_regions,
        target_regions_by_group,
    )

    expected_keys = []

    for group_name, phantoms in resolved_groups.items():
        group_target_regions = target_regions_by_group[group_name]

        for phantom in phantoms:
            phantom_target_regions = filter_target_regions_for_sex(
                group_target_regions,
                phantom.sex,
            )

            for source_region in source_regions:
                source_key = source_region["key"]

                for target_region in phantom_target_regions:
                    expected_keys.append(
                        (
                            group_name,
                            phantom.code,
                            source_key,
                            target_region["name"],
                        )
                    )

    expected_keys = set(expected_keys)
    existing_keys = set(existing)
    missing_keys = sorted(expected_keys - existing_keys)

    return existing_keys, missing_keys


def ask_existing_results_action(
    existing_count: int,
) -> str:
    """
    Ask whether to skip or redo calculations when all expected CSVs exist.

    Returns:
        "skip" or "redo"
    """

    print()
    print("=" * 70)
    print("EXISTING CLD RESULTS")
    print("=" * 70)
    print(
        f"All {existing_count:,} expected CLD CSV results are already present."
    )
    print()
    print("[S] Skip recalculation and use the existing CSV results")
    print("[R] Redo all CLD calculations and overwrite the existing CSVs")

    while True:
        choice = input("\nChoose S or R: ").strip().upper()

        if choice == "S":
            return "skip"

        if choice == "R":
            return "redo"

        print("Invalid choice. Enter S to skip or R to redo.")


def load_cld_csv(
    csv_file: Path,
    phantom: PhantomInfo,
    source_key: str,
    target_region: str,
) -> dict:
    """
    Load one existing CLD CSV into the same structure used by new results.
    """

    required_columns = {
        "Distance (mm)",
        "Relative Number",
        "Mean Chord Length (mm)",
        "Standard Deviation (mm)",
    }

    df = clean_columns(pd.read_csv(csv_file))

    missing_columns = sorted(
        required_columns - set(df.columns)
    )

    if missing_columns:
        raise ValueError(
            f"Existing CLD CSV is missing required column(s): "
            f"{missing_columns}\n{csv_file}"
        )

    if df.empty:
        raise ValueError(
            f"Existing CLD CSV is empty:\n{csv_file}"
        )

    distances_mm = pd.to_numeric(
        df["Distance (mm)"],
        errors="coerce",
    ).to_numpy(dtype=float)

    relative_number = pd.to_numeric(
        df["Relative Number"],
        errors="coerce",
    ).to_numpy(dtype=float)

    if (
        not np.isfinite(distances_mm).all()
        or not np.isfinite(relative_number).all()
    ):
        raise ValueError(
            f"Existing CLD CSV contains invalid distance or "
            f"relative-number values:\n{csv_file}"
        )

    mean_values = pd.to_numeric(
        df["Mean Chord Length (mm)"],
        errors="coerce",
    ).dropna()

    std_values = pd.to_numeric(
        df["Standard Deviation (mm)"],
        errors="coerce",
    ).dropna()

    if mean_values.empty or std_values.empty:
        raise ValueError(
            f"Existing CLD CSV does not contain mean/std values "
            f"in its metadata row:\n{csv_file}"
        )

    mean_mm = float(mean_values.iloc[0])
    std_mm = float(std_values.iloc[0])

    return {
        "phantom": phantom,
        "distance_mm": distances_mm,
        "relative_number": relative_number,
        "mean_mm": mean_mm,
        "std_mm": std_mm,
    }


def load_existing_cld_results(
    existing_files: dict[tuple[str, str, str, str], Path],
    resolved_groups: dict[str, list[PhantomInfo]],
) -> dict:
    """
    Load all existing CLD CSVs into the plotting/calculation result structure.
    """

    phantom_lookup = {
        phantom.code: phantom
        for phantoms in resolved_groups.values()
        for phantom in phantoms
    }

    cld_results = {
        group_name: {}
        for group_name in resolved_groups
    }

    for (
        group_name,
        phantom_code,
        source_key,
        target_name,
    ), csv_file in existing_files.items():

        phantom = phantom_lookup[phantom_code]

        result = load_cld_csv(
            csv_file,
            phantom,
            source_key,
            target_name,
        )

        cld_results[group_name].setdefault(source_key, {})
        cld_results[group_name][source_key].setdefault(target_name, {})
        cld_results[group_name][source_key][target_name][phantom_code] = result

    return cld_results


# ======================================================================
# PLOTTING
# ======================================================================
def plot_x_axis_limit(
    results: list[dict],
) -> float:
    """Return an automatic x-axis maximum rounded up to the nearest 50 mm."""

    finite_maxima = []
    for result in results:
        distances = np.asarray(result["distance_mm"], dtype=float)
        finite_distances = distances[np.isfinite(distances)]
        if finite_distances.size:
            finite_maxima.append(float(finite_distances.max()))

    if not finite_maxima:
        return PLOT_X_ROUNDING_MM

    data_max = max(finite_maxima)
    return math.ceil(data_max / PLOT_X_ROUNDING_MM) * PLOT_X_ROUNDING_MM


def plot_selected_clds(
    plot_pairs: list[tuple[int, str]],
    source_names: dict[str, str],
    cld_results: dict,
    output_file: Path,
    n_rows: int,
    n_cols: int,
):
    """
    Create one multi-panel PDF containing the user-selected source/target
    CLD pairs.

    Phantom colors and line styles are shared across every subplot. The
    phantom legend is placed once above the whole figure, while each subplot
    has its own statistics legend so the line sample remains directly
    associated with its mean ± standard deviation.
    """

    # Build the plotting order directly from the phantom tuples in the
    # registry. This makes the registry the authoritative source of legend
    # order and phantom styling.
    ordered_phantoms = []

    for group_code in PHANTOM_GROUPS_TO_COMPARE:
        group_name = registry.PHANTOM_GROUPS[group_code].display_name

        for phantom_code in registry.get_phantom_group(group_code):
            phantom = registry.get_phantom(phantom_code)
            ordered_phantoms.append(
                (
                    group_name,
                    phantom_code,
                    phantom,
                )
            )

    if not ordered_phantoms or not plot_pairs:
        return

    # Use a slightly wider overall figure so the panels are not cramped
    # against the left/right edges, while keeping the vertical layout compact.
    panel_width = FIGURE_WIDTH
    panel_height = FIGURE_WIDTH * 0.70
    figure, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(
            panel_width * n_cols,
            panel_height * n_rows,
        ),
        squeeze=False,
    )
    axes = axes.ravel()

    # Assign colors and line styles automatically from Matplotlib's built-in
    # cycles. The order follows the registry's phantom tuples.
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color",
        [],
    )

    if not color_cycle:
        color_cycle = [
            f"C{i}"
            for i in range(len(ordered_phantoms))
        ]

    line_style_cycle = [
        "-",
        "--",
        "-.",
        ":",
    ]

    phantom_colors = {
        phantom_code: color_cycle[
            index % len(color_cycle)
        ]
        for index, (_, phantom_code, _) in enumerate(
            ordered_phantoms
        )
    }

    phantom_line_styles = {
        phantom_code: line_style_cycle[
            index % len(line_style_cycle)
        ]
        for index, (_, phantom_code, _) in enumerate(
            ordered_phantoms
        )
    }

    # Draw every selected source-target pair.
    used_handles = []
    used_labels = []

    for axis, (source_key, target_region) in zip(
        axes,
        plot_pairs,
    ):
        source_name = source_names.get(
            source_key,
            f"Organ {source_key}",
        )

        subplot_handles = []
        statistics_labels = []
        subplot_results = []

        for (
            group_name,
            phantom_code,
            phantom,
        ) in ordered_phantoms:

            # Gonads is a display target. For each phantom, plot the
            # sex-specific CSV that actually represents that phantom:
            # female -> Ovaries, male -> Testes. Do not plot a mixed
            # Ovaries+Testes distribution under the Gonads label.
            csv_target_region = target_region
            if target_region == "Gonads":
                csv_target_region = (
                    "Ovaries" if phantom.sex == "AF" else "Testes"
                )

            result = (
                cld_results
                .get(group_name, {})
                .get(source_key, {})
                .get(csv_target_region, {})
                .get(phantom_code)
            )

            if result is None:
                continue

            line, = axis.plot(
                result["distance_mm"],
                result["relative_number"],
                color=phantom_colors[phantom_code],
                linestyle=phantom_line_styles[phantom_code],
                linewidth=1.2,
            )

            subplot_handles.append(line)
            subplot_results.append(result)
            statistics_labels.append(
                f"{result['mean_mm']:.2f} ± "
                f"{result['std_mm']:.2f} mm"
            )

            if phantom_code not in used_labels:
                used_handles.append(line)
                used_labels.append(phantom_code)

        # Show the target region inside each panel. If multiple source
        # organs are selected, include the source name as well so every
        # panel remains unambiguous.
        if len({
            source_key
            for source_key, _ in plot_pairs
        }) == 1:
            panel_label = target_region
        else:
            panel_label = f"{source_name} → {target_region}"

        axis.text(
            0.02,
            0.96,
            panel_label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )

        # Set the x-axis automatically from the largest sampled distance,
        # rounded upward to the nearest 50 mm.
        axis.set_xlim(
            0.0,
            plot_x_axis_limit(subplot_results),
        )

        axis.grid(
            alpha=0.18,
            linewidth=0.5,
        )

        # Keep the line samples in this legend so readers can directly
        # associate each mean ± SD with the corresponding curve.
        if subplot_handles:
            axis.legend(
                subplot_handles,
                statistics_labels,
                fontsize=5.5,
                frameon=False,
                loc="best",
                handletextpad=0.35,
                handlelength=1.8,
                borderaxespad=0.2,
            )

    # Hide unused panels.
    for axis in axes[len(plot_pairs):]:
        axis.set_visible(False)

    # Determine the source title for the whole figure.
    # If all selected panels use the same source organ, show one concise
    # title such as "Source: Liver". If multiple source organs are selected,
    # list them together so the title remains unambiguous.
    selected_source_keys = list(dict.fromkeys(
        source_key
        for source_key, _ in plot_pairs
    ))
    selected_source_names = [
        source_names.get(
            source_key,
            f"Organ {source_key}",
        )
        for source_key in selected_source_keys
    ]

    source_title = "Source: " + ", ".join(selected_source_names)

    # Put the source title at the top of the entire figure.
    figure.suptitle(
        source_title,
        fontsize=10,
        y=0.995,
    )

    # One shared phantom legend for the entire figure, outside the plots
    # and directly below the source title.
    if used_handles:
        figure.legend(
            used_handles,
            used_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.970),
            ncol=len(used_labels),
            fontsize=6,
            frameon=False,
            handletextpad=0.25,
            columnspacing=0.8,
        )

    # One shared x/y label for the entire multi-panel figure.
    figure.supxlabel(
        "Distance (mm)",
        fontsize=7,
        y=0.015,
    )
    figure.supylabel(
        "Relative Number",
        fontsize=7,
        x=0.015,
    )

    # Compact spacing: reserve only the space needed for the title, shared
    # phantom legend, and common axis labels.
    figure.subplots_adjust(
        left=0.0975,
        right=0.975,
        bottom=0.060,
        top=0.925,
        wspace=0.2,
        hspace=0.2,
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pdf_file = output_file.with_suffix(".pdf")

    figure.savefig(
        pdf_file,
        format="pdf",
        bbox_inches=None,
    )

    plt.close(figure)


def ask_plot_layout() -> tuple[int, int]:
    """Ask the user for the subplot layout, e.g. 3x2 or 2x4."""

    print()
    print("=" * 70)
    print("PLOTTING FORMAT")
    print("=" * 70)
    print("Enter the number of rows and columns.")
    print("Examples: 3x2, 2x4, 1x3")

    while True:
        choice = input("\nPlotting format [rows x columns]: ").strip().lower()
        match = re.fullmatch(r"(\d+)\s*[x×]\s*(\d+)", choice)

        if match is None:
            print("Invalid format. Enter something like 3x2 or 2x4.")
            continue

        n_rows = int(match.group(1))
        n_cols = int(match.group(2))

        if n_rows < 1 or n_cols < 1:
            print("Rows and columns must both be at least 1.")
            continue

        return n_rows, n_cols


def available_plot_pairs(
    cld_results: dict,
    source_keys: list[str],
    source_names: dict[str, str],
    target_regions: list[dict],
) -> list[tuple[int, str]]:
    """Return source-target pairs for which at least one CLD CSV exists."""

    pairs = []

    for source_key in source_keys:
        for region in target_regions:
            target_name = region["name"]

            found = False

            for group_name, group_results in cld_results.items():
                if target_name == "Gonads":
                    # A Gonads plot is available when at least one
                    # sex-appropriate Ovaries/Testes CSV exists.
                    entries = {}

                    for phantom_code, result in (
                        group_results
                        .get(source_key, {})
                        .get("Ovaries", {})
                        .items()
                    ):
                        if result["phantom"].sex == "AF":
                            entries[phantom_code] = result

                    for phantom_code, result in (
                        group_results
                        .get(source_key, {})
                        .get("Testes", {})
                        .items()
                    ):
                        if result["phantom"].sex == "AM":
                            entries[phantom_code] = result
                else:
                    entries = (
                        group_results
                        .get(source_key, {})
                        .get(target_name, {})
                    )

                if entries:
                    found = True
                    break

            if found:
                pairs.append(
                    (
                        source_key,
                        target_name,
                    )
                )

    return pairs


def select_plot_pairs(
    pairs: list[tuple[int, str]],
    source_names: dict[str, str],
    n_rows: int,
    n_cols: int,
) -> list[tuple[int, str]]:
    """
    Ask which generated source-target CLDs should be included in the PDF.

    Only pairs with available CSV-backed CLD results are shown.
    """

    if not pairs:
        return []

    capacity = n_rows * n_cols

    print()
    print("=" * 70)
    print("PLOT SELECTION")
    print("=" * 70)
    print(
        "Select which source-target pairs to include in the PDF."
    )
    print(
        f"Your {n_rows}x{n_cols} layout can contain up to "
        f"{capacity} plots."
    )

    for index, (source_key, target_name) in enumerate(
        pairs,
        start=1,
    ):
        source_name = source_names.get(
            source_key,
            f"Organ {source_key}",
        )
        print(
            f"[{index}] {source_name} ({source_key}) → "
            f"{target_name}"
        )

    print("[A] All available source-target pairs")

    while True:
        choice = input(
            "\nSelect plot(s): "
        ).strip().upper()

        if choice == "A":
            selected = pairs
        else:
            try:
                indices = [
                    int(value.strip()) - 1
                    for value in choice.split(",")
                ]

                if not indices:
                    raise ValueError

                if any(
                    index < 0 or index >= len(pairs)
                    for index in indices
                ):
                    raise ValueError

                selected = list(
                    dict.fromkeys(
                        pairs[index]
                        for index in indices
                    )
                )

            except ValueError:
                print(
                    "Invalid selection. Enter plot numbers "
                    "separated by commas, or A."
                )
                continue

        if len(selected) > capacity:
            print(
                f"You selected {len(selected)} plots, but the "
                f"{n_rows}x{n_cols} layout only has {capacity} panels."
            )
            print(
                "Choose a larger plotting format or select fewer plots."
            )
            continue

        return selected


def create_selected_plot_pdf(
    cld_results: dict,
    source_keys: list[str],
    source_names: dict[str, str],
    target_regions: list[dict],
    output_dir: Path,
) -> None:
    """Ask for the layout/pairs and create one selected multi-panel PDF."""

    pairs = available_plot_pairs(
        cld_results,
        source_keys,
        source_names,
        target_regions,
    )

    if not pairs:
        print("No CSV-backed CLD results are available for plotting.")
        return

    n_rows, n_cols = ask_plot_layout()

    selected_pairs = select_plot_pairs(
        pairs,
        source_names,
        n_rows,
        n_cols,
    )

    if not selected_pairs:
        print("No plots selected.")
        return

    figure_name = (
        f"selected_clds_{n_rows}x{n_cols}.pdf"
    )

    figure_base = (
        output_dir
        / "figures"
        / figure_name
    )

    plot_selected_clds(
        selected_pairs,
        source_names,
        cld_results,
        figure_base,
        n_rows,
        n_cols,
    )

    print()
    print("Saved selected journal figure:")
    print(f"  {figure_base}")


# ======================================================================
# USER INTERFACE
# ======================================================================

def select_source_organs(
    source_keys: list[str],
    source_names: dict[str, str],
) -> list[str]:

    if not ASK_SOURCE_SELECTION:
        return source_keys

    print()
    print("=" * 70)
    print("SOURCE ORGAN SELECTION")
    print("=" * 70)

    for index, source_key in enumerate(
        source_keys,
        start=1,
    ):
        name = source_names.get(
            source_key,
            "Unknown",
        )

        print(
            f"[{index}] {source_key} ({name})"
        )

    print("[A] All source organs")

    while True:

        choice = input(
            "\nSelect source organ(s): "
        ).strip().upper()

        if choice == "A":
            return source_keys

        try:
            indices = [
                int(value.strip()) - 1
                for value in choice.split(",")
            ]

            if not indices:
                raise ValueError

            if any(
                index < 0
                or index >= len(source_keys)
                for index in indices
            ):
                raise ValueError

            selected = list(
                dict.fromkeys(
                    source_keys[index]
                    for index in indices
                )
            )

            return selected

        except ValueError:
            print(
                "Invalid selection. "
                "Enter numbers separated by commas, "
                "or A."
            )


# ======================================================================
# PHANTOM RESOLUTION
# ======================================================================

def resolve_phantoms(
    registry,
) -> dict[str, list[PhantomInfo]]:
    """
    Resolve phantom groups through b_phantom_registry.py.

    The registry is the authoritative source for NODE/ELE/MATERIAL
    file paths. This function only converts PhantomSpec objects into
    the local PhantomInfo representation used by the CLD calculation.
    """

    resolved = {}

    for group_code in PHANTOM_GROUPS_TO_COMPARE:

        phantom_codes = registry.get_phantom_group(
            group_code
        )

        group = registry.PHANTOM_GROUPS[group_code]

        phantoms = []

        for phantom_code in phantom_codes:

            spec = registry.get_phantom(
                phantom_code
            )

            phantoms.append(
                PhantomInfo(
                    code=spec.code,
                    display_name=spec.display_name,
                    sex=spec.sex,
                    node_file=Path(spec.node_file),
                    element_file=Path(spec.element_file),
                    material_file=Path(spec.material_file),
                )
            )

        resolved[group.display_name] = phantoms

    return resolved

def validate_phantom_files(
    resolved_groups: dict[str, list[PhantomInfo]],
) -> None:
    """Check all NODE/ELE/MATERIAL files before starting sampling."""

    missing = []

    for group_name, phantoms in resolved_groups.items():

        for phantom in phantoms:

            for label, filename in (
                ("NODE", phantom.node_file),
                ("ELE", phantom.element_file),
                ("MATERIAL", phantom.material_file),
            ):

                if not filename.is_file():
                    missing.append(
                        f"{phantom.display_name}: "
                        f"{label} -> {filename}"
                    )

    if missing:
        raise FileNotFoundError(
            "The following phantom files could not be found:\n"
            + "\n".join(
                f"  - {item}"
                for item in missing
            )
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    selected_font = configure_journal_font()

    print()
    print(f"Figure font: {selected_font}")
    print("=" * 70)
    print("CHORD-LENGTH DISTRIBUTION CALCULATION")
    print("=" * 70)

    source_csv = Path(config.SOURCE_CSV_CLD)

    output_dir = (
        Path(config.RESULTS_DIR)
        / OUTPUT_DIR_NAME
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    n_threads = ask_threads()

    print()
    print(f"SOURCE_CSV:")
    print(f"  {source_csv}")

    print()
    print(f"Point pairs per CLD: {N_SAMPLES:,}")
    print(f"Bin width: {BIN_WIDTH_MM:g} mm")
    print(f"Threads: {n_threads}")

    # --------------------------------------------------------------
    # Load source and target definitions
    # --------------------------------------------------------------

    source_regions = load_source_regions(
        source_csv
    )

    # Each phantom group uses its own target-region mapping from the registry.
    # The mappings do not need to be identical.
    target_regions_by_group_code: dict[str, list[dict]] = {}
    for group_code in PHANTOM_GROUPS_TO_COMPARE:
        target_region_file = registry.get_target_region_file(group_code)
        target_regions_by_group_code[group_code] = merge_gonads_target_regions(
            load_target_regions(Path(target_region_file))
        )

    # Resolve groups so target-region definitions can use the same display
    # names as cld_results.
    resolved_groups = resolve_phantoms(registry)
    target_regions_by_group = {
        registry.PHANTOM_GROUPS[group_code].display_name:
        target_regions_by_group_code[group_code]
        for group_code in PHANTOM_GROUPS_TO_COMPARE
    }

    # Union of target-region names is used only for plotting/summary ordering.
    target_regions = []
    target_names_seen = set()
    for regions in target_regions_by_group.values():
        for region in regions:
            if region["name"] not in target_names_seen:
                target_regions.append(region)
                target_names_seen.add(region["name"])

    # SOURCE_CSV is authoritative for source-region names and IDs.
    # Legacy one-ID source CSVs use the organ database only to recover the
    # human-readable organ name for plotting/summary output.
    source_names: dict[str, str] = {}
    source_regions_by_key = {}

    for source_region in source_regions:
        source_key = source_region["key"]
        source_name = source_region["name"]

        if not source_name:
            source_id = source_region["ids"][0]
            source_name = None

            for phantom_organs in organ_database.ORGANS.values():
                for organ_id, organ_data in phantom_organs.items():
                    if int(organ_id) == source_id:
                        source_name = str(organ_data["name"])
                        break
                if source_name is not None:
                    break

            if source_name is None:
                source_name = f"Organ {source_id}"

        source_names[source_key] = source_name
        source_regions_by_key[source_key] = source_region

    # All source regions are always calculated. The user chooses only which
    # CSV-backed source-target pairs to include in the PDF later.
    selected_source_keys = [
        source_region["key"]
        for source_region in source_regions
    ]

    # --------------------------------------------------------------
    # Validate resolved phantom files
    # --------------------------------------------------------------

    validate_phantom_files(
        resolved_groups
    )

    # --------------------------------------------------------------
    # Check existing CLD CSV results BEFORE loading the meshes.
    # --------------------------------------------------------------

    existing_keys, missing_keys = inspect_cld_csv_status(
        output_dir,
        resolved_groups,
        source_regions,
        target_regions_by_group,
    )

    total_expected = len(existing_keys) + len(missing_keys)

    if total_expected == 0:
        print()
        print("No valid source-target CLD combinations were found.")
        return

    print()
    print("=" * 70)
    print("CLD CSV STATUS")
    print("=" * 70)
    print(f"Expected CLD CSVs: {total_expected:,}")
    print(f"Already present:   {len(existing_keys):,}")
    print(f"Missing:            {len(missing_keys):,}")

    if not missing_keys:
        action = ask_existing_results_action(
            len(existing_keys)
        )

        if action == "skip":
            # Use all existing CSVs to regenerate the figures/summary
            # without repeating the expensive Monte Carlo sampling.
            existing_files = existing_cld_csvs(
                output_dir,
                resolved_groups,
                source_regions,
                target_regions_by_group,
            )

            cld_results = load_existing_cld_results(
                existing_files,
                resolved_groups,
            )

            summary_rows = []

            for (
                group_name,
                phantom_code,
                source_key,
                target_name,
            ), csv_file in existing_files.items():

                result = (
                    cld_results
                    .get(group_name, {})
                    .get(source_key, {})
                    .get(target_name, {})
                    .get(phantom_code)
                )

                if result is None:
                    continue

                phantom = result["phantom"]
                source_region = source_regions_by_key[source_key]
                source_name = source_names[source_key]

                summary_rows.append(
                    {
                        "Phantom": phantom.display_name,
                        "Phantom Code": phantom.code,
                        "Sex": (
                            "Male"
                            if phantom.sex == "AM"
                            else "Female"
                        ),
                        "Source Region": source_name,
                        "Source Acronym": source_region["acronym"],
                        "Source Organ ID(s)": "_".join(
                            str(organ_id)
                            for organ_id in source_region["ids"]
                        ),
                        "Target Region": target_name,
                        "Mean Chord Length (mm)": result["mean_mm"],
                        "Standard Deviation (mm)": result["std_mm"],
                        "Number of Point Pairs": N_SAMPLES,
                        "Bin Width (mm)": BIN_WIDTH_MM,
                    }
                )

            # Regenerate only the user-selected multi-panel PDF from
            # the existing CSV results.
            create_selected_plot_pdf(
                cld_results,
                selected_source_keys,
                source_names,
                target_regions,
                output_dir,
            )

            summary_file = (
                output_dir
                / "csv"
                / "cld_summary.csv"
            )

            save_summary_csv(
                summary_rows,
                summary_file,
            )

            print()
            print("=" * 70)
            print("EXISTING CLD RESULTS USED")
            print("=" * 70)
            print(f"CSV files: {output_dir / 'csv'}")
            print(f"Journal figures: {output_dir / 'figures'}")
            print(f"Summary:   {summary_file}")

            return

        # Redo means every expected result is recalculated.
        missing_keys = sorted(
            existing_keys | set(missing_keys)
        )

    else:
        print()
        print(
            f"Only the {len(missing_keys):,} missing CLD CSV(s) "
            "will be calculated."
        )

    # --------------------------------------------------------------
    # Load existing CSVs first. These are retained and combined with
    # newly calculated CLDs for plotting and the summary.
    # --------------------------------------------------------------

    existing_files = existing_cld_csvs(
        output_dir,
        resolved_groups,
        source_regions,
        target_regions_by_group,
    )

    # In redo mode, do not load the existing results because every
    # expected combination will be replaced.
    if not missing_keys or set(missing_keys) != set(existing_files):
        existing_files = {
            key: path
            for key, path in existing_files.items()
            if key not in set(missing_keys)
        }

    cld_results = load_existing_cld_results(
        existing_files,
        resolved_groups,
    )

    # --------------------------------------------------------------
    # Load all meshes once.
    # --------------------------------------------------------------

    meshes = {}
    samplers = {}
    source_samplers = {}

    for group_name, phantoms in resolved_groups.items():

        print()
        print("-" * 70)
        print(f"Loading phantom group: {group_name}")
        print("-" * 70)

        for phantom in phantoms:

            print(
                f"  Loading {phantom.display_name}..."
            )

            mesh = load_mesh(
                phantom.node_file,
                phantom.element_file,
            )

            densities = read_material_densities(
                phantom.material_file
            )

            tetra_masses = calculate_tetrahedron_masses(
                mesh,
                densities,
            )

            group_code_for_phantom = next(
                group_code
                for group_code in PHANTOM_GROUPS_TO_COMPARE
                if phantom.code in registry.get_phantom_group(group_code)
            )
            group_name_for_phantom = registry.PHANTOM_GROUPS[
                group_code_for_phantom
            ].display_name
            phantom_target_regions = filter_target_regions_for_sex(
                target_regions_by_group[group_name_for_phantom],
                phantom.sex,
            )

            region_samplers = build_region_samplers(
                mesh,
                tetra_masses,
                phantom_target_regions,
            )

            meshes[phantom.code] = mesh
            samplers[phantom.code] = region_samplers

            # Pre-build samplers for every source region so that
            # density/volume processing is not repeated for every
            # target region.
            source_samplers[phantom.code] = {}

            mesh_organ_ids = set(
                int(organ_id)
                for organ_id in np.unique(mesh.organ_ids)
            )

            for source_key in selected_source_keys:
                source_region = source_regions_by_key[source_key]

                # Apply the same sex-specific handling used by target
                # regions. This also allows a compound source such as
                # Gonads to resolve to Testes in a male phantom or
                # Ovaries in a female phantom.
                phantom_source_regions = filter_target_regions_for_sex(
                    [source_region],
                    phantom.sex,
                )

                if not phantom_source_regions:
                    continue

                source_region_for_phantom = phantom_source_regions[0]
                missing_ids = sorted(
                    set(source_region_for_phantom["ids"])
                    - mesh_organ_ids
                )

                if missing_ids:
                    print(
                        f"    [SKIP] {phantom.display_name}: source region "
                        f"'{source_region['name']}' contains organ ID(s) "
                        f"not present in mesh: {missing_ids}"
                    )
                    continue

                source_samplers[phantom.code][source_key] = (
                    build_region_samplers(
                        mesh,
                        tetra_masses,
                        [source_region_for_phantom],
                    )[source_region_for_phantom["name"]]
                )

            print(
                f"    Nodes: {len(mesh.nodes_cm):,}"
            )
            print(
                f"    Tetrahedra: {len(mesh.tetrahedra):,}"
            )

    # --------------------------------------------------------------
    # Calculate only missing CLDs.
    # --------------------------------------------------------------

    missing_key_set = set(missing_keys)

    for source_key in selected_source_keys:

        source_region = source_regions_by_key[source_key]
        source_name = source_names[source_key]
        source_ids_text = "_".join(
            str(organ_id)
            for organ_id in source_region["ids"]
        )

        print()
        print("=" * 70)
        print(
            f"SOURCE: {source_name} "
            f"(IDs {source_ids_text})"
        )
        print("=" * 70)

        for group_name, phantoms in resolved_groups.items():

            print()
            print(
                f"Phantom group: {group_name}"
            )

            group_target_regions = target_regions_by_group[group_name]

            # Ensure every target dictionary exists, including those
            # loaded from pre-existing CSVs.
            for target_region in group_target_regions:
                target_name = target_region["name"]
                cld_results.setdefault(
                    group_name,
                    {},
                ).setdefault(
                    source_key,
                    {},
                ).setdefault(
                    target_name,
                    {},
                )

            for target_region in group_target_regions:

                target_name = target_region["name"]

                print(
                    f"  Target: {target_name}"
                )

                for phantom in phantoms:

                    key = (
                        group_name,
                        phantom.code,
                        source_key,
                        target_name,
                    )

                    # Existing result is retained unless redo/missing
                    # logic explicitly marked this key for calculation.
                    if key not in missing_key_set:
                        continue

                    mesh = meshes[
                        phantom.code
                    ]

                    phantom_samplers = samplers[
                        phantom.code
                    ]

                    # The source sampler was built from the complete
                    # compound source region. If any required source IDs
                    # were absent in this phantom, the sampler is omitted.
                    source_sampler = (
                        source_samplers
                        .get(phantom.code, {})
                        .get(source_key)
                    )

                    if source_sampler is None:
                        print(
                            f"    [SKIP] {phantom.display_name}: "
                            f"source region '{source_name}' "
                            f"({source_ids_text}) is not available."
                        )
                        continue

                    if target_name not in phantom_samplers:
                        # Sex-specific target region is not present
                        # in this phantom, so skip this combination.
                        continue

                    target_sampler = phantom_samplers[
                        target_name
                    ]

                    # Create a deterministic but distinct seed
                    # for each source/group/target/phantom.
                    stable_key = (
                        f"{source_key}|"
                        f"{group_name}|"
                        f"{target_name}|"
                        f"{phantom.code}|"
                        f"{RANDOM_SEED}"
                    )

                    seed = int.from_bytes(
                        hashlib.sha256(
                            stable_key.encode("utf-8")
                        ).digest()[:4],
                        byteorder="little",
                        signed=False,
                    )

                    (
                        distances_mm,
                        relative_number,
                        mean_mm,
                        std_mm,
                    ) = calculate_cld(
                        mesh,
                        source_sampler,
                        target_sampler,
                        N_SAMPLES,
                        n_threads,
                        seed,
                    )

                    cld_results.setdefault(
                        group_name,
                        {},
                    ).setdefault(
                        source_key,
                        {},
                    ).setdefault(
                        target_name,
                        {},
                    )[phantom.code] = {
                        "phantom": phantom,
                        "distance_mm": distances_mm,
                        "relative_number": relative_number,
                        "mean_mm": mean_mm,
                        "std_mm": std_mm,
                    }

                    # Save CLD CSV.
                    csv_file = cld_csv_path(
                        output_dir,
                        phantom.code,
                        source_region,
                        target_name,
                    )

                    save_cld_csv(
                        csv_file,
                        phantom,
                        source_region,
                        target_name,
                        distances_mm,
                        relative_number,
                        mean_mm,
                        std_mm,
                    )

                    print(
                        f"    {phantom.display_name}: "
                        f"{mean_mm:.2f} ± {std_mm:.2f} mm"
                    )

    # --------------------------------------------------------------
    # Build the summary from every calculated/available CLD.
    # Plotting is handled separately so only user-selected pairs
    # are rendered in the PDF.
    # --------------------------------------------------------------

    summary_rows = []

    for source_key in selected_source_keys:

        source_region = source_regions_by_key[source_key]
        source_name = source_names[source_key]

        for region in target_regions:
            target_name = region["name"]

            # Add every available phantom result for this
            # source/target pair to the summary.
            for group_name, phantoms in resolved_groups.items():
                if target_name not in {
                    region["name"]
                    for region in target_regions_by_group[group_name]
                }:
                    continue

                for phantom in phantoms:
                    result = (
                        cld_results
                        .get(group_name, {})
                        .get(source_key, {})
                        .get(target_name, {})
                        .get(phantom.code)
                    )

                    if result is None:
                        continue

                    summary_rows.append(
                        {
                            "Phantom": phantom.display_name,
                            "Phantom Code": phantom.code,
                            "Sex": (
                                "Male"
                                if phantom.sex == "AM"
                                else "Female"
                            ),
                            "Source Region": source_name,
                            "Source Acronym": source_region["acronym"],
                            "Source Organ ID(s)": "_".join(
                                str(organ_id)
                                for organ_id in source_region["ids"]
                            ),
                            "Target Region": target_name,
                            "Mean Chord Length (mm)": result["mean_mm"],
                            "Standard Deviation (mm)": result["std_mm"],
                            "Number of Point Pairs": N_SAMPLES,
                            "Bin Width (mm)": BIN_WIDTH_MM,
                        }
                    )

    # --------------------------------------------------------------
    # Save summary.
    # --------------------------------------------------------------

    summary_file = (
        output_dir
        / "csv"
        / "cld_summary.csv"
    )

    save_summary_csv(
        summary_rows,
        summary_file,
    )

    # Ask which of the generated CSV-backed CLDs should be included in
    # the multi-panel PDF. All CLDs have already been calculated/saved.
    create_selected_plot_pdf(
        cld_results,
        selected_source_keys,
        source_names,
        target_regions,
        output_dir,
    )

    print()
    print("=" * 70)
    print("CLD CALCULATION COMPLETE")
    print("=" * 70)
    print(f"CSV files: {output_dir / 'csv'}")
    print(f"Journal figures: {output_dir / 'figures'}")
    print(f"Summary:   {summary_file}")


if __name__ == "__main__":
    main()
