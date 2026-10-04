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
    - one vector PDF and one 1000-dpi RGB TIFF per source-organ/
      target-region combination, with all selected phantom groups
      combined into a single plot

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

# Figures are written as vector PDF and high-resolution TIFF, so use a
# non-interactive backend suitable for headless systems.
# This also avoids Qt-specific rendering problems on headless systems.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Journal figure typography.
# Elsevier recommends Arial (or Helvetica) for artwork.  This script uses
# Arial and embeds the TrueType font in the PDF when an Arial installation is
# available on the machine running the script.
JOURNAL_FONT = "Arial"

plt.rcParams.update({
    "font.family": JOURNAL_FONT,
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,

    # Embed TrueType fonts in vector PDF output.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
from matplotlib.lines import Line2D
from matplotlib import font_manager

def validate_journal_font() -> None:
    """Ensure the required journal font is installed before plotting."""
    try:
        font_manager.findfont(
            JOURNAL_FONT,
            fallback_to_default=False,
        )
    except ValueError as exc:
        raise RuntimeError(
            f"{JOURNAL_FONT!r} was not found. Install Arial on the "
            "machine running this script, or change JOURNAL_FONT to "
            "'Helvetica' if that is the font available in your "
            "environment."
        ) from exc


import b_config.a_config as config
from b_config import b_phantom_registry as registry
from c_database import b_organ_database as organ_database

# ======================================================================
# USER SETTINGS
# ======================================================================

# Registry group codes to compare.
# Phantom file paths are NOT hard-coded here.
PHANTOM_GROUPS_TO_COMPARE = [
    "MRCP_AF_AM_Filipino",
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
# Elsevier's general artwork guidance gives 140 mm as a typical 1.5-column
# width.  This is a good default for the four-phantom CLD plots while
# keeping the figure compact enough for a manuscript.
FIGURE_WIDTH_MM = 140.0
FIGURE_WIDTH = FIGURE_WIDTH_MM / 25.4
FIGURE_HEIGHT = 2.8

# Resolution for raster artwork.  Elsevier's artwork guidance specifies
# 1000 dpi for bitmapped line drawings.  The PDF output is vector and does
# not have a raster DPI limitation.
TIFF_DPI = 1000

# Prevent very tall multi-row figures from becoming enormous raster images.
# The row height is reduced automatically when there are many target regions.

# Ask the user which source organs to calculate.
ASK_SOURCE_SELECTION = True


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

def load_source_organs(source_csv: Path) -> list[int]:
    """
    Load source-organ IDs from SOURCE_CSV.

    The existing project uses the column:
        source_organ_ID
    """

    if not source_csv.is_file():
        raise FileNotFoundError(
            f"SOURCE_CSV was not found:\n{source_csv}"
        )

    df = clean_columns(pd.read_csv(source_csv))

    if "source_organ_ID" not in df.columns:
        raise ValueError(
            f"Column 'source_organ_ID' was not found in:\n{source_csv}\n"
            f"Available columns: {list(df.columns)}"
        )

    source_ids = []

    for value in df["source_organ_ID"].dropna():
        try:
            source_id = int(float(str(value).strip()))
        except ValueError as exc:
            raise ValueError(
                f"Invalid source organ ID '{value}' in {source_csv}."
            ) from exc

        if source_id not in source_ids:
            source_ids.append(source_id)

    if not source_ids:
        raise ValueError(
            f"No source-organ IDs were found in {source_csv}."
        )

    return source_ids


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


# ======================================================================
# SEX-SPECIFIC TARGET REGIONS
# ======================================================================

# Organ IDs that do not exist in the opposite-sex phantom.
# These are skipped rather than treated as errors when building
# target-region samplers.
SEX_EXCLUDED_ORGAN_IDS = {
    "AF": {
        11500,       # Prostate
        12900, 13000, # Testes
    },
    "AM": {
        11100, 11200, # Ovaries
        13900,        # Uterus/cervix
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

    if NODE_UNIT == "mm":
        return coordinates / 10.0

    if NODE_UNIT == "m":
        return coordinates * 100.0

    raise ValueError(
        "NODE_UNIT must be 'mm', 'cm', or 'm'."
    )


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
    source_id: int,
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

    df = pd.DataFrame(
        {
            "Phantom": phantom.display_name,
            "Phantom Code": phantom.code,
            "Sex": (
                "Male"
                if phantom.sex == "AM"
                else "Female"
            ),
            "Source Organ ID": source_id,
            "Target Region": target_region,
            "Distance (mm)": distances_mm,
            "Relative Number": relative_number,
            "Mean Chord Length (mm)": mean_mm,
            "Standard Deviation (mm)": std_mm,
        }
    )

    # The phantom/source/target metadata and the summary statistics are
    # constant for the entire CLD.  Keep them only on the first row so
    # the CSV is easier to read without changing the distance distribution.
    repeated_columns = [
        "Phantom",
        "Phantom Code",
        "Sex",
        "Source Organ ID",
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
                "Source Organ ID",
                "Target Region",
                "Phantom",
            ]
        ).reset_index(drop=True)

    df.to_csv(
        output_file,
        index=False,
    )


# ======================================================================
# PLOTTING
# ======================================================================

def plot_source_clds(
    source_id: int,
    source_name: str,
    target_region: str,
    cld_results: dict,
    output_file: Path,
):
    """
    Create one journal-ready figure for one source-organ/target-region combination.

    All selected phantom groups are combined into one plot.

    The legend order follows the `phantoms=()` tuples in
    registry.PHANTOM_GROUPS exactly.  The phantom group display names are
    not used as legend entries.

    Colors and line styles are assigned automatically from Matplotlib's
    configured cycles, so adding phantoms does not require
    sex-specific plotting logic.
    """

    # Build the plotting order directly from the phantom tuples in the
    # registry.  This makes `PHANTOM_GROUPS[group_code].phantoms` the
    # authoritative source for the legend order.
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

    if not ordered_phantoms:
        return

    figure, axis = plt.subplots(
        1,
        1,
        figsize=(
            FIGURE_WIDTH,
            FIGURE_HEIGHT,
        ),
    )

    # Use Matplotlib's configured color cycle instead of hard-coded colors.
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color",
        [],
    )

    if not color_cycle:
        color_cycle = ["C0"]

    # Use Matplotlib's available line-style definitions instead of
    # hard-coding styles for male/female.
    line_styles = [
        linestyle
        for linestyle in Line2D.lineStyles
        if linestyle != "None"
    ]

    if not line_styles:
        line_styles = ["-"]

    # Combine colors and line styles so additional phantoms can be added
    # without requiring new sex-specific style logic.
    style_cycle = [
        (color, linestyle)
        for color in color_cycle
        for linestyle in line_styles
    ]

    for index, (group_name, phantom_code, phantom) in enumerate(
        ordered_phantoms
    ):
        entries = cld_results.get(
            group_name,
            {},
        ).get(
            target_region,
            {},
        )

        # Results are keyed by phantom code, matching the registry tuple.
        result = entries.get(phantom_code)

        if result is None:
            continue

        color, linestyle = style_cycle[
            index % len(style_cycle)
        ]

        label = (
            f"{phantom.display_name}"
            f" ({result['mean_mm']:.2f} "
            f"± {result['std_mm']:.2f} mm)"
        )

        axis.plot(
            result["distance_mm"],
            result["relative_number"],
            color=color,
            linestyle=linestyle,
            linewidth=1.2,
            label=label,
        )

    # Put the target-region name inside the plot rather than in the
    # figure title.
    axis.text(
        0.02,
        0.96,
        target_region,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )

    axis.set_xlabel(
        "Distance (mm)"
    )

    axis.set_ylabel(
        "Relative Number"
    )

    axis.grid(
        alpha=0.18,
        linewidth=0.5,
    )

    axis.legend(
        fontsize=8,
        frameon=False,
        loc="best",
    )

    figure.tight_layout()

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # The PDF is the preferred manuscript/vector output: lines, axes,
    # labels, and legend remain resolution-independent in LaTeX.
    pdf_file = output_file.with_suffix(".pdf")

    # The TIFF is the raster submission alternative.  These CLD graphs are
    # line-art figures, so use 1000 dpi rather than the 300 dpi used for
    # photographs/halftones.  RGB is retained for the colored curves.
    tiff_file = output_file.with_suffix(".tiff")

    figure.savefig(
        pdf_file,
        format="pdf",
        bbox_inches=None,
    )

    figure.savefig(
        tiff_file,
        dpi=TIFF_DPI,
        format="tiff",
        bbox_inches=None,
        facecolor="white",
    )

    plt.close(figure)


# ======================================================================
# USER INTERFACE
# ======================================================================

def select_source_organs(
    source_ids: list[int],
    source_names: dict[int, str],
) -> list[int]:

    if not ASK_SOURCE_SELECTION:
        return source_ids

    print()
    print("=" * 70)
    print("SOURCE ORGAN SELECTION")
    print("=" * 70)

    for index, source_id in enumerate(
        source_ids,
        start=1,
    ):
        name = source_names.get(
            source_id,
            "Unknown",
        )

        print(
            f"[{index}] {source_id} ({name})"
        )

    print("[A] All source organs")

    while True:

        choice = input(
            "\nSelect source organ(s): "
        ).strip().upper()

        if choice == "A":
            return source_ids

        try:
            indices = [
                int(value.strip()) - 1
                for value in choice.split(",")
            ]

            if not indices:
                raise ValueError

            if any(
                index < 0
                or index >= len(source_ids)
                for index in indices
            ):
                raise ValueError

            selected = list(
                dict.fromkeys(
                    source_ids[index]
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

    validate_journal_font()

    print()
    print("=" * 70)
    print("CHORD-LENGTH DISTRIBUTION CALCULATION")
    print("=" * 70)

    source_csv = Path(config.SOURCE_CSV)

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

    source_ids = load_source_organs(
        source_csv
    )

    # Resolve the target-region file through the phantom registry.
    target_region_file = registry.get_target_region_file(
        PHANTOM_GROUPS_TO_COMPARE[0]
    )

    # All comparison groups must use the same target-region mapping.
    for group_code in PHANTOM_GROUPS_TO_COMPARE:
        group_target_file = registry.get_target_region_file(
            group_code
        )

        if Path(group_target_file) != Path(target_region_file):
            raise ValueError(
                "The selected phantom groups use different "
                "target-region files. The current CLD calculation "
                "expects one common source-target mapping."
            )

    target_regions = load_target_regions(
        Path(target_region_file)
    )

    # Source-organ names are taken from c_database.b_organ_database.py.
    # SOURCE_CSV remains the authoritative source of the selected organ IDs.
    source_names: dict[int, str] = {}

    for phantom_organs in organ_database.ORGANS.values():
        for organ_id, organ_data in phantom_organs.items():
            source_names.setdefault(
                int(organ_id),
                str(organ_data["name"]),
            )

    selected_source_ids = select_source_organs(
        source_ids,
        source_names,
    )

    # --------------------------------------------------------------
    # Resolve phantom groups
    # --------------------------------------------------------------

    resolved_groups = resolve_phantoms(registry)

    # Both phantom groups use the user-supplied target-region mapping.
    # This is intentional: SOURCE_CSV and TARGET_REGION_FILE are the
    # only CSV inputs required by the CLD script.
    validate_phantom_files(
        resolved_groups
    )

    # --------------------------------------------------------------
    # Load all meshes once
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

            phantom_target_regions = filter_target_regions_for_sex(
                target_regions,
                phantom.sex,
            )

            region_samplers = build_region_samplers(
                mesh,
                tetra_masses,
                phantom_target_regions,
            )

            meshes[phantom.code] = mesh
            samplers[phantom.code] = region_samplers

            # Pre-build samplers for every source organ so that
            # density/volume processing is not repeated for every
            # target region.
            source_samplers[phantom.code] = {}
            for source_id in source_ids:
                if source_id not in mesh.organ_ids:
                    continue

                source_samplers[phantom.code][source_id] = (
                    build_region_samplers(
                        mesh,
                        tetra_masses,
                        [
                            {
                                "name": f"Source_{source_id}",
                                "ids": (source_id,),
                            }
                        ],
                    )[f"Source_{source_id}"]
                )

            print(
                f"    Nodes: {len(mesh.nodes_cm):,}"
            )
            print(
                f"    Tetrahedra: {len(mesh.tetrahedra):,}"
            )

    # --------------------------------------------------------------
    # Calculate CLDs
    # --------------------------------------------------------------

    summary_rows = []

    # Results structure:
    #
    #   group_name
    #       target_region
    #           phantom_code
    #
    # The phantom code is used as the key so the plotting order can follow
    # the `phantoms=()` tuple in registry.PHANTOM_GROUPS exactly.
    cld_results = {
        group_name: {}
        for group_name in resolved_groups
    }

    for source_id in selected_source_ids:

        source_name = source_names.get(
            source_id,
            f"Organ {source_id}",
        )

        print()
        print("=" * 70)
        print(
            f"SOURCE: {source_name} "
            f"(ID {source_id})"
        )
        print("=" * 70)

        for group_name, phantoms in resolved_groups.items():

            print()
            print(
                f"Phantom group: {group_name}"
            )

            cld_results[group_name] = {}

            for target_region in target_regions:

                target_name = target_region["name"]

                cld_results[group_name][
                    target_name
                ] = {}

                print(
                    f"  Target: {target_name}"
                )

                for phantom in phantoms:

                    mesh = meshes[
                        phantom.code
                    ]

                    phantom_samplers = samplers[
                        phantom.code
                    ]

                    # SOURCE_CSV contains individual source organ IDs.
                    if source_id not in mesh.organ_ids:
                        print(
                            f"    [SKIP] {phantom.display_name}: "
                            f"source organ ID {source_id} "
                            "not present in mesh."
                        )
                        continue

                    source_sampler = source_samplers[
                        phantom.code
                    ][source_id]

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
                        f"{source_id}|"
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

                    cld_results[group_name][
                        target_name
                    ][phantom.code] = {
                        "phantom": phantom,
                        "distance_mm": distances_mm,
                        "relative_number": relative_number,
                        "mean_mm": mean_mm,
                        "std_mm": std_mm,
                    }

                    # --------------------------------------------------
                    # Save CLD CSV
                    # --------------------------------------------------

                    csv_name = (
                        f"{safe_filename(phantom.code)}_"
                        f"source_{source_id}_"
                        f"target_{safe_filename(target_name)}_"
                        f"cld.csv"
                    )

                    save_cld_csv(
                        output_dir / "csv" / csv_name,
                        phantom,
                        source_id,
                        target_name,
                        distances_mm,
                        relative_number,
                        mean_mm,
                        std_mm,
                    )

                    summary_rows.append(
                        {
                            "Phantom": phantom.display_name,
                            "Phantom Code": phantom.code,
                            "Sex": (
                                "Male"
                                if phantom.sex == "AM"
                                else "Female"
                            ),
                            "Source Organ ID": source_id,
                            "Source Organ": source_name,
                            "Target Region": target_name,
                            "Mean Chord Length (mm)": mean_mm,
                            "Standard Deviation (mm)": std_mm,
                            "Number of Point Pairs": N_SAMPLES,
                            "Bin Width (mm)": BIN_WIDTH_MM,
                        }
                    )

                    print(
                        f"    {phantom.display_name}: "
                        f"{mean_mm:.2f} ± {std_mm:.2f} mm"
                    )

        # --------------------------------------------------------------
        # Save one plot per source-organ/target-region combination.
        # --------------------------------------------------------------

        for region in target_regions:
            target_name = region["name"]

            figure_name = (
                f"source_{source_id}_"
                f"target_{safe_filename(target_name)}.png"
            )

            # The .png suffix here is only used as a temporary base name;
            # plot_source_clds replaces it with .pdf and .tiff.
            figure_base = output_dir / "figures" / figure_name

            plot_source_clds(
                source_id,
                source_name,
                target_name,
                cld_results,
                figure_base,
            )

            print(
                f"\nSaved journal figures:\n"
                f"  {figure_base.with_suffix('.pdf')}\n"
                f"  {figure_base.with_suffix('.tiff')}"
            )

    # --------------------------------------------------------------
    # Save summary
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

    print()
    print("=" * 70)
    print("CLD CALCULATION COMPLETE")
    print("=" * 70)
    print(f"CSV files: {output_dir / 'csv'}")
    print(f"Journal figures: {output_dir / 'figures'}")
    print(f"Summary:   {summary_file}")


if __name__ == "__main__":
    main()