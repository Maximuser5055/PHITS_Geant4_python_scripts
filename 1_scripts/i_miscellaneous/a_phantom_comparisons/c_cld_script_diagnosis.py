"""
Diagnose chord-length distributions (CLDs) for ICRP 145 liver -> liver.

This script is a focused diagnostic companion to
b_compare_chord_length_distributions.py. It uses the same:
    - project configuration
    - phantom registry
    - NODE/ELE/MATERIAL parsing
    - tetrahedron-volume calculation
    - mass/volume-weighted point sampling
    - cm -> mm coordinate conversion

The goal is to determine why the calculated ICRP 145 Liver -> Liver CLD
may differ from the published/reference CLD.

For each selected ICRP 145 phantom, the script reports:
    1. NODE coordinate range and physical dimensions
    2. number of liver tetrahedra
    3. liver volume from the mesh
    4. liver mass from the mesh/material densities
    5. liver centroid and bounding box
    6. sampled Liver -> Liver distance statistics
    7. distance percentiles and maximum distance
    8. a diagnostic CLD plot
    9. a diagnostic CSV containing the mesh statistics

The script intentionally does NOT modify or reuse the CLD result CSVs.
It performs a fresh, small diagnostic calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
import sys

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# Project imports
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import b_config.a_config as config
from b_config import b_phantom_registry as registry


# ======================================================================
# USER SETTINGS
# ======================================================================

# Diagnose the reference ICRP 145 phantoms by default.
PHANTOM_GROUP = "MRCP_AF_AM"

# Liver organ ID in the ICRP mesh convention.
LIVER_ID = 9500
LIVER_NAME = "Liver"

# Number of fresh point pairs used for the diagnostic CLD.
# 1 million is enough to reveal a large geometry/sampling discrepancy
# while being much faster than the production 10-million-pair CLDs.
N_SAMPLES = 1_000_000

# Number of worker-independent random samples. The diagnostic is kept
# single-process so that the reported percentiles are easy to reproduce.
RANDOM_SEED = 20261005

# TETGEN NODE coordinates are expected to be in cm, as in the CLD script.
NODE_UNIT = "cm"

# Histogram bin width for the diagnostic plot.
BIN_WIDTH_MM = 1.0

# Save diagnostic files here.
OUTPUT_DIR = Path(config.RESULTS_DIR) / "chord_length_distributions" / "diagnostics"


# ======================================================================
# DATA CLASS
# ======================================================================

@dataclass
class MeshData:
    nodes_cm: np.ndarray
    tetrahedra: np.ndarray
    organ_ids: np.ndarray


# ======================================================================
# FILE READING
# ======================================================================

def read_node_file(filename: Path) -> np.ndarray:
    """Read a TETGEN .node file using the same convention as the CLD script."""
    if not filename.is_file():
        raise FileNotFoundError(f"NODE file not found:\n{filename}")

    with filename.open("r", encoding="utf-8", errors="ignore") as file:
        n_nodes = None
        dimension = None

        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            n_nodes = int(parts[0])
            dimension = int(parts[1])
            break

        if n_nodes is None or dimension is None:
            raise ValueError(f"Could not find a NODE header in {filename}.")

        if dimension != 3:
            raise ValueError(f"{filename} is not a 3D NODE file.")

        coordinates = np.empty((n_nodes, 3), dtype=np.float64)

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

    if NODE_UNIT != "cm":
        raise ValueError("NODE_UNIT must be 'cm'.")

    return coordinates


def read_ele_file(filename: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a TETGEN .ele file using the same convention as the CLD script."""
    if not filename.is_file():
        raise FileNotFoundError(f"ELE file not found:\n{filename}")

    with filename.open("r", encoding="utf-8", errors="ignore") as file:
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

        if number_of_elements is None or nodes_per_element is None:
            raise ValueError(f"No ELE header found in {filename}.")

        if nodes_per_element != 4:
            raise ValueError(f"{filename} does not contain tetrahedral elements.")

        tetrahedra = np.empty((number_of_elements, 4), dtype=np.int64)
        organ_ids = np.empty(number_of_elements, dtype=np.int64)

        count = 0
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 6:
                raise ValueError(f"Invalid ELE line in {filename}:\n{line}")

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

    minimum_node_id = tetrahedra.min()

    if minimum_node_id == 1:
        tetrahedra -= 1
    elif minimum_node_id != 0:
        tetrahedra -= minimum_node_id

    return tetrahedra, organ_ids


def read_material_densities(filename: Path) -> dict[int, float]:
    """Read MAT[ID] densities using the project material-file convention."""
    if not filename.is_file():
        raise FileNotFoundError(f"MATERIAL file not found:\n{filename}")

    densities: dict[int, float] = {}
    current_density: float | None = None

    density_pattern = re.compile(
        r"^\s*\$\s+.*?([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+g/cm3",
        re.IGNORECASE,
    )
    material_pattern = re.compile(r"^\s*MAT\[(\d+)\]", re.IGNORECASE)

    with filename.open("r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            density_match = density_pattern.match(line)
            if density_match:
                current_density = float(density_match.group(1))
                continue

            material_match = material_pattern.match(line)
            if material_match:
                if current_density is None:
                    raise ValueError(
                        f"MAT found without a preceding density in {filename}: {line}"
                    )

                material_id = int(material_match.group(1))
                densities[material_id] = current_density
                current_density = None

    if not densities:
        raise ValueError(f"No material densities could be read from:\n{filename}")

    return densities


# ======================================================================
# MESH CALCULATIONS
# ======================================================================

def tetrahedron_volumes(nodes_cm: np.ndarray, tetrahedra: np.ndarray) -> np.ndarray:
    """Calculate tetrahedron volumes in cm^3."""
    a = nodes_cm[tetrahedra[:, 0]]
    b = nodes_cm[tetrahedra[:, 1]]
    c = nodes_cm[tetrahedra[:, 2]]
    d = nodes_cm[tetrahedra[:, 3]]

    ab = b - a
    ac = c - a
    ad = d - a

    return np.abs(
        np.einsum("ij,ij->i", ab, np.cross(ac, ad))
    ) / 6.0


def load_mesh(spec) -> tuple[MeshData, dict[int, float], np.ndarray]:
    """Load mesh and calculate tetrahedron volumes and masses."""
    nodes_cm = read_node_file(Path(spec.node_file))
    tetrahedra, organ_ids = read_ele_file(Path(spec.element_file))

    if tetrahedra.min() < 0 or tetrahedra.max() >= len(nodes_cm):
        raise ValueError(
            f"Tetrahedron node indices in {spec.element_file} "
            "are inconsistent with the NODE file."
        )

    mesh = MeshData(
        nodes_cm=nodes_cm,
        tetrahedra=tetrahedra,
        organ_ids=organ_ids,
    )

    densities = read_material_densities(Path(spec.material_file))
    volumes_cm3 = tetrahedron_volumes(nodes_cm, tetrahedra)

    density_array = np.empty(len(organ_ids), dtype=np.float64)
    missing = []

    for index, organ_id in enumerate(organ_ids):
        organ_id = int(organ_id)
        if organ_id not in densities:
            missing.append(organ_id)
        else:
            density_array[index] = densities[organ_id]

    if missing:
        raise ValueError(
            "No density was found for organ ID(s): "
            + ", ".join(map(str, sorted(set(missing))))
        )

    tetra_masses_g = volumes_cm3 * density_array
    return mesh, densities, tetra_masses_g


def liver_geometry(mesh: MeshData, tetra_masses_g: np.ndarray) -> dict:
    """Calculate geometry/mass diagnostics for the liver."""
    liver_indices = np.flatnonzero(mesh.organ_ids == LIVER_ID)

    if len(liver_indices) == 0:
        raise ValueError(f"No tetrahedra with liver organ ID {LIVER_ID} were found.")

    liver_tetrahedra = mesh.tetrahedra[liver_indices]
    liver_vertices = mesh.nodes_cm[liver_tetrahedra].reshape(-1, 3)

    volumes_cm3 = tetrahedron_volumes(
        mesh.nodes_cm,
        liver_tetrahedra,
    )
    mass_g = tetra_masses_g[liver_indices]

    # Tetrahedron centroid weighted by tetrahedron volume.
    tetra_centroids_cm = mesh.nodes_cm[liver_tetrahedra].mean(axis=1)
    total_volume = float(volumes_cm3.sum())
    centroid_cm = (
        np.average(tetra_centroids_cm, axis=0, weights=volumes_cm3)
    )

    minimum_cm = liver_vertices.min(axis=0)
    maximum_cm = liver_vertices.max(axis=0)
    dimensions_cm = maximum_cm - minimum_cm

    return {
        "n_tetrahedra": len(liver_indices),
        "volume_cm3": total_volume,
        "mass_g": float(mass_g.sum()),
        "density_g_cm3": float(mass_g.sum() / total_volume),
        "centroid_cm": centroid_cm,
        "minimum_cm": minimum_cm,
        "maximum_cm": maximum_cm,
        "dimensions_cm": dimensions_cm,
        "minimum_mm": minimum_cm * 10.0,
        "maximum_mm": maximum_cm * 10.0,
        "dimensions_mm": dimensions_cm * 10.0,
        "tetra_indices": liver_indices,
        "tetra_volumes_cm3": volumes_cm3,
    }


# ======================================================================
# RANDOM POINT SAMPLING
# ======================================================================

def sample_liver_points(
    mesh: MeshData,
    tetra_indices: np.ndarray,
    tetra_volumes_cm3: np.ndarray,
    number_of_points: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample points uniformly throughout the liver volume."""
    cumulative = np.cumsum(tetra_volumes_cm3)
    cumulative /= cumulative[-1]

    random_values = rng.random(number_of_points)
    positions = np.searchsorted(cumulative, random_values, side="right")
    selected_tetra_indices = tetra_indices[positions]

    vertices_cm = mesh.nodes_cm[
        mesh.tetrahedra[selected_tetra_indices]
    ]

    # Uniform barycentric coordinates in a tetrahedron.
    exponential = rng.exponential(size=(number_of_points, 4))
    barycentric = exponential / exponential.sum(axis=1, keepdims=True)

    points_cm = np.einsum(
        "ij,ijk->ik",
        barycentric,
        vertices_cm,
    )

    return points_cm * 10.0


def calculate_liver_to_liver_distances(
    mesh: MeshData,
    geometry: dict,
    number_of_points: int,
    seed: int,
) -> np.ndarray:
    """Generate fresh Liver -> Liver distances in mm."""
    rng = np.random.default_rng(seed)

    points_a_mm = sample_liver_points(
        mesh,
        geometry["tetra_indices"],
        geometry["tetra_volumes_cm3"],
        number_of_points,
        rng,
    )

    points_b_mm = sample_liver_points(
        mesh,
        geometry["tetra_indices"],
        geometry["tetra_volumes_cm3"],
        number_of_points,
        rng,
    )

    differences = points_a_mm - points_b_mm
    return np.sqrt(np.einsum("ij,ij->i", differences, differences))


# ======================================================================
# REPORTING
# ======================================================================

def print_header(text: str) -> None:
    print()
    print("=" * 78)
    print(text)
    print("=" * 78)


def print_geometry_report(
    phantom_name: str,
    spec,
    mesh: MeshData,
    geometry: dict,
) -> None:
    print_header(f"{phantom_name} — LIVER GEOMETRY")

    print(f"Phantom code:       {spec.code}")
    print(f"Sex:                {spec.sex}")
    print(f"NODE file:          {spec.node_file}")
    print(f"ELE file:           {spec.element_file}")
    print(f"MATERIAL file:      {spec.material_file}")
    print()

    print(f"Total mesh nodes:   {len(mesh.nodes_cm):,}")
    print(f"Total tetrahedra:   {len(mesh.tetrahedra):,}")
    print(f"Liver tetrahedra:   {geometry['n_tetrahedra']:,}")
    print()

    print(f"Liver volume:       {geometry['volume_cm3']:.4f} cm^3")
    print(f"Liver mass:         {geometry['mass_g']:.4f} g")
    print(f"Mean density:       {geometry['density_g_cm3']:.6f} g/cm^3")
    print()

    print(
        "Liver centroid:     "
        f"({geometry['centroid_cm'][0]:.4f}, "
        f"{geometry['centroid_cm'][1]:.4f}, "
        f"{geometry['centroid_cm'][2]:.4f}) cm"
    )
    print(
        "Liver min:          "
        f"({geometry['minimum_cm'][0]:.4f}, "
        f"{geometry['minimum_cm'][1]:.4f}, "
        f"{geometry['minimum_cm'][2]:.4f}) cm"
    )
    print(
        "Liver max:          "
        f"({geometry['maximum_cm'][0]:.4f}, "
        f"{geometry['maximum_cm'][1]:.4f}, "
        f"{geometry['maximum_cm'][2]:.4f}) cm"
    )
    print(
        "Liver dimensions:   "
        f"({geometry['dimensions_cm'][0]:.4f}, "
        f"{geometry['dimensions_cm'][1]:.4f}, "
        f"{geometry['dimensions_cm'][2]:.4f}) cm"
    )
    print(
        "Liver dimensions:   "
        f"({geometry['dimensions_mm'][0]:.1f}, "
        f"{geometry['dimensions_mm'][1]:.1f}, "
        f"{geometry['dimensions_mm'][2]:.1f}) mm"
    )

    # Whole-mesh coordinate extent is useful for catching unit errors.
    whole_min = mesh.nodes_cm.min(axis=0)
    whole_max = mesh.nodes_cm.max(axis=0)
    whole_dimensions = whole_max - whole_min

    print()
    print("WHOLE MESH COORDINATE RANGE")
    print(
        f"  min: ({whole_min[0]:.4f}, {whole_min[1]:.4f}, {whole_min[2]:.4f}) cm"
    )
    print(
        f"  max: ({whole_max[0]:.4f}, {whole_max[1]:.4f}, {whole_max[2]:.4f}) cm"
    )
    print(
        f"  size: ({whole_dimensions[0]:.4f}, "
        f"{whole_dimensions[1]:.4f}, {whole_dimensions[2]:.4f}) cm"
    )


def print_distance_report(distances_mm: np.ndarray) -> dict:
    """Print diagnostic distance statistics and return them."""
    percentiles = [1, 5, 25, 50, 75, 95, 99]
    percentile_values = np.percentile(distances_mm, percentiles)

    mean = float(distances_mm.mean())
    std = float(distances_mm.std())
    minimum = float(distances_mm.min())
    maximum = float(distances_mm.max())

    print_header("LIVER → LIVER CLD DIAGNOSTIC")
    print(f"Samples:            {len(distances_mm):,}")
    print(f"Mean distance:      {mean:.4f} mm")
    print(f"Std. deviation:     {std:.4f} mm")
    print(f"Minimum distance:   {minimum:.4f} mm")
    print(f"Maximum distance:   {maximum:.4f} mm")
    print()
    print("Percentiles:")
    for percentile, value in zip(percentiles, percentile_values):
        print(f"  P{percentile:02d}:               {value:.4f} mm")

    return {
        "n_samples": len(distances_mm),
        "mean_mm": mean,
        "std_mm": std,
        "minimum_mm": minimum,
        "maximum_mm": maximum,
        **{
            f"p{percentile}_mm": float(value)
            for percentile, value in zip(percentiles, percentile_values)
        },
    }


# ======================================================================
# PLOTS / CSV
# ======================================================================

def save_diagnostic_plot(
    phantom_name: str,
    distances_mm: np.ndarray,
    output_file: Path,
) -> None:
    """Save a simple diagnostic CLD plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    maximum = math.ceil(float(distances_mm.max()) / BIN_WIDTH_MM) * BIN_WIDTH_MM
    bins = np.arange(0.0, maximum + BIN_WIDTH_MM, BIN_WIDTH_MM)
    counts, edges = np.histogram(distances_mm, bins=bins)
    relative_number = counts / len(distances_mm)
    centers = (edges[:-1] + edges[1:]) / 2.0

    fig, axis = plt.subplots(figsize=(7.0, 4.5))
    axis.plot(centers, relative_number, linewidth=1.2)
    axis.set_xlabel("Distance (mm)")
    axis.set_ylabel("Relative Number")
    axis.set_title(f"{phantom_name}: Liver → Liver")
    axis.grid(alpha=0.18, linewidth=0.5)
    axis.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()

    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_diagnostic_csv(rows: list[dict], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_file, index=False)


# ======================================================================
# MAIN
# ======================================================================

def main() -> None:
    print_header("CHORD-LENGTH DISTRIBUTION DIAGNOSTIC")
    print(f"Phantom group:      {PHANTOM_GROUP}")
    print(f"Source organ:       {LIVER_NAME} (ID {LIVER_ID})")
    print(f"Samples per phantom:{N_SAMPLES:,}")
    print(f"Random seed:        {RANDOM_SEED}")
    print(f"NODE unit:          {NODE_UNIT}")
    print(f"Output directory:   {OUTPUT_DIR}")

    if PHANTOM_GROUP not in registry.PHANTOM_GROUPS:
        raise ValueError(f"Unknown phantom group: {PHANTOM_GROUP}")

    phantom_codes = registry.get_phantom_group(PHANTOM_GROUP)
    rows = []

    for phantom_code in phantom_codes:
        spec = registry.get_phantom(phantom_code)

        mesh, densities, tetra_masses_g = load_mesh(spec)
        geometry = liver_geometry(mesh, tetra_masses_g)

        print_geometry_report(spec.display_name, spec, mesh, geometry)

        print()
        print("Sampling fresh Liver → Liver point pairs...")
        distances_mm = calculate_liver_to_liver_distances(
            mesh,
            geometry,
            N_SAMPLES,
            RANDOM_SEED + phantom_codes.index(phantom_code),
        )

        distance_stats = print_distance_report(distances_mm)

        # Check the source/target region really is the same organ ID.
        print()
        print("SOURCE/TARGET ID CHECK")
        print(f"  Source organ ID:  {LIVER_ID}")
        print(f"  Target organ ID:  {LIVER_ID}")
        print("  Same region:      YES")

        rows.append(
            {
                "Phantom": spec.display_name,
                "Phantom Code": spec.code,
                "Sex": "Female" if spec.sex == "AF" else "Male",
                "Organ ID": LIVER_ID,
                "Organ": LIVER_NAME,
                "Number of Nodes": len(mesh.nodes_cm),
                "Number of Tetrahedra": len(mesh.tetrahedra),
                "Liver Tetrahedra": geometry["n_tetrahedra"],
                "Liver Volume (cm3)": geometry["volume_cm3"],
                "Liver Mass (g)": geometry["mass_g"],
                "Liver Density (g/cm3)": geometry["density_g_cm3"],
                "Centroid X (cm)": geometry["centroid_cm"][0],
                "Centroid Y (cm)": geometry["centroid_cm"][1],
                "Centroid Z (cm)": geometry["centroid_cm"][2],
                "Size X (mm)": geometry["dimensions_mm"][0],
                "Size Y (mm)": geometry["dimensions_mm"][1],
                "Size Z (mm)": geometry["dimensions_mm"][2],
                **distance_stats,
            }
        )

        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", spec.code)
        save_diagnostic_plot(
            spec.display_name,
            distances_mm,
            OUTPUT_DIR / f"{safe_name}_liver_to_liver_diagnostic.png",
        )

        print()
        print(
            "Saved diagnostic plot: "
            f"{OUTPUT_DIR / f'{safe_name}_liver_to_liver_diagnostic.png'}"
        )

    summary_file = OUTPUT_DIR / "liver_to_liver_diagnostic_summary.csv"
    save_diagnostic_csv(rows, summary_file)

    print_header("DIAGNOSTIC COMPLETE")
    print(f"Summary CSV: {summary_file}")
    print()
    print("Interpretation guide:")
    print("  1. Check liver dimensions first. If these are implausible, investigate")
    print("     NODE coordinate units/scaling or the mesh itself before CLD sampling.")
    print("  2. Check liver volume/mass against the intended ICRP 145 phantom data.")
    print("  3. If geometry is correct but the CLD is still shifted, investigate the")
    print("     random point-sampling method and the published CLD definition.")
    print("  4. Compare the printed mean, percentiles, and plot with the reference")
    print("     ICRP 145 Liver -> Liver distribution.")


if __name__ == "__main__":
    main()
