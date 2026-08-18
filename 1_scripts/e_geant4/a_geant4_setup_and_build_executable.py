# THis script sets up the Geant4 environment and builds the project using CMake and Make.

# Import necessary libraries
import os
import subprocess
from pathlib import Path
import json
import re
import b_config.a_config as config

# Configurations
executable = config.GEANT4_EXECUTABLE_FILE
geant4_include_dir = config.INCLUDE_DIR
fluence_source_types = config.FLUENCE_SOURCE_TYPES
energy_bins = config.ENERGY_BINS
energy_min = config.ENERGY_MIN
energy_max = config.ENERGY_MAX

def find_geant4make():
    """Search for geant4make.sh"""
    home = Path.home()

    for root, dirs, _ in os.walk(home):
        if "geant4make" in dirs:
            candidate = Path(root) / "geant4make" / "geant4make.sh"
            if candidate.is_file():
                return candidate

    raise FileNotFoundError("Could not find geant4make.sh")

def executable_exists():
    return executable.is_file()

def geant4_change_fluence_settings(params):

    photon_fluence_file = (geant4_include_dir / "TETPSPhotonFluence.hh")

    text = photon_fluence_file.read_text()

    enable_photon_fluence = (params["source_type"].lower()in fluence_source_types)

    text = re.sub(
        r"static constexpr G4bool ENABLE_PHOTON_FLUENCE\s*=\s*[^;]+;",
        f"static constexpr G4bool ENABLE_PHOTON_FLUENCE = "
        f"{str(enable_photon_fluence).lower()};",
        text
    )

    text = re.sub(
        r"static constexpr G4int nEnergyBins\s*=\s*[^;]+;",
        f"static constexpr G4int nEnergyBins = {energy_bins};",
        text
    )

    text = re.sub(
        r"static constexpr G4double Emin\s*=\s*[^;]+;",
        f"static constexpr G4double Emin = {energy_min};",
        text
    )

    text = re.sub(
        r"static constexpr G4double Emax\s*=\s*[^;]+;",
        f"static constexpr G4double Emax = {energy_max};",
        text
    )

    # Only write if the contents actually changed
    if text != photon_fluence_file.read_text():
        photon_fluence_file.write_text(text)

    return photon_fluence_file

def source_snapshot(project_dir):
    snapshot = {}

    for p in project_dir.rglob("*"):

        if (
            p.name == "CMakeLists.txt"
            or p.suffix in {".cc", ".hh", ".cpp", ".hpp"}
        ):
            snapshot[str(p.relative_to(project_dir))] = p.stat().st_mtime

    return snapshot

def build_geant4(params):

    geant4_change_fluence_settings(params)
    geant4make_path = find_geant4make()
    project_dir = config.INTERNAL_DIR.resolve()
    executable = project_dir / "build" / "Internal"
    build_dir = project_dir / "build"
    build_dir.mkdir(parents=True, exist_ok=True)

    snapshot_file = build_dir / ".cmake_snapshot.json"

    current = source_snapshot(project_dir)

    need_cmake = False
    if not snapshot_file.exists():
        print("First build.")
        need_cmake = True

    else:
        previous = json.loads(snapshot_file.read_text())

        if previous != current:
            print("\nSource tree changed (file edited, added, or removed).")
            need_cmake = True

    if need_cmake:

        cmd = (
            f'source "{geant4make_path}" && '
            f'cmake .. && '
            f'make'
        )

    else:

        print("\nNo source tree changes. Running make only.")

        cmd = (
            f'source "{geant4make_path}" && '
            f'make'
        )

    subprocess.run(
        ["bash", "-c", cmd],
        cwd=build_dir,
        check=True
    )

    # Save new snapshot after successful build
    snapshot_file.write_text(json.dumps(current, indent=2))

    executable = build_dir / "Internal"

    return executable