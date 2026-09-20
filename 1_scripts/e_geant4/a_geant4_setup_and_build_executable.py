# THis script sets up the Geant4 environment and builds the project using CMake and Make.

# Import necessary libraries
import os
import subprocess
from pathlib import Path
import json
import re
import b_config.a_config as config
from b_config.b_phantom_registry import PHANTOMS, PHANTOM_GROUPS, get_skeletal_ids

# Configurations
executable = config.GEANT4_EXECUTABLE_FILE
geant4_include_dir = config.INCLUDE_DIR
fluence_source_types = config.FLUENCE_SOURCE_TYPES
energy_bins = config.ENERGY_BINS
energy_min = config.ENERGY_MIN
energy_max = config.ENERGY_MAX
internal_dir = config.INTERNAL_DIR

tet_model_import_file = config.SRC_DIR / "TETModelImport.cc"
tet_run_action_file = config.SRC_DIR / "TETRunAction.cc"

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

def geant4_change_phantom_family(params):

    phantom_selection = params["phantom"]

    # Convert an individual phantom selection into its group
    if phantom_selection in PHANTOM_GROUPS:
        phantom_group = phantom_selection

    else:
        phantom_group = next(
            (
                group_code
                for group_code, group_spec in PHANTOM_GROUPS.items()
                if phantom_selection in group_spec.phantoms
            ),
            None
        )

        if phantom_group is None:
            raise ValueError(
                f"Phantom '{phantom_selection}' does not belong "
                f"to a registered phantom group."
            )

    male_phantom = next(
        phantom_spec.cell_file.stem
        for phantom_code in PHANTOM_GROUPS[phantom_group].phantoms
        if (phantom_spec := PHANTOMS[phantom_code]).sex == "AM"
    )

    female_phantom = next(
        phantom_spec.cell_file.stem
        for phantom_code in PHANTOM_GROUPS[phantom_group].phantoms
        if (phantom_spec := PHANTOMS[phantom_code]).sex == "AF"
    )

    text = tet_model_import_file.read_text()

    old_pattern = (
        r'// set phantom name\s*'
        r'if\(!isAF\)\s*phantomName\s*=\s*"[^"]+"\s*;\s*'
        r'else\s*phantomName\s*=\s*"[^"]+"\s*;'
    )

    new_code = (
        f'// set phantom name\n'
        f'if(!isAF) phantomName = "{male_phantom}";\n'
        f'else      phantomName = "{female_phantom}";'
    )

    new_text, count = re.subn(
        old_pattern,
        new_code,
        text
    )

    if count != 1:
        raise RuntimeError("Could not uniquely update the phantom selection in TETModelImport.cc.")

    if new_text != text:
        tet_model_import_file.write_text(new_text)

    return tet_model_import_file

def geant4_change_skeletal_ids(params):

    phantom_selection = params["phantom"]
    skeletal_ids = get_skeletal_ids(phantom_selection)

    text = tet_run_action_file.read_text()

    old_pattern = (
        r"// SKELETAL_IDS_BEGIN.*?"
        r"// SKELETAL_IDS_END"
    )

    new_code = (
        "// SKELETAL_IDS_BEGIN\n"
        "std::vector<G4int> skeletalIDs =\n"
        "{\n"
        + "".join(
            f"        {skeletal_id},\n"
            for skeletal_id in skeletal_ids
        )
        + "    };\n"
        "// SKELETAL_IDS_END"
    )

    new_text, count = re.subn(
        old_pattern,
        new_code,
        text,
        flags=re.DOTALL
    )

    if count != 1:
        raise RuntimeError(
            "Could not uniquely update the skeletal IDs "
            "in TETRunAction.cc."
        )

    if new_text != text:
        tet_run_action_file.write_text(new_text)

    return tet_run_action_file

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
    geant4_change_phantom_family(params)
    geant4_change_skeletal_ids(params)
    
    geant4make_path = find_geant4make()
    project_dir = internal_dir.resolve()
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
            print("\nSource tree changed (.hh and .cc files were edited, added, or removed).")
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