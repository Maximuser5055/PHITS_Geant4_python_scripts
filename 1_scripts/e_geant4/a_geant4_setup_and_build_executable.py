# THis script sets up the Geant4 environment and builds the project using CMake and Make.

# Import necessary libraries
import os
import subprocess
from pathlib import Path
import json
import b_config.a_config as config

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
    executable = config.INTERNAL_DIR / "build" / "Internal"
    return executable.is_file()

def source_snapshot(project_dir):
    snapshot = {}

    for p in project_dir.rglob("*"):

        if (
            p.name == "CMakeLists.txt"
            or p.suffix in {".cc", ".hh", ".cpp", ".hpp"}
        ):
            snapshot[str(p.relative_to(project_dir))] = p.stat().st_mtime

    return snapshot

def build_geant4():

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

        print("No source tree changes. Running make only.")

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