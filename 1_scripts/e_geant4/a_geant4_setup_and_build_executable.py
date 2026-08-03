# THis script sets up the Geant4 environment and builds the project using CMake and Make.

# temporary
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import necessary libraries
import os
import subprocess
from pathlib import Path
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
    return config.INTERNAL_FOLDER.is_file()

def build_geant4(project_dir: Path, geant4make_path: Path, rebuild: bool = False):

    executable = project_dir / "build" / "Internal"

    # Skip rebuilding only if the user didn't request it
    if executable.exists() and not rebuild:

        executable_time = executable.stat().st_mtime

        newest_source = max(
            p.stat().st_mtime
            for p in project_dir.rglob("*")
            if p.name == "CMakeLists.txt"
            or p.suffix in {".cc", ".hh", ".cpp", ".hpp"}
        )

        if executable_time >= newest_source:
            print("Executable is already up to date.")
            return executable

        print("Source files have changed. Rebuilding...")

    elif executable.exists() and rebuild:
        print("Rebuilding executable...")

    else:
        print("Executable not found. Building...")

    build_dir = project_dir / "build"
    build_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f'source "{geant4make_path}" && '
        f'cmake .. && '
        f'make'
    )

    subprocess.run(
        ["bash", "-c", cmd],
        cwd=build_dir,
        check=True
    )

    print(f"Build complete: {executable}")

    return executable

geant4make_path = find_geant4make()

project_dir = config.INTERNAL_FOLDER.resolve()

build_geant4(project_dir, geant4make_path)
print(project_dir)