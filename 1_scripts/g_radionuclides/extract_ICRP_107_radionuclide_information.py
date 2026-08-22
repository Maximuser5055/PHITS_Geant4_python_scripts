"""
ICRP-107 RAD/BET extractor and viewer.

Extracts:
    ICRP-107.RAD
        - nuclide name
        - physical half-life
        - half-life units
        - number of data records
        - ICODE
        - radiation type
        - yield
        - energy (MeV)
        - JCODE

    ICRP-107.BET
        - nuclide name
        - number of data records
        - energy (MeV)
        - yield

The two files are combined into one CSV.  A summary is printed and saved
to a text file.  The user can then select radionuclides to inspect and
optionally generate a log-log emission plot.

Important:
    The ICRP-107 files are fixed-width files.  Do NOT parse them only with
    str.split(), because some BET records contain no whitespace between
    the F7.0 and E10.0 fields.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


# ---------------------------------------------------------------------------
# ICODE / JCODE definitions from ICRP-107 Table 3
# ---------------------------------------------------------------------------

RADIATION_TYPES = {
    "1": "Gamma rays",
    "2": "X rays",
    "3": "Annihilation photons",
    "4": "Beta+ particles",
    "5": "Beta- particles",
    "6": "IC electrons",
    "7": "Auger electrons",
    "8": "Alpha particles",
    "9": "Alpha recoil nuclei",
    "10": "Fission fragments",
    "11": "Neutrons",
}

JCODE_TYPES = {
    "G": "Gamma rays",
    "PG": "Prompt gamma rays",
    "DG": "Delayed gamma rays",
    "X": "X rays",
    "AQ": "Annihilation photons",
    "B+": "Beta+ particles",
    "B-": "Beta- particles",
    "BD": "Delayed beta particles",
    "IE": "IC electrons",
    "AE": "Auger electrons",
    "A": "Alpha particles",
    "AR": "Alpha recoil nuclei",
    "FF": "Fission fragments",
    "N": "Neutrons",
}

CSV_FIELDS = [
    "source_file (1=RAD, 2=BET)",
    "nuclide",
    "physical_half_life",
    "half_life_units",
    "number_of_data_records",
    "icode",
    "jcode",
    "energy_mev",
    "yield (/nt for RAD file, beta/MeV/nt for BET file)",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean_float(value: str) -> float:
    """Convert an ICRP scientific-notation field to float."""
    value = value.strip().replace("D", "E").replace("d", "e")
    return float(value)


def safe_filename(name: str) -> str:
    """Make a radionuclide name safe for use as a filename."""
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)


def prompt_path(prompt: str, default: Path) -> Path:
    text = input(f"{prompt} [{default}]: ").strip()
    return Path(text).expanduser() if text else default


# ---------------------------------------------------------------------------
# RAD parser
# ---------------------------------------------------------------------------

def parse_rad(rad_path: Path) -> tuple[list[dict], dict]:
    """
    Parse ICRP-107.RAD.

    Fixed-width layout:
        Nuclide      A7
        T1/2         E11.0
        Time unit    A2
        N            I9

        ICODE        A2
        Yield        E12.0
        Energy       E12.0
        JCODE        A3
    """
    records = []
    radionuclides = []
    total_declared_records = 0
    radiation_counter = Counter()
    half_life_unit_counter = Counter()

    with rad_path.open("r", encoding="ascii", errors="replace") as f:
        line_number = 0

        while True:
            header = f.readline()
            if not header:
                break
            line_number += 1

            header = header.rstrip("\r\n")
            if not header.strip():
                continue

            if len(header) < 29:
                raise ValueError(
                    f"Malformed RAD header at line {line_number}: {header!r}"
                )

            # A7 + E11.0 + A2 + I9 = 29 characters
            nuclide = header[0:7].strip()
            half_life_text = header[7:18].strip()
            half_life_unit = header[18:20].strip()
            n_text = header[20:29].strip()

            try:
                number_of_records = int(n_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid RAD record count at line {line_number}: {n_text!r}"
                ) from exc

            physical_half_life = clean_float(half_life_text)

            radionuclides.append(nuclide)
            total_declared_records += number_of_records
            half_life_unit_counter[half_life_unit] += 1

            for _ in range(number_of_records):
                data = f.readline()
                line_number += 1

                if not data:
                    raise ValueError(
                        f"Unexpected end of RAD file while reading {nuclide}."
                    )

                data = data.rstrip("\r\n")

                if len(data) < 29:
                    raise ValueError(
                        f"Malformed RAD data at line {line_number}: {data!r}"
                    )

                # A2 + E12.0 + E12.0 + A3 = 29 characters
                icode = data[0:2].strip()
                yield_value = clean_float(data[2:14])
                energy_mev = clean_float(data[14:26])
                jcode = data[26:29].strip()

                radiation_type = RADIATION_TYPES.get(
                    icode, f"Unknown ICODE {icode}"
                )

                radiation_counter[radiation_type] += 1

                records.append(
                    {
                        "source_file (1=RAD, 2=BET)": 1,
                        "nuclide": nuclide,
                        "physical_half_life": physical_half_life,
                        "half_life_units": half_life_unit,
                        "number_of_data_records": number_of_records,
                        "icode": icode,
                        "radiation_type": radiation_type,
                        "jcode": jcode,
                        "energy_mev": energy_mev,
                        "yield (/nt for RAD file, beta/MeV/nt for BET file)": yield_value
                    }
                )

    summary = {
        "radionuclides": radionuclides,
        "total_declared_records": total_declared_records,
        "radiation_counter": radiation_counter,
        "half_life_unit_counter": half_life_unit_counter,
    }

    return records, summary


# ---------------------------------------------------------------------------
# BET parser
# ---------------------------------------------------------------------------

def parse_bet(bet_path: Path) -> tuple[list[dict], dict]:
    """
    Parse ICRP-107.BET.

    Fixed-width layout:
        Nuclide      A7
        N            I10

        Energy       F7.0
        Number       E10.0

    The fixed-width parsing is intentional.  For example:
        10.000001.711E-03
    represents:
        energy = 10.00000 MeV
        yield  = 1.711E-03
    """
    records = []
    radionuclides = []
    total_declared_records = 0

    with bet_path.open("r", encoding="ascii", errors="replace") as f:
        line_number = 0

        while True:
            header = f.readline()
            if not header:
                break
            line_number += 1

            header = header.rstrip("\r\n")
            if not header.strip():
                continue

            if len(header) < 17:
                raise ValueError(
                    f"Malformed BET header at line {line_number}: {header!r}"
                )

            # A7 + I10 = 17 characters
            nuclide = header[0:7].strip()
            n_text = header[7:17].strip()

            try:
                number_of_records = int(n_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid BET record count at line {line_number}: {n_text!r}"
                ) from exc

            radionuclides.append(nuclide)
            total_declared_records += number_of_records

            for _ in range(number_of_records):
                data = f.readline()
                line_number += 1

                if not data:
                    raise ValueError(
                        f"Unexpected end of BET file while reading {nuclide}."
                    )

                data = data.rstrip("\r\n")

                # F7.0 + E10.0 = 17 characters
                if len(data) < 17:
                    raise ValueError(
                        f"Malformed BET data at line {line_number}: {data!r}"
                    )

                energy_mev = clean_float(data[0:7])
                yield_value = clean_float(data[7:17])

                records.append(
                    {
                        "source_file (1=RAD, 2=BET)": 2,
                        "nuclide": nuclide,
                        "physical_half_life": "",
                        "half_life_units": "",
                        "number_of_data_records": number_of_records,
                        "icode": "5",
                        "radiation_type": "Beta- particles",
                        "jcode": "",
                        "energy_mev": energy_mev,
                        "yield (/nt for RAD file, beta/MeV/nt for BET file)": yield_value
                    }
                )

    summary = {
        "radionuclides": radionuclides,
        "total_declared_records": total_declared_records,
    }

    return records, summary


# ---------------------------------------------------------------------------
# CSV + summary
# ---------------------------------------------------------------------------

def write_csv(records: list[dict], output_csv: Path) -> None:
    with output_csv.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def build_summary(
    rad_summary: dict,
    bet_summary: dict,
    records: list[dict],
) -> str:
    rad_set = set(rad_summary["radionuclides"])
    bet_set = set(bet_summary["radionuclides"])
    combined_set = rad_set | bet_set

    source_counter = Counter(r["source_file (1=RAD, 2=BET)"] for r in records)

    radiation_counter = Counter(
        r["radiation_type"]
        for r in records
        if r["source_file (1=RAD, 2=BET)"] == 1
    )

    lines = [
        "=" * 78,
        "ICRP-107 RAD/BET EXTRACTION SUMMARY",
        "=" * 78,
        "",
        "INPUT FILES",
        "-" * 78,
        f"RAD file radionuclides:       {len(rad_set):,}",
        f"BET file radionuclides:       {len(bet_set):,}",
        f"Combined unique radionuclides: {len(combined_set):,}",
        f"Radionuclides in both files:  {len(rad_set & bet_set):,}",
        f"RAD-only radionuclides:       {len(rad_set - bet_set):,}",
        f"BET-only radionuclides:       {len(bet_set - rad_set):,}",
        "",
        "DATA RECORDS",
        "-" * 78,
        f"RAD declared data records:    {rad_summary['total_declared_records']:,}",
        f"BET declared data records:    {bet_summary['total_declared_records']:,}",
        f"Combined CSV rows:             {len(records):,}",
        f"RAD rows in CSV:               {source_counter.get(1, 0):,}",
        f"BET rows in CSV:               {source_counter.get(2, 0):,}",
        "",
        "RAD RADIATION TYPES",
        "-" * 78,
    ]

    for radiation_type, count in radiation_counter.most_common():
        lines.append(f"{radiation_type:<30} {count:>12,}")

    lines.extend(
        [
            "",
            "RAD HALF-LIFE UNITS",
            "-" * 78,
        ]
    )

    for unit, count in sorted(
        rad_summary["half_life_unit_counter"].items()
    ):
        lines.append(f"{unit:<30} {count:>12,}")

    lines.extend(
        [
            "",
            "ICODE / JCODE DEFINITIONS",
            "-" * 90,
            "1  G   Gamma rays",
            "1  PG  Prompt gamma rays",
            "1  DG  Delayed gamma rays",
            "2  X   X rays",
            "3  AQ  Annihilation photons",
            "4  B+  Beta+ particles",
            "5  B-  Beta- particles",
            "5  BD  Delayed beta particles",
            "6  IE  IC electrons",
            "7  AE  Auger electrons",
            "8  A   Alpha particles",
            "9  AR  Alpha recoil nuclei",
            "10 FF  Fission fragments",
            "11 N   Neutrons",
            "",
            "=" * 90,
        ]
    )

    return "\n".join(lines)


def print_summary(summary_text: str) -> None:
    print("\n" + summary_text + "\n")


# ---------------------------------------------------------------------------
# Individual radionuclide viewer
# ---------------------------------------------------------------------------

def get_records_for_nuclide(
    records: list[dict], nuclide: str
) -> list[dict]:
    return [
        r for r in records
        if r["nuclide"].upper() == nuclide.upper()
    ]


def print_radionuclide(
    records: list[dict], nuclide: str
) -> None:
    selected = get_records_for_nuclide(records, nuclide)

    if not selected:
        print(f"\nNo records found for {nuclide}.")
        return

    rad = [r for r in selected if r["source_file (1=RAD, 2=BET)"] == 1]
    bet = [r for r in selected if r["source_file (1=RAD, 2=BET)"] == 2]

    print("\n" + "=" * 100)
    print(f"RADIONUCLIDE: {nuclide}")
    print("=" * 100)

    if rad:
        first = rad[0]
        print("\n[RAD HEADER]")
        print(f"  Nuclide:               {first['nuclide']}")
        print(f"  Physical half-life:    {first['physical_half_life']}")
        print(f"  Half-life units:       {first['half_life_units']}")
        print(f"  Number of data records:{first['number_of_data_records']:,}")

        print("\n[RAD DATA]")
        print(
            f"{'ICODE':>6} {'Radiation type':<25} {'JCODE':<5} "
            f"{'Yield':>15} {'Energy (MeV)':>15}"
        )
        print("-" * 100)

        for r in rad:
            print(
                f"{r['icode']:>6} "
                f"{r['radiation_type']:<25} "
                f"{r['jcode']:<5} "
                f"{r['yield']:>15.6E} "
                f"{r['energy_mev']:>15.6E}"
            )
    else:
        print("\nNo RAD records for this radionuclide.")

    if bet:
        first = bet[0]
        print("\n[BET HEADER]")
        print(f"  Nuclide:               {first['nuclide']}")
        print(f"  Number of data records:{first['number_of_data_records']:,}")

        print("\n[BET DATA]")
        print(
            f"{'Energy (MeV)':>15} "
            f"{'Yield / MeV / transformation':>35}"
        )
        print("-" * 55)

        for r in bet:
            print(
                f"{r['energy_mev']:>15.6E} "
                f"{r['yield']:>35.6E}"
            )
    else:
        print("\nNo BET records for this radionuclide.")

    print("=" * 100)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_radionuclide(
    records: list[dict],
    nuclide: str,
    output_directory: Path,
) -> Path | None:
    if plt is None:
        print(
            "\nMatplotlib is not installed. "
            "Install it with: pip install matplotlib"
        )
        return None

    selected = get_records_for_nuclide(records, nuclide)

    if not selected:
        print(f"No records found for {nuclide}.")
        return None

    # Only positive x/y values can appear on a log-log plot.
    positive = [
        r for r in selected
        if r["energy_mev"] > 0 and r["yield"] > 0
    ]

    if not positive:
        print(f"No positive energy/yield data available for {nuclide}.")
        return None

    rad = [r for r in positive if r["source_file (1=RAD, 2=BET)"] == 1]
    bet = [r for r in positive if r["source_file (1=RAD, 2=BET)"] == 2]

    fig, ax = plt.subplots(figsize=(10, 7))

    # RAD: discrete emission lines
    rad_by_type = defaultdict(list)
    for r in rad:
        rad_by_type[r["radiation_type"]].append(r)

    for radiation_type, rows in sorted(rad_by_type.items()):
        rows = sorted(rows, key=lambda x: x["energy_mev"])
        ax.scatter(
            [r["energy_mev"] for r in rows],
            [r["yield"] for r in rows],
            label=radiation_type,
            s=28,
            alpha=0.75,
        )

    # BET: continuous beta spectrum
    if bet:
        bet_rows = sorted(bet, key=lambda x: x["energy_mev"])
        ax.plot(
            [r["energy_mev"] for r in bet_rows],
            [r["yield"] for r in bet_rows],
            label="Beta- particles",
            linewidth=1.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Emission yield")

    ax.set_title(f"ICRP-107 Emission Spectrum: {nuclide}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output_directory.mkdir(parents=True, exist_ok=True)
    graph_path = output_directory / (
        f"{safe_filename(nuclide)}_emission_loglog.png"
    )

    fig.savefig(graph_path, dpi=300, bbox_inches="tight")
    print(f"\nGraph saved to:\n  {graph_path}")

    try:
        plt.show()
    except Exception as exc:
        print(f"Could not display the graph interactively: {exc}")
        print("The PNG file was still saved successfully.")

    plt.close(fig)
    return graph_path


# ---------------------------------------------------------------------------
# Interactive selection
# ---------------------------------------------------------------------------

def interactive_viewer(
    records: list[dict],
    output_directory: Path,
) -> None:
    answer = input(
        "\nDo you want to view individual radionuclides? [y/N]: "
    ).strip().lower()

    if answer not in {"y", "yes"}:
        print("\nIndividual radionuclide viewing skipped.")
        return

    available = sorted(
        {r["nuclide"] for r in records},
        key=lambda x: x.upper(),
    )

    print(
        "\nEnter radionuclide names separated by commas "
        "(e.g. Ac-223, Na-22)."
    )
    print("Enter 'all' to view every radionuclide.")
    print("Enter 'list' to print the available radionuclides.")

    while True:
        selection = input("\nRadionuclide selection: ").strip()

        if selection.lower() == "list":
            print("\nAvailable radionuclides:")
            print(", ".join(available))
            continue

        if selection.lower() == "all":
            print(
                f"\nWARNING: This will print {len(available):,} "
                "radionuclides and may produce a very large terminal output."
            )
            confirm = input("Continue? [y/N]: ").strip().lower()
            if confirm not in {"y", "yes"}:
                continue
            selected_names = available
        else:
            selected_names = [
                item.strip()
                for item in selection.split(",")
                if item.strip()
            ]

        if not selected_names:
            print("No radionuclides selected.")
            continue

        # Preserve user order while matching case-insensitively.
        lookup = {name.upper(): name for name in available}
        normalized = []

        for name in selected_names:
            key = name.upper()
            if key not in lookup:
                print(f"Radionuclide not found: {name}")
            else:
                normalized.append(lookup[key])

        if not normalized:
            print("No valid radionuclides selected.")
            continue

        for name in normalized:
            print_radionuclide(records, name)

            graph_answer = input(
                f"\nDo you want a log-log emission graph for {name}? [y/N]: "
            ).strip().lower()

            if graph_answer in {"y", "yes"}:
                plot_radionuclide(records, name, output_directory)

        another = input(
            "\nView another radionuclide selection? [y/N]: "
        ).strip().lower()

        if another not in {"y", "yes"}:
            break


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 90)
    print("ICRP-107 RAD/BET EXTRACTOR")
    print("=" * 90)

    IRCP_107_dir = Path("/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/ICRP_107")

    default_rad = IRCP_107_dir / "ICRP-107.RAD"
    default_bet = IRCP_107_dir / "ICRP-107.BET"

    rad_path = prompt_path("Path to .RAD file", default_rad)
    bet_path = prompt_path("Path to .BET file", default_bet)

    if not rad_path.exists():
        raise FileNotFoundError(f"RAD file not found: {rad_path}")

    if not bet_path.exists():
        raise FileNotFoundError(f"BET file not found: {bet_path}")

    output_directory = Path("/home/clarence/Geant4_SAF_Calculations/PHITS_Geant4_python_scripts/5_other_input_files/ICRP_107")
    output_directory.mkdir(parents=True, exist_ok=True)

    csv_path = output_directory / "ICRP-107_combined.csv"
    summary_path = output_directory / "ICRP-107_summary.txt"
    graph_directory = output_directory / "graphs"

    print("\nReading RAD file...")
    rad_records, rad_summary = parse_rad(rad_path)
    print(
        f"  Extracted {len(rad_summary['radionuclides']):,} radionuclides "
        f"and {len(rad_records):,} RAD records."
    )

    print("\nReading BET file...")
    bet_records, bet_summary = parse_bet(bet_path)
    print(
        f"  Extracted {len(bet_summary['radionuclides']):,} radionuclides "
        f"and {len(bet_records):,} BET records."
    )

    # Keep RAD first, then BET.
    combined_records = rad_records + bet_records

    print("\nWriting combined CSV...")
    write_csv(combined_records, csv_path)

    summary_text = build_summary(
        rad_summary,
        bet_summary,
        combined_records,
    )
    summary_path.write_text(summary_text + "\n", encoding="utf-8")

    print_summary(summary_text)

    print("OUTPUT FILES")
    print("-" * 78)
    print(f"Combined CSV: {csv_path.resolve()}")
    print(f"Summary TXT:  {summary_path.resolve()}")
    print(f"Graph folder: {graph_directory.resolve()}")

    interactive_viewer(
        combined_records,
        graph_directory,
    )

    print("\nDone.")


if __name__ == "__main__":
    main()