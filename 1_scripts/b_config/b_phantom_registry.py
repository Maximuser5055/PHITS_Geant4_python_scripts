# Add phantoms in this registry to make them available for use in the simulation framework.

from dataclasses import dataclass
from pathlib import Path

import b_config.a_config as config

@dataclass(frozen=True)
class PhantomSpec:
    code: str
    display_name: str
    family: str
    sex: str
    cell_file: Path
    material_file: Path
    node_file: Path
    element_file: Path
    skeletal_mass_columns: tuple[str, str] | None = None


PHANTOMS = {
    "MRCP_AF": PhantomSpec(
        code="MRCP_AF",
        display_name="Adult ICRP 145 Female",
        family="MRCP",
        sex="AF",
        cell_file=config.ROOT / "2_phits/phantoms/MRCP-AF.cell",
        material_file=config.ROOT / "2_phits/phantoms/MRCP-AF.material",
        node_file=config.ROOT / "2_phits/phantoms/MRCP-AF.node",
        element_file=config.ROOT / "2_phits/phantoms/MRCP-AF.ele",
        skeletal_mass_columns=(
            "Ref_AF_Marrow_Mass(g)",
            "Ref_AF_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AM": PhantomSpec(
        code="MRCP_AM",
        display_name="Adult ICRP 145 Male",
        family="MRCP",
        sex="AM",
        cell_file=config.ROOT / "2_phits/phantoms/MRCP-AM.cell",
        material_file=config.ROOT / "2_phits/phantoms/MRCP-AM.material",
        node_file=config.ROOT / "2_phits/phantoms/MRCP-AM.node",
        element_file=config.ROOT / "2_phits/phantoms/MRCP-AM.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AF_Resized": PhantomSpec(
        code="F_H150W55",
        display_name="Adult ICRP 145 Female Resized to Filipino",
        family="MRCP_Resized",
        sex="AF",
        cell_file=config.ROOT / "2_phits/phantoms/F_H150W55.cell",
        material_file=config.ROOT / "2_phits/phantoms/F_H150W55.material",
        node_file=config.ROOT / "2_phits/phantoms/F_H150W55.node",
        element_file=config.ROOT / "2_phits/phantoms/F_H150W55.ele",
        skeletal_mass_columns=(
            "Ref_AF_Marrow_Mass(g)",
            "Ref_AF_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AM_Resized": PhantomSpec(
        code="M_H165W65",
        display_name="Adult ICRP 145 Male Resized to Filipino",
        family="MRCP_Resized",
        sex="AM",
        cell_file=config.ROOT / "2_phits/phantoms/M_H165W65.cell",
        material_file=config.ROOT / "2_phits/phantoms/M_H165W65.material",
        node_file=config.ROOT / "2_phits/phantoms/M_H165W65.node",
        element_file=config.ROOT / "2_phits/phantoms/M_H165W65.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),

    "MFCP_AF": PhantomSpec(
        code="MFCP_AF",
        display_name="Adult Filipino Female",
        family="MFCP",
        sex="AF",
        cell_file=config.ROOT / "2_phits/phantoms/MFCP-AF.cell",
        material_file=config.ROOT / "2_phits/phantoms/MFCP-AF.material",
        node_file=config.ROOT / "2_phits/phantoms/MFCP-AF.node",
        element_file=config.ROOT / "2_phits/phantoms/MFCP-AF.ele",
        skeletal_mass_columns=(
            "Ref_AF_Marrow_Mass(g)",
            "Ref_AF_Endosteum_Mass(g)",
        ),
    ),

    "MFCP_AM": PhantomSpec(
        code="MFCP_AM",
        display_name="Adult Filipino Male",
        family="MFCP",
        sex="AM",
        cell_file=config.ROOT / "2_phits/phantoms/MFCP-AM.cell",
        material_file=config.ROOT / "2_phits/phantoms/MFCP-AM.material",
        node_file=config.ROOT / "2_phits/phantoms/MFCP-AM.node",
        element_file=config.ROOT / "2_phits/phantoms/MFCP-AM.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),
}