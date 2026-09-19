# Add phantoms in this registry to make them available for use in the simulation framework.

from dataclasses import dataclass
from pathlib import Path

import b_config.a_config as config

def get_phantom(code: str) -> PhantomSpec:
    try:
        return PHANTOMS[code]
    except KeyError:
        raise ValueError(
            f"Unknown phantom: {code}"
        )

def get_phantom_group(code: str) -> tuple[str, ...]:

    if code in PHANTOM_GROUPS:
        return PHANTOM_GROUPS[code].phantoms

    if code in PHANTOMS:
        return (code,)

    raise ValueError(
        f"Unknown phantom or phantom group: {code}"
    )

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

@dataclass(frozen=True)
class PhantomGroupSpec:
    code: str
    display_name: str
    phantoms: tuple[str, ...]

PHANTOMS = {
    "MRCP_AF": PhantomSpec(
        code="MRCP_AF",
        display_name="Adult ICRP 145 Female",
        family="MRCP",
        sex="AF",
        cell_file=config.INCLUDE_FILES_DIR / "MRCP-AF.cell",
        material_file=config.INCLUDE_FILES_DIR / "MRCP-AF.material",
        node_file=config.INCLUDE_FILES_DIR / "MRCP-AF.node",
        element_file=config.INCLUDE_FILES_DIR / "MRCP-AF.ele",
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
        cell_file=config.INCLUDE_FILES_DIR / "MRCP-AM.cell",
        material_file=config.INCLUDE_FILES_DIR / "MRCP-AM.material",
        node_file=config.INCLUDE_FILES_DIR / "MRCP-AM.node",
        element_file=config.INCLUDE_FILES_DIR / "MRCP-AM.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AF_H150W55": PhantomSpec(
        code="F_H150W55",
        display_name="Adult ICRP 145 Female Resized to Filipino",
        family="MRCP_Resized",
        sex="AF",
        cell_file=config.INCLUDE_FILES_DIR / "F_H150W55.cell",
        material_file=config.INCLUDE_FILES_DIR / "F_H150W55.material",
        node_file=config.INCLUDE_FILES_DIR / "F_H150W55.node",
        element_file=config.INCLUDE_FILES_DIR / "F_H150W55.ele",
        skeletal_mass_columns=(
            "Ref_AF_Marrow_Mass(g)",
            "Ref_AF_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AM_H165W65": PhantomSpec(
        code="M_H165W65",
        display_name="Adult ICRP 145 Male Resized to Filipino",
        family="MRCP_Resized",
        sex="AM",
        cell_file=config.INCLUDE_FILES_DIR / "M_H165W65.cell",
        material_file=config.INCLUDE_FILES_DIR / "M_H165W65.material",
        node_file=config.INCLUDE_FILES_DIR / "M_H165W65.node",
        element_file=config.INCLUDE_FILES_DIR / "M_H165W65.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AF_Filipino": PhantomSpec(
        code="MRCP_AF_Filipino",
        display_name="Adult Filipino Female",
        family="MRCP_Filipino",
        sex="AF",
        cell_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AF.cell",
        material_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AF.material",
        node_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AF.node",
        element_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AF.ele",
        skeletal_mass_columns=(
            "Ref_AF_Marrow_Mass(g)",
            "Ref_AF_Endosteum_Mass(g)",
        ),
    ),

    "MRCP_AM_Filipino": PhantomSpec(
        code="MRCP_AM_Filipino",
        display_name="Adult Filipino Male",
        family="MRCP_Filipino",
        sex="AM",
        cell_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AM.cell",
        material_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AM.material",
        node_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AM.node",
        element_file=config.INCLUDE_FILES_DIR / "Filipino-MRCP-AM.ele",
        skeletal_mass_columns=(
            "Ref_AM_Marrow_Mass(g)",
            "Ref_AM_Endosteum_Mass(g)",
        ),
    ),
}

PHANTOM_GROUPS = {

    "MRCP_AF_AM": PhantomGroupSpec(
        code="MRCP_AF_AM",
        display_name="Both ICRP 145 Phantoms",
        phantoms=(
            "MRCP_AF",
            "MRCP_AM",
        ),
    ),

    "MRCP_AF_AM_Filipino_Resized": PhantomGroupSpec(
        code="MRCP_AF_AM_Filipino_Resized",
        display_name="Both ICRP 145 Phantoms Resized to Filipino",
        phantoms=(
            "MRCP_6AF_H150W55",
            "MRCP_AM_H15W65",
        ),
    ),

    "MRCP_AF_AM_Filipino": PhantomGroupSpec(
        code="MRCP_AF_AM_Filipino",
        display_name="Both Filipino MRCP Phantoms",
        phantoms=(
            "MRCP_AF_Filipino",
            "MRCP_AM_Filipino",
        ),
    ),
}