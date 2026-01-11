"""
PDB structure validation for MM/GBSA calculations.

This module uses BioPython to parse PDB files and check for:
- Missing backbone atoms (N, CA, C, O) in standard residues
- Presence of metal atoms
- Automatic ligand detection (exactly one required)

LIMITATION: Side-chain heavy atoms are NOT validated. This is a deliberate
design choice to keep validation fast and avoid false positives from
crystallographic artifacts.

All checks are explicit and fail fast with detailed error messages.
"""

from pathlib import Path
from typing import List, Tuple

from Bio.PDB import PDBParser, Structure
from Bio.PDB.Residue import Residue
from loguru import logger

from rescore.config import METAL_ELEMENTS, STANDARD_AMINO_ACIDS, SOLVENT_RESIDUES
from rescore.validation.errors import (
    MissingBackboneAtomsError,
    NoLigandDetectedError,
    MultipleLigandsError,
    MetalDetectedError,
    AltlocOccupancyError,
)


def validate_pdb_complex(pdb_path: Path) -> str:
    """
    Validate a protein-ligand complex PDB file for MM/GBSA compatibility.

    This function automatically detects ligands and runs all validation checks
    in sequence. Raises the first error encountered.

    Ligand detection logic:
    - Standard amino acids (ALA, GLY, etc.) → protein
    - Water and common ions (HOH, NA, CL, etc.) → solvent
    - Everything else → ligand

    Args:
        pdb_path: Path to PDB file containing protein-ligand complex

    Returns:
        Three-letter residue name of the detected ligand

    Raises:
        NoLigandDetectedError: If no ligand molecules are detected
        MultipleLigandsError: If multiple ligand molecules are detected
        MissingBackboneAtomsError: If backbone atoms are missing in standard residues
        MetalDetectedError: If metal atoms are present
        AltlocOccupancyError: If alternate locations or partial occupancy detected
        FileNotFoundError: If PDB file does not exist
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    logger.info(f"Validating PDB complex: {pdb_path}")

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))

    # Run validation checks (fail fast)
    check_missing_backbone_atoms(structure)
    check_for_metals(structure)
    check_altloc_and_occupancy(structure)
    ligand_resname = detect_ligands(structure)

    logger.success(f"PDB validation passed: {pdb_path}")
    logger.info(f"Detected ligand: {ligand_resname}")
    
    return ligand_resname


def check_missing_backbone_atoms(structure: Structure) -> None:
    """
    Check for missing backbone atoms (N, CA, C, O) in standard amino acid residues.

    TERMINUS HANDLING: C-terminal residues are allowed to have OT1/OT2 instead of O.
    This is standard PDB nomenclature for terminal oxygens.

    Validation rules:
    - Internal residues: Must have N, CA, C, O
    - C-terminal residues: Must have N, CA, C, and either O OR (OT1 AND OT2)

    LIMITATION: This function only validates backbone atoms. Side-chain heavy atoms
    are NOT checked in this version. This is a deliberate limitation to keep
    validation fast and deterministic.

    Args:
        structure: BioPython Structure object

    Raises:
        MissingBackboneAtomsError: If backbone atoms are missing in standard residues
    """
    missing_residues: List[Tuple[str, int, str, str]] = []  # chain, resnum, resname, reason

    for model in structure:
        for chain in model:
            # Get all standard amino acid residues in this chain
            chain_residues = [
                res for res in chain 
                if res.get_resname() in STANDARD_AMINO_ACIDS
            ]
            
            if not chain_residues:
                continue
            
            for idx, residue in enumerate(chain_residues):
                resname = residue.get_resname()
                is_c_terminal = (idx == len(chain_residues) - 1)

                # Check if residue is marked as disordered
                if residue.is_disordered():
                    missing_residues.append(
                        (chain.id, residue.id[1], resname, "disordered residue")
                    )
                    continue

                # Get present atoms
                present_atoms = {atom.get_name() for atom in residue}

                # Common backbone atoms for all residues
                required_core = {"N", "CA", "C"}
                
                if is_c_terminal:
                    # C-terminal: require N, CA, C and either O or (OT1 and OT2)
                    missing_core = required_core - present_atoms
                    
                    if missing_core:
                        missing_residues.append(
                            (chain.id, residue.id[1], resname, 
                             f"C-terminal missing {missing_core}")
                        )
                        logger.warning(
                            f"Missing core backbone atoms {missing_core} in "
                            f"C-terminal {resname} {chain.id}:{residue.id[1]}"
                        )
                    elif "O" not in present_atoms:
                        # O is missing, check for terminal oxygens
                        terminal_oxygens = {"OT1", "OT2"}
                        if not terminal_oxygens.issubset(present_atoms):
                            missing = terminal_oxygens - present_atoms
                            missing_residues.append(
                                (chain.id, residue.id[1], resname,
                                 f"C-terminal missing O or {missing}")
                            )
                            logger.warning(
                                f"C-terminal {resname} {chain.id}:{residue.id[1]} "
                                f"missing O and incomplete OT1/OT2 (missing {missing})"
                            )
                else:
                    # Internal residue: require N, CA, C, O
                    required_all = required_core | {"O"}
                    missing = required_all - present_atoms
                    
                    if missing:
                        missing_residues.append(
                            (chain.id, residue.id[1], resname,
                             f"internal missing {missing}")
                        )
                        logger.warning(
                            f"Missing backbone atoms {missing} in "
                            f"internal {resname} {chain.id}:{residue.id[1]}"
                        )

    if missing_residues:
        error_msg = (
            f"Found {len(missing_residues)} residue(s) with missing backbone atoms:\n"
        )
        for chain_id, res_num, resname, reason in missing_residues[:5]:  # Show first 5
            error_msg += f"  - {resname} {chain_id}:{res_num} ({reason})\n"
        if len(missing_residues) > 5:
            error_msg += f"  ... and {len(missing_residues) - 5} more\n"
        error_msg += (
            "\nMM/GBSA requires complete backbone atoms.\n"
            "Internal residues: N, CA, C, O\n"
            "C-terminal residues: N, CA, C, and either O or (OT1 and OT2)\n\n"
            "Note: Side-chain atoms are not validated in this check.\n"
            "Use modeling tools to reconstruct missing atoms before proceeding."
        )
        raise MissingBackboneAtomsError(error_msg)


def check_for_metals(structure: Structure) -> None:
    """
    Check for presence of metal atoms in the structure.

    Metal coordination is not supported by GAFF2, so we reject any
    structure containing metal atoms.

    IMPORTANT: This function skips metal detection in standard amino acid residues
    to avoid false positives from mislabeled element columns (e.g., "CA" element
    for alpha carbon atoms). Only HETATM records and non-standard residues are checked.

    Args:
        structure: BioPython Structure object

    Raises:
        MetalDetectedError: If any metal atoms are found in non-standard residues
    """
    detected_metals: List[Tuple[str, str, int]] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()

                # Skip standard amino acids - they cannot contain actual metals
                # This prevents false positives from mislabeled element columns
                if resname in STANDARD_AMINO_ACIDS:
                    continue

                # Check if residue name matches known metals
                if resname in METAL_ELEMENTS:
                    detected_metals.append(
                        (resname, chain.id, residue.id[1])
                    )
                    continue

                # Also check individual atom elements in non-standard residues
                for atom in residue:
                    element = atom.element.upper()
                    if element in METAL_ELEMENTS:
                        detected_metals.append(
                            (element, chain.id, residue.id[1])
                        )

    if detected_metals:
        error_msg = f"Detected {len(detected_metals)} metal atom(s):\n"
        for element, chain_id, res_num in detected_metals[:10]:  # Show first 10
            error_msg += f"  - {element} in chain {chain_id} residue {res_num}\n"
        if len(detected_metals) > 10:
            error_msg += f"  ... and {len(detected_metals) - 10} more\n"
        error_msg += (
            "\nMetal coordination is not supported by GAFF2. "
            "Remove metal atoms or use a different parameterization tool."
        )
        raise MetalDetectedError(error_msg)


def check_altloc_and_occupancy(structure: Structure) -> None:
    """
    Check for alternate locations (ALTLOC) and partial occupancy.

    AmberTools requires single-conformation structures with full occupancy.
    Alternate locations and partial occupancies break tleap, antechamber,
    and MM/GBSA calculations.

    NOTE: Occupancy values of 0.00 are ignored, as these typically indicate
    missing/unspecified occupancy in improperly formatted PDB files and are
    not relevant for AmberTools parameterization.

    Args:
        structure: BioPython Structure object

    Raises:
        AltlocOccupancyError: If any alternate locations or partial occupancy detected
    """
    altloc_atoms: List[Tuple[str, str, int, str, str]] = []  # (chain, resname, resnum, atom, altloc)
    partial_occupancy_atoms: List[Tuple[str, str, int, str, float]] = []  # (chain, resname, resnum, atom, occ)

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                res_id = residue.id[1]

                for atom in residue:
                    atom_name = atom.get_name()
                    altloc = atom.get_altloc()
                    occupancy = atom.get_occupancy()

                    # Check for alternate location
                    if altloc != ' ' and altloc != '':
                        altloc_atoms.append(
                            (chain.id, resname, res_id, atom_name, altloc)
                        )

                    # Check for partial occupancy
                    # NOTE: 0.0 occupancy is treated as "unspecified" and ignored
                    # Only flag occupancies between 0.0 (exclusive) and 1.0 (exclusive)
                    if 0.0 < occupancy < 1.0:
                        partial_occupancy_atoms.append(
                            (chain.id, resname, res_id, atom_name, occupancy)
                        )

    # Build error message if issues detected
    if altloc_atoms or partial_occupancy_atoms:
        error_msg = ""

        if altloc_atoms:
            error_msg += f"Detected {len(altloc_atoms)} atom(s) with alternate locations:\n"
            for chain_id, resname, res_num, atom_name, altloc in altloc_atoms[:10]:
                error_msg += f"  - {resname} {chain_id}:{res_num} atom {atom_name} (altloc '{altloc}')\n"
            if len(altloc_atoms) > 10:
                error_msg += f"  ... and {len(altloc_atoms) - 10} more\n"
            error_msg += "\n"

        if partial_occupancy_atoms:
            error_msg += f"Detected {len(partial_occupancy_atoms)} atom(s) with partial occupancy:\n"
            for chain_id, resname, res_num, atom_name, occ in partial_occupancy_atoms[:10]:
                error_msg += f"  - {resname} {chain_id}:{res_num} atom {atom_name} (occupancy {occ:.2f})\n"
            if len(partial_occupancy_atoms) > 10:
                error_msg += f"  ... and {len(partial_occupancy_atoms) - 10} more\n"
            error_msg += "\n"

        error_msg += (
            "AmberTools requires single-conformation structures with full occupancy.\n"
            "Alternate locations and partial occupancies break parameterization and MM/GBSA.\n\n"
            "Resolution options:\n"
            "  1. Use a crystallographic refinement tool to select one conformer\n"
            "  2. Remove alternate location indicators and set occupancy to 1.0\n"
            "  3. Manually edit the PDB file to resolve conformational ambiguity\n\n"
            "Tools: phenix.pdbtools, PyMOL (remove alt, alter occupancy), or pdb-tools."
        )

        raise AltlocOccupancyError(error_msg)


def detect_ligands(structure: Structure) -> str:
    """
    Automatically detect ligand molecules in a protein-ligand complex.

    Detection logic:
    - Standard amino acids (from STANDARD_AMINO_ACIDS) → classified as protein
    - Water and common solvents/ions (from SOLVENT_RESIDUES) → classified as solvent
    - Everything else → classified as ligand

    Exactly one ligand molecule is required. Multiple copies or zero ligands
    are rejected explicitly.

    Args:
        structure: BioPython Structure object

    Returns:
        Three-letter residue name of the single detected ligand

    Raises:
        NoLigandDetectedError: If no ligand molecules are detected
        MultipleLigandsError: If multiple ligand molecules are detected
    """
    # Collect all ligand residues (resname, chain, resnum)
    ligands: List[Tuple[str, str, int]] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()

                # Skip standard protein residues
                if resname in STANDARD_AMINO_ACIDS:
                    continue

                # Skip water and common solvents/ions
                if resname in SOLVENT_RESIDUES:
                    continue

                # Everything else is a potential ligand
                ligands.append((resname, chain.id, residue.id[1]))

    # Enforce exactly one ligand
    if len(ligands) == 0:
        error_msg = (
            "No ligand molecules detected in structure.\n\n"
            "Ligand detection assumes:\n"
            "  - Standard amino acids → protein\n"
            "  - Water/ions (HOH, NA, CL, etc.) → solvent\n"
            "  - Everything else → ligand\n\n"
            "If your structure contains a ligand, ensure it has a residue name\n"
            "distinct from standard amino acids and common solvents."
        )
        raise NoLigandDetectedError(error_msg)

    if len(ligands) > 1:
        # Group by residue name to provide detailed error message
        resname_groups: dict[str, List[Tuple[str, int]]] = {}
        for resname, chain_id, res_num in ligands:
            if resname not in resname_groups:
                resname_groups[resname] = []
            resname_groups[resname].append((chain_id, res_num))

        error_msg = f"Detected {len(ligands)} ligand molecule(s):\n\n"
        for resname, locations in resname_groups.items():
            error_msg += f"  {resname} ({len(locations)} copy/copies):\n"
            for chain_id, res_num in locations[:5]:  # Show first 5 per type
                error_msg += f"    - Chain {chain_id} residue {res_num}\n"
            if len(locations) > 5:
                error_msg += f"    ... and {len(locations) - 5} more\n"

        error_msg += (
            "\nOnly single-ligand systems are supported.\n"
            "Remove duplicate ligands or process them separately.\n"
            "If multiple residues are part of the same ligand, they must share\n"
            "the same three-letter residue name."
        )
        raise MultipleLigandsError(error_msg)

    # Exactly one ligand found
    ligand_resname, chain_id, res_num = ligands[0]
    logger.info(
        f"Detected ligand: {ligand_resname} (chain {chain_id}, residue {res_num})"
    )
    return ligand_resname
