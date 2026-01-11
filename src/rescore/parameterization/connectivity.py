"""
Ligand extraction with CONECT record preservation for docked complexes.

This module handles text-based extraction of ligands from PDB files while
preserving bond connectivity information encoded in CONECT records.
"""

from pathlib import Path
from typing import Dict, List, Set, Tuple

from loguru import logger

from rescore.parameterization.errors import LigandExtractionError


def extract_ligand_simple(
    pdb_path: Path,
    ligand_resname: str,
    output_pdb: Path,
) -> int:
    """
    Extract ligand from PDB file using simple text parsing (no CONECT records).

    This is a fallback extraction method that doesn't require or preserve CONECT
    records. Useful when docking software doesn't generate connectivity information.

    WARNING: Without CONECT records, antechamber must infer bond orders from
    geometry alone, which often fails for complex ligands, non-standard bonding,
    or aromatic systems.

    Args:
        pdb_path: Path to complex PDB file
        ligand_resname: Three-letter residue name of ligand to extract
        output_pdb: Path where extracted ligand PDB will be written

    Returns:
        Number of atoms extracted

    Raises:
        LigandExtractionError: If no ligand atoms found or file write fails
    """
    logger.info(f"Extracting ligand '{ligand_resname}' using simple text parsing")

    # Read PDB and extract ligand ATOM/HETATM lines
    try:
        with open(pdb_path, 'r') as f:
            pdb_lines = f.readlines()
    except Exception as e:
        raise LigandExtractionError(f"Failed to read PDB file {pdb_path}: {e}")

    # Extract ligand atoms
    ligand_atoms: List[str] = []
    original_serials: List[int] = []

    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        # PDB format: columns 18-20 for residue name (1-indexed, 0-indexed is 17-20)
        resname = line[17:20].strip()
        if resname == ligand_resname:
            ligand_atoms.append(line)
            # Extract serial number (columns 7-11)
            serial = int(line[6:11].strip())
            original_serials.append(serial)

    if not ligand_atoms:
        raise LigandExtractionError(
            f"No atoms found for ligand '{ligand_resname}' in {pdb_path}"
        )

    logger.debug(f"Extracted {len(ligand_atoms)} ligand atoms")

    # Renumber atoms starting from 1
    renumbered_lines: List[str] = []
    for new_serial, line in enumerate(ligand_atoms, start=1):
        # Replace serial number (columns 7-11, right-justified)
        new_line = f"{line[:6]}{new_serial:5d}{line[11:]}"
        renumbered_lines.append(new_line)

    # Write output PDB
    try:
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        with open(output_pdb, 'w') as f:
            for line in renumbered_lines:
                f.write(line)
            f.write(f"TER   {len(renumbered_lines) + 1:5d}      {ligand_resname}   1\n")
            f.write("END   \n")
    except Exception as e:
        raise LigandExtractionError(f"Failed to write ligand PDB to {output_pdb}: {e}")

    logger.success(f"Ligand extracted: {output_pdb} ({len(ligand_atoms)} atoms, no CONECT)")
    return len(ligand_atoms)


def extract_ligand_with_connectivity(
    pdb_path: Path,
    ligand_resname: str,
    output_pdb: Path,
) -> None:
    """
    Extract ligand from PDB file preserving CONECT records with renumbering.

    Performs text-based extraction without BioPython to preserve exact formatting
    and connectivity information. Extracts ATOM/HETATM records for the ligand,
    extracts corresponding CONECT records, renumbers atoms starting from 1,
    and updates all CONECT references.

    STRICT REQUIREMENT: Input PDB must contain CONECT records for the ligand.
    If no CONECT records are found, this function raises an exception.

    Args:
        pdb_path: Path to complex PDB file with CONECT records
        ligand_resname: Three-letter residue name of ligand to extract
        output_pdb: Path where extracted ligand PDB will be written

    Raises:
        LigandExtractionError: If no ligand atoms found, no CONECT records found,
                               or CONECT references invalid atoms
    """
    logger.info(f"Extracting ligand '{ligand_resname}' with connectivity from {pdb_path}")

    # Read entire PDB file
    try:
        with open(pdb_path, 'r') as f:
            pdb_lines = f.readlines()
    except Exception as e:
        raise LigandExtractionError(f"Failed to read PDB file {pdb_path}: {e}")

    # Step 1: Extract ligand ATOM/HETATM lines and build serial number mapping
    ligand_atom_lines = []
    original_serials = []  # List of original serial numbers in order
    old_to_new_serial: Dict[int, int] = {}  # Mapping from old serial to new serial

    for line in pdb_lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        # Extract residue name (columns 18-20, right-justified with spaces)
        resname = line[17:20].strip()
        
        if resname != ligand_resname:
            continue

        # Extract original serial number (columns 7-11)
        try:
            original_serial = int(line[6:11].strip())
        except ValueError:
            logger.warning(f"Invalid atom serial number in line: {line.strip()}")
            continue

        original_serials.append(original_serial)
        new_serial = len(original_serials)  # 1-indexed
        old_to_new_serial[original_serial] = new_serial
        ligand_atom_lines.append(line)

    if not ligand_atom_lines:
        raise LigandExtractionError(
            f"No atoms found for ligand '{ligand_resname}' in {pdb_path}"
        )

    logger.debug(f"Extracted {len(ligand_atom_lines)} ligand atoms")
    logger.debug(f"Original serial numbers: {min(original_serials)} to {max(original_serials)}")

    # Step 2: Extract CONECT records that involve ligand atoms
    ligand_serials_set = set(original_serials)
    conect_records = []

    for line in pdb_lines:
        if not line.startswith("CONECT"):
            continue

        # Parse CONECT record
        # Format: CONECT atom1 atom2 atom3 atom4 ...
        # Columns: 7-11 (atom1), 12-16 (atom2), 17-21 (atom3), 22-26 (atom4), ...
        try:
            parts = line[6:].split()
            serials = [int(s) for s in parts]
        except ValueError:
            logger.warning(f"Invalid CONECT record: {line.strip()}")
            continue

        if not serials:
            continue

        # First serial is the "from" atom
        from_serial = serials[0]
        
        # Only keep CONECT if the from_serial is a ligand atom
        if from_serial not in ligand_serials_set:
            continue

        # Filter bonded atoms to only include ligand atoms
        bonded_to_ligand = [s for s in serials[1:] if s in ligand_serials_set]
        
        if bonded_to_ligand:
            conect_records.append((from_serial, bonded_to_ligand))

    if not conect_records:
        raise LigandExtractionError(
            f"No CONECT records found for ligand '{ligand_resname}' in {pdb_path}\n\n"
            "CONECT records are required for proper bond order assignment in antechamber.\n\n"
            "Docked PDB files must contain connectivity information.\n"
            "If your docking software did not generate CONECT records:\n"
            "  1. Re-dock with connectivity output enabled\n"
            "  2. Use a molecular editor to add bonds and save with CONECT\n"
            "  3. Convert from MOL2/SDF format (which has explicit bonds)\n\n"
            "Without CONECT records, antechamber will fail or produce incorrect topology."
        )

    logger.debug(f"Extracted {len(conect_records)} CONECT records")

    # Step 3: Renumber atoms in ATOM/HETATM lines
    renumbered_atom_lines = []
    
    for idx, line in enumerate(ligand_atom_lines, start=1):
        # Renumber atom serial (columns 7-11, right-justified)
        new_line = line[:6] + f"{idx:5d}" + line[11:]
        renumbered_atom_lines.append(new_line)

    # Step 4: Renumber CONECT records
    renumbered_conect_lines = []
    
    for from_serial, bonded_serials in conect_records:
        # Map old serials to new serials
        new_from = old_to_new_serial[from_serial]
        new_bonded = [old_to_new_serial[s] for s in bonded_serials]
        
        # Build CONECT line
        # Format: "CONECT" + serials in 5-char fields, right-justified
        conect_line = "CONECT"
        conect_line += f"{new_from:5d}"
        for bonded_serial in new_bonded:
            conect_line += f"{bonded_serial:5d}"
        conect_line += "\n"
        
        renumbered_conect_lines.append(conect_line)

    # Step 5: Write output PDB
    try:
        with open(output_pdb, 'w') as f:
            # Write atom records
            f.writelines(renumbered_atom_lines)
            # Write CONECT records
            f.writelines(renumbered_conect_lines)
            # Write END record
            f.write("END\n")
        
        logger.success(
            f"Ligand extracted to {output_pdb}: "
            f"{len(renumbered_atom_lines)} atoms, {len(renumbered_conect_lines)} CONECT records"
        )
        
    except Exception as e:
        raise LigandExtractionError(
            f"Failed to write ligand PDB to {output_pdb}: {e}"
        )


def validate_conect_records(pdb_path: Path) -> Tuple[int, int]:
    """
    Validate CONECT records in a PDB file.

    Counts total CONECT records and checks for dangling references
    (CONECT lines referencing non-existent atom serial numbers).

    Args:
        pdb_path: Path to PDB file to validate

    Returns:
        Tuple of (num_conect_records, num_atoms_referenced)

    Raises:
        LigandExtractionError: If CONECT references invalid atoms
    """
    logger.debug(f"Validating CONECT records in {pdb_path}")

    # Read file and collect atom serials
    atom_serials = set()
    conect_lines = []

    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        serial = int(line[6:11].strip())
                        atom_serials.add(serial)
                    except ValueError:
                        continue
                elif line.startswith("CONECT"):
                    conect_lines.append(line)
    except Exception as e:
        raise LigandExtractionError(f"Failed to read PDB file {pdb_path}: {e}")

    # Validate CONECT references
    referenced_serials = set()
    invalid_references = []

    for line in conect_lines:
        try:
            parts = line[6:].split()
            serials = [int(s) for s in parts]
            
            for serial in serials:
                referenced_serials.add(serial)
                if serial not in atom_serials:
                    invalid_references.append(serial)
        except ValueError:
            continue

    if invalid_references:
        raise LigandExtractionError(
            f"CONECT records reference non-existent atoms: {invalid_references}\n"
            f"Available atom serials: {min(atom_serials)} to {max(atom_serials)}\n"
            f"File: {pdb_path}"
        )

    logger.debug(
        f"Validation passed: {len(conect_lines)} CONECT records, "
        f"{len(referenced_serials)} unique atoms referenced"
    )

    return len(conect_lines), len(referenced_serials)
