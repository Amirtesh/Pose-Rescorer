"""
Coordinate transfer from docked PDB poses to MOL2/SDF chemistry files.

This module handles transferring docked coordinates from PDB files to ligand
chemistry files (MOL2/SDF) that contain complete bond information. This separation
is critical because:

1. PDB files from docking contain optimized coordinates but no bond orders
2. Antechamber requires bond information for GAFF2 parameterization
3. Inferring bonds from PDB geometry is unreliable and scientifically incorrect
4. MOL2/SDF files contain explicit bond orders from chemical preparation

The workflow matches heavy atoms by element and order, transferring only coordinates
while preserving all bond topology from the chemistry file.
"""

from pathlib import Path
from typing import Dict, List, Tuple

from loguru import logger

from mmgbsa.parameterization.errors import LigandExtractionError


def parse_mol2_atoms(mol2_path: Path) -> List[Dict]:
    """
    Parse atom records from MOL2 file.

    Args:
        mol2_path: Path to MOL2 file

    Returns:
        List of atom dictionaries with keys: serial, name, x, y, z, type, element

    Raises:
        LigandExtractionError: If MOL2 parsing fails
    """
    try:
        with open(mol2_path, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        raise LigandExtractionError(f"Failed to read MOL2 file {mol2_path}: {e}")

    atoms = []
    in_atom_section = False

    for line in lines:
        line = line.strip()
        
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue
        elif line.startswith("@<TRIPOS>"):
            in_atom_section = False
            continue
        
        if in_atom_section and line:
            parts = line.split()
            if len(parts) >= 6:
                serial = int(parts[0])
                name = parts[1]
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                atom_type = parts[5]
                
                # Extract element from atom type (e.g., "C.3" -> "C", "O.3" -> "O")
                element = atom_type.split('.')[0]
                
                atoms.append({
                    'serial': serial,
                    'name': name,
                    'x': x,
                    'y': y,
                    'z': z,
                    'type': atom_type,
                    'element': element,
                })

    if not atoms:
        raise LigandExtractionError(f"No atoms found in MOL2 file: {mol2_path}")

    logger.debug(f"Parsed {len(atoms)} atoms from MOL2")
    return atoms


def parse_pdb_ligand_atoms(pdb_path: Path, ligand_resname: str) -> List[Dict]:
    """
    Parse ligand atom coordinates from PDB file.

    Args:
        pdb_path: Path to docked PDB file
        ligand_resname: Three-letter residue name of ligand

    Returns:
        List of atom dictionaries with keys: serial, name, x, y, z, element

    Raises:
        LigandExtractionError: If no ligand atoms found
    """
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        raise LigandExtractionError(f"Failed to read PDB file {pdb_path}: {e}")

    atoms = []
    
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        
        # PDB format parsing
        resname = line[17:20].strip()
        if resname != ligand_resname:
            continue
        
        serial = int(line[6:11].strip())
        name = line[12:16].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        element = line[76:78].strip() if len(line) > 76 else ""
        
        # If element not in columns 76-78, try to infer from atom name
        if not element:
            # Remove digits and take first 1-2 characters
            clean_name = ''.join(c for c in name if not c.isdigit())
            element = clean_name[0:2] if len(clean_name) > 1 and clean_name[1].islower() else clean_name[0]
        
        atoms.append({
            'serial': serial,
            'name': name,
            'x': x,
            'y': y,
            'z': z,
            'element': element.upper(),
        })

    if not atoms:
        raise LigandExtractionError(
            f"No atoms found for ligand '{ligand_resname}' in PDB: {pdb_path}"
        )

    logger.debug(f"Parsed {len(atoms)} ligand atoms from PDB")
    return atoms


def transfer_coordinates(
    pdb_path: Path,
    ligand_resname: str,
    mol2_path: Path,
    output_mol2: Path,
) -> None:
    """
    Transfer docked coordinates from PDB to MOL2 file.

    Matches heavy atoms (non-hydrogen) between PDB and MOL2 by element and order,
    then transfers the docked coordinates from PDB to MOL2 while preserving all
    bond information and hydrogen atoms from MOL2.

    Strategy:
    1. Extract heavy atoms from both PDB and MOL2 in file order
    2. Validate element sequences match exactly
    3. Transfer coordinates from PDB heavy atoms to MOL2 heavy atoms
    4. Keep all MOL2 hydrogens with original coordinates
    5. Preserve all bond records unchanged

    Args:
        pdb_path: Path to docked complex PDB file
        ligand_resname: Three-letter residue name of ligand
        mol2_path: Path to ligand MOL2 file with complete chemistry
        output_mol2: Path for output MOL2 with docked coordinates

    Raises:
        LigandExtractionError: If atom counts mismatch, elements don't match,
                               or coordinate transfer fails
    """
    logger.info(f"Transferring coordinates from PDB to MOL2")
    logger.debug(f"  PDB source: {pdb_path}")
    logger.debug(f"  MOL2 chemistry: {mol2_path}")
    
    # Parse both files
    pdb_atoms = parse_pdb_ligand_atoms(pdb_path, ligand_resname)
    mol2_atoms = parse_mol2_atoms(mol2_path)
    
    # Extract heavy atoms (non-H) from both
    pdb_heavy = [a for a in pdb_atoms if a['element'] != 'H']
    mol2_heavy = [a for a in mol2_atoms if a['element'] != 'H']
    
    logger.info(f"PDB: {len(pdb_atoms)} atoms ({len(pdb_heavy)} heavy)")
    logger.info(f"MOL2: {len(mol2_atoms)} atoms ({len(mol2_heavy)} heavy)")
    
    # Validate heavy atom counts match
    if len(pdb_heavy) != len(mol2_heavy):
        raise LigandExtractionError(
            f"Heavy atom count mismatch:\n"
            f"  PDB has {len(pdb_heavy)} heavy atoms\n"
            f"  MOL2 has {len(mol2_heavy)} heavy atoms\n"
            f"Coordinate transfer requires identical heavy atom count.\n"
            f"Ensure the docked PDB and MOL2 chemistry file represent the same molecule."
        )
    
    # Validate element sequences match
    pdb_elements = [a['element'] for a in pdb_heavy]
    mol2_elements = [a['element'] for a in mol2_heavy]
    
    if pdb_elements != mol2_elements:
        # Show first few mismatches
        mismatches = []
        for i, (p, m) in enumerate(zip(pdb_elements, mol2_elements), 1):
            if p != m:
                mismatches.append(f"  Position {i}: PDB={p}, MOL2={m}")
                if len(mismatches) >= 5:
                    mismatches.append("  ...")
                    break
        
        raise LigandExtractionError(
            f"Heavy atom element sequence mismatch:\n"
            f"  PDB order: {' '.join(pdb_elements[:10])}{'...' if len(pdb_elements) > 10 else ''}\n"
            f"  MOL2 order: {' '.join(mol2_elements[:10])}{'...' if len(mol2_elements) > 10 else ''}\n"
            f"\nMismatches:\n" + "\n".join(mismatches) + "\n\n"
            f"This indicates the PDB and MOL2 represent different molecules or atom orderings.\n"
            f"Ensure both files are for the same ligand with consistent atom ordering."
        )
    
    logger.success("Heavy atom validation passed - elements match")
    
    # Build coordinate mapping: mol2_serial -> (x, y, z)
    coord_map = {}
    for pdb_atom, mol2_atom in zip(pdb_heavy, mol2_heavy):
        coord_map[mol2_atom['serial']] = (pdb_atom['x'], pdb_atom['y'], pdb_atom['z'])
        logger.debug(
            f"  Map MOL2 atom {mol2_atom['serial']} ({mol2_atom['element']}) -> "
            f"PDB coords ({pdb_atom['x']:.3f}, {pdb_atom['y']:.3f}, {pdb_atom['z']:.3f})"
        )
    
    logger.info(f"Mapped {len(coord_map)} heavy atom coordinates")
    
    # Read original MOL2 file
    with open(mol2_path, 'r') as f:
        mol2_lines = f.readlines()
    
    # Write new MOL2 with updated coordinates
    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_mol2, 'w') as f:
        in_atom_section = False
        
        for line in mol2_lines:
            stripped = line.strip()
            
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                f.write(line)
                continue
            elif stripped.startswith("@<TRIPOS>"):
                in_atom_section = False
                f.write(line)
                continue
            
            if in_atom_section and stripped:
                parts = line.split()
                if len(parts) >= 6:
                    serial = int(parts[0])
                    
                    # If this atom has mapped coordinates, use them
                    if serial in coord_map:
                        x, y, z = coord_map[serial]
                        # Preserve MOL2 formatting: reconstruct line with new coords
                        new_line = f"{parts[0]:>7} {parts[1]:>4} {x:>9.4f} {y:>9.4f} {z:>9.4f}"
                        if len(parts) > 5:
                            new_line += " " + " ".join(parts[5:])
                        f.write(new_line + "\n")
                    else:
                        # Hydrogen - keep original coordinates
                        f.write(line)
                else:
                    f.write(line)
            else:
                f.write(line)
    
    logger.success(f"Coordinates transferred to: {output_mol2}")
    logger.info(f"  Heavy atoms updated: {len(coord_map)}")
    logger.info(f"  Hydrogens preserved: {len(mol2_atoms) - len(mol2_heavy)}")
