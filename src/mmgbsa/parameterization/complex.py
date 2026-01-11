"""
Complex assembly module for MM/GBSA preparation.

This module combines protein and ligand topologies into a single-frame complex
suitable for MM/GBSA single-frame rescoring. It does NOT perform:
- Solvation
- Minimization
- MD simulation
- Protonation state inference

The philosophy is strict separation of concerns:
- Protein chemistry: from protein.prmtop
- Ligand chemistry: from ligand.mol2 + ligand.frcmod
- Pose geometry: from docked complex PDB

This ensures reproducible, deterministic complex preparation.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from loguru import logger

from mmgbsa.config import PROTEIN_FF, LIGAND_FF
from mmgbsa.parameterization.errors import ParameterizationError


class ComplexAssemblyError(ParameterizationError):
    """Raised when complex assembly fails."""
    pass


def validate_inputs(
    protein_dir: Path,
    ligand_dir: Path,
    complex_pdb: Path,
) -> Tuple[Path, Path, Path, Path]:
    """
    Validate that all required input files exist.
    
    Args:
        protein_dir: Directory containing protein parameters
        ligand_dir: Directory containing ligand parameters
        complex_pdb: Path to docked complex PDB file
        
    Returns:
        Tuple of (protein_prmtop, protein_inpcrd, ligand_mol2, ligand_frcmod)
        
    Raises:
        ComplexAssemblyError: If any required files are missing
    """
    logger.info("Validating input files...")
    
    # Check protein files
    protein_prmtop = protein_dir / "protein.prmtop"
    protein_inpcrd = protein_dir / "protein.inpcrd"
    
    if not protein_prmtop.exists():
        raise ComplexAssemblyError(
            f"Missing protein topology: {protein_prmtop}\n"
            f"Run: mmgbsa prep-protein <receptor.pdb> -o {protein_dir}"
        )
    
    if not protein_inpcrd.exists():
        raise ComplexAssemblyError(
            f"Missing protein coordinates: {protein_inpcrd}\n"
            f"Run: mmgbsa prep-protein <receptor.pdb> -o {protein_dir}"
        )
    
    # Check ligand files
    ligand_mol2 = ligand_dir / "ligand.mol2"
    ligand_frcmod = ligand_dir / "ligand.frcmod"
    
    if not ligand_mol2.exists():
        raise ComplexAssemblyError(
            f"Missing ligand structure: {ligand_mol2}\n"
            f"Run: mmgbsa parameterize <ligand.mol2> -o {ligand_dir}"
        )
    
    if not ligand_frcmod.exists():
        raise ComplexAssemblyError(
            f"Missing ligand parameters: {ligand_frcmod}\n"
            f"Run: mmgbsa parameterize <ligand.mol2> -o {ligand_dir}"
        )
    
    # Check complex PDB
    if not complex_pdb.exists():
        raise ComplexAssemblyError(
            f"Docked complex PDB not found: {complex_pdb}"
        )
    
    logger.success("All input files validated")
    return protein_prmtop, protein_inpcrd, ligand_mol2, ligand_frcmod


def extract_ligand_coords_from_pdb(
    complex_pdb: Path,
    ligand_resname: str = "LIG",
) -> List[Tuple[str, float, float, float]]:
    """
    Extract ligand atom coordinates from docked complex PDB.
    
    Args:
        complex_pdb: Path to docked complex PDB
        ligand_resname: Residue name of ligand (default: LIG)
        
    Returns:
        List of (atom_name, x, y, z) tuples in PDB order
        
    Raises:
        ComplexAssemblyError: If ligand not found or parsing fails
    """
    logger.info(f"Extracting ligand coordinates from {complex_pdb.name}...")
    logger.info(f"  Looking for residue: {ligand_resname}")
    
    coords = []
    
    with open(complex_pdb, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            
            # Parse PDB ATOM/HETATM record
            resname = line[17:20].strip()
            
            if resname == ligand_resname:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((atom_name, x, y, z))
    
    if not coords:
        raise ComplexAssemblyError(
            f"Ligand residue '{ligand_resname}' not found in {complex_pdb}\n"
            f"Check that the complex PDB contains the ligand"
        )
    
    logger.success(f"Extracted {len(coords)} ligand atoms from {complex_pdb.name}")
    return coords


def update_mol2_coordinates(
    ligand_mol2: Path,
    docked_coords: List[Tuple[str, float, float, float]],
    output_mol2: Path,
) -> None:
    """
    Update MOL2 file with docked coordinates from complex PDB.
    
    This preserves ligand chemistry (atom types, bonds, charges) from the
    parameterized MOL2 while using the docked pose geometry.
    
    Args:
        ligand_mol2: Original parameterized ligand MOL2
        docked_coords: List of (atom_name, x, y, z) from complex PDB
        output_mol2: Path where updated MOL2 will be written
        
    Raises:
        ComplexAssemblyError: If atom count mismatch or parsing fails
    """
    logger.info("Updating ligand MOL2 with docked coordinates...")
    logger.info(f"  Input MOL2:  {ligand_mol2.name}")
    logger.info(f"  Output MOL2: {output_mol2.name}")
    
    # Read original MOL2
    with open(ligand_mol2, "r") as f:
        mol2_lines = f.readlines()
    
    # Build coordinate lookup by atom name
    coord_map = {name: (x, y, z) for name, x, y, z in docked_coords}
    
    # Find @<TRIPOS>ATOM section
    atom_start = None
    atom_end = None
    
    for i, line in enumerate(mol2_lines):
        if "@<TRIPOS>ATOM" in line:
            atom_start = i + 1
        elif atom_start is not None and line.startswith("@<TRIPOS>"):
            atom_end = i
            break
    
    if atom_start is None:
        raise ComplexAssemblyError(
            f"Invalid MOL2 file: no @<TRIPOS>ATOM section in {ligand_mol2}"
        )
    
    if atom_end is None:
        atom_end = len(mol2_lines)
    
    # Count atoms in MOL2
    mol2_atom_count = atom_end - atom_start
    pdb_atom_count = len(docked_coords)
    
    if mol2_atom_count != pdb_atom_count:
        raise ComplexAssemblyError(
            f"Atom count mismatch:\n"
            f"  MOL2 has {mol2_atom_count} atoms\n"
            f"  Complex PDB has {pdb_atom_count} ligand atoms\n"
            f"Ensure ligand parameterization and docking used the same structure"
        )
    
    # Update coordinates in MOL2 atom lines
    updated_lines = mol2_lines[:atom_start]
    
    for line in mol2_lines[atom_start:atom_end]:
        parts = line.split()
        if len(parts) < 9:
            updated_lines.append(line)
            continue
        
        # MOL2 format: atom_id atom_name x y z atom_type ...
        atom_name = parts[1]
        
        if atom_name not in coord_map:
            # Try matching without numbers (e.g., C1 -> C, H2A -> H)
            base_name = re.sub(r'\d+', '', atom_name)
            matching = [name for name in coord_map.keys() if re.sub(r'\d+', '', name) == base_name]
            
            if len(matching) == 1:
                atom_name = matching[0]
            else:
                raise ComplexAssemblyError(
                    f"Cannot match MOL2 atom '{atom_name}' to complex PDB atoms.\n"
                    f"Available PDB atoms: {', '.join(coord_map.keys())}"
                )
        
        x, y, z = coord_map[atom_name]
        
        # Rebuild line with updated coordinates
        # Format: %7d %8s %9.4f %9.4f %9.4f %8s %5s %8s %9.4f
        updated_line = (
            f"{parts[0]:>7s} {parts[1]:>8s} "
            f"{x:9.4f} {y:9.4f} {z:9.4f} "
            f"{' '.join(parts[5:])}\n"
        )
        updated_lines.append(updated_line)
    
    # Append remaining sections
    updated_lines.extend(mol2_lines[atom_end:])
    
    # Write updated MOL2
    with open(output_mol2, "w") as f:
        f.writelines(updated_lines)
    
    logger.success(f"Updated {len(docked_coords)} atom coordinates in {output_mol2.name}")


def assemble_complex(
    protein_prmtop: Path,
    protein_inpcrd: Path,
    ligand_mol2: Path,
    ligand_frcmod: Path,
    output_dir: Path,
) -> None:
    """
    Assemble protein-ligand complex using tleap.
    
    This combines:
    - Protein topology (ff14SB)
    - Ligand topology (GAFF2)
    - Ligand coordinates from MOL2
    
    Note: We use ambpdb to convert protein coordinates back to PDB format
    so tleap can reload and combine with ligand.
    
    No solvation, minimization, or MD is performed.
    
    Args:
        protein_prmtop: Protein topology file
        protein_inpcrd: Protein coordinates
        ligand_mol2: Ligand MOL2 with coordinates
        ligand_frcmod: Ligand force field modifications
        output_dir: Directory for output files
        
    Raises:
        ComplexAssemblyError: If tleap fails
    """
    logger.info("Assembling protein-ligand complex with tleap...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Convert protein inpcrd to PDB using ambpdb
    protein_pdb = output_dir / "protein.pdb"
    logger.info("Converting protein coordinates to PDB format...")
    
    try:
        with open(protein_pdb, "w") as pdb_out:
            result = subprocess.run(
                ["ambpdb", "-p", str(protein_prmtop), "-c", str(protein_inpcrd)],
                stdout=pdb_out,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        
        if result.returncode != 0 or not protein_pdb.exists():
            raise ComplexAssemblyError(
                f"ambpdb failed to convert protein coordinates:\n{result.stderr}"
            )
        
        logger.success(f"Protein PDB generated: {protein_pdb}")
        
    except FileNotFoundError:
        raise ComplexAssemblyError(
            "ambpdb command not found. Ensure AmberTools is properly installed."
        )
    
    # Step 2: Generate tleap input script
    tleap_input = output_dir / "tleap_complex.in"
    complex_prmtop = output_dir / "complex.prmtop"
    complex_inpcrd = output_dir / "complex.inpcrd"
    
    # Use relative paths within output directory
    rel_ligand_frcmod = Path("..") / ligand_frcmod.relative_to(output_dir.parent)
    rel_ligand_mol2 = Path("..") / ligand_mol2.relative_to(output_dir.parent)
    
    tleap_script = f"""# Complex assembly for MM/GBSA
# Generated by mmgbsa.parameterization.complex

# Load force fields
source leaprc.protein.{PROTEIN_FF}
source leaprc.{LIGAND_FF}

# Load ligand parameters
loadamberparams {rel_ligand_frcmod}

# Load protein from PDB (coordinates from inpcrd conversion)
protein = loadpdb protein.pdb

# Load ligand with coordinates from MOL2
ligand = loadmol2 {rel_ligand_mol2}

# Combine protein + ligand
complex = combine {{protein ligand}}

# Check structure
check complex

# Save complex topology
saveamberparm complex complex.prmtop complex.inpcrd

# Exit
quit
"""
    
    with open(tleap_input, "w") as f:
        f.write(tleap_script)
    
    logger.debug(f"Generated tleap input: {tleap_input}")
    
    # Step 3: Run tleap from output directory
    logger.info("Running tleap on complex...")
    
    result = subprocess.run(
        ["tleap", "-f", "tleap_complex.in"],
        cwd=output_dir,
        capture_output=True,
        text=True,
    )
    
    # Save log
    log_file = output_dir / "tleap_complex.log"
    with open(log_file, "w") as f:
        f.write(result.stdout)
        if result.stderr:
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)
    
    logger.debug(f"tleap log saved: {log_file}")
    
    # Check for output files
    if not complex_prmtop.exists() or not complex_inpcrd.exists():
        # Parse errors from log
        fatal_errors = []
        warnings = []
        
        for line in result.stdout.splitlines():
            if "FATAL" in line or "ERROR" in line:
                fatal_errors.append(line.strip())
            elif "Warning" in line or "WARNING" in line:
                warnings.append(line.strip())
        
        error_msg = f"tleap failed to generate complex topology\n"
        error_msg += f"Log file: {log_file}\n\n"
        
        if fatal_errors:
            error_msg += "Errors:\n"
            for err in fatal_errors[:10]:  # Show first 10
                error_msg += f"  {err}\n"
            if len(fatal_errors) > 10:
                error_msg += f"  ... and {len(fatal_errors) - 10} more errors\n"
        
        raise ComplexAssemblyError(error_msg)
    
    # Check for warnings
    warnings = [line for line in result.stdout.splitlines() if "Warning" in line]
    if warnings:
        logger.warning(f"tleap warnings ({len(warnings)}):")
        for warn in warnings[:5]:
            logger.warning(f"  {warn.strip()}")
        if len(warnings) > 5:
            logger.warning(f"  ... and {len(warnings) - 5} more warnings")
        logger.warning(f"Check log file: {log_file}")
    
    logger.success("Complex assembly complete")
    logger.info(f"  Topology:    {complex_prmtop}")
    logger.info(f"  Coordinates: {complex_inpcrd}")


def prepare_complex(
    protein_dir: Path,
    ligand_dir: Path,
    output_dir: Path,
) -> None:
    """
    Prepare protein-ligand complex for MM/GBSA calculation.
    
    This is a strict, single-frame complex assembly workflow:
    1. Validate all input files exist
    2. Assemble complex using tleap (protein + ligand)
    
    Ligand coordinates are taken directly from the parameterized ligand.mol2.
    No solvation, minimization, or MD is performed.
    
    Args:
        protein_dir: Directory containing protein.prmtop and protein.inpcrd
        ligand_dir: Directory containing ligand.mol2 and ligand.frcmod
        output_dir: Directory for output files
        
    Raises:
        ComplexAssemblyError: If any step fails
    """
    logger.info("Preparing protein-ligand complex for MM/GBSA")
    logger.info(f"  Protein: {protein_dir}")
    logger.info(f"  Ligand:  {ligand_dir}")
    logger.info(f"  Output:  {output_dir}")
    
    # Step 1: Validate inputs (simplified - no complex_pdb needed)
    protein_prmtop = protein_dir / "protein.prmtop"
    protein_inpcrd = protein_dir / "protein.inpcrd"
    ligand_mol2 = ligand_dir / "ligand.mol2"
    ligand_frcmod = ligand_dir / "ligand.frcmod"
    
    if not protein_prmtop.exists():
        raise ComplexAssemblyError(
            f"Missing protein topology: {protein_prmtop}\n"
            f"Run: mmgbsa prep-protein <receptor.pdb> -o {protein_dir}"
        )
    
    if not protein_inpcrd.exists():
        raise ComplexAssemblyError(
            f"Missing protein coordinates: {protein_inpcrd}\n"
            f"Run: mmgbsa prep-protein <receptor.pdb> -o {protein_dir}"
        )
    
    if not ligand_mol2.exists():
        raise ComplexAssemblyError(
            f"Missing ligand structure: {ligand_mol2}\n"
            f"Run: mmgbsa parameterize <ligand.mol2> -o {ligand_dir}"
        )
    
    if not ligand_frcmod.exists():
        raise ComplexAssemblyError(
            f"Missing ligand parameters: {ligand_frcmod}\n"
            f"Run: mmgbsa parameterize <ligand.mol2> -o {ligand_dir}"
        )
    
    logger.success("All input files validated")
    
    # Step 2: Assemble complex with tleap
    output_dir.mkdir(parents=True, exist_ok=True)
    
    assemble_complex(
        protein_prmtop,
        protein_inpcrd,
        ligand_mol2,
        ligand_frcmod,
        output_dir,
    )
    
    logger.success("Complex preparation complete")
    logger.info(f"\nOutput files in {output_dir}:")
    logger.info(f"  - complex.prmtop")
    logger.info(f"  - complex.inpcrd")
    logger.info(f"  - tleap_complex.in (tleap input script)")
    logger.info(f"  - tleap_complex.log (tleap output)")
