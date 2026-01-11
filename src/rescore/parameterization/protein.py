"""
Protein preparation module for MM/GBSA calculations.

This module handles strict, reproducible protein preparation using AmberTools tleap
with the ff14SB force field. The philosophy is fail-fast and explicit:

- Accept only protein-only PDB files (no ligands, metals, waters, ions)
- Validate all residues are standard amino acids
- Do NOT infer protonation states
- Do NOT add or remove atoms
- Let tleap fail loudly if structure has issues

This ensures reproducibility and prevents silent errors that reviewers would catch.
"""

import subprocess
from pathlib import Path
from typing import Optional, Set

from Bio.PDB import PDBParser
from loguru import logger

from rescore.config import PROTEIN_FF, STANDARD_AMINO_ACIDS
from rescore.parameterization.errors import ParameterizationError


class ProteinPreparationError(ParameterizationError):
    """Raised when protein preparation fails."""
    pass


def run_pdb4amber(
    input_pdb: Path,
    output_pdb: Path,
) -> None:
    """
    Run pdb4amber to normalize PDB file for Amber compatibility.
    
    pdb4amber performs critical preprocessing:
    - Converts non-standard atom names (HN → H)
    - Removes alternate conformations
    - Fixes residue naming
    - Handles terminal residues
    - Removes waters and ions
    
    Args:
        input_pdb: Path to original PDB file
        output_pdb: Path where normalized PDB will be written
        
    Raises:
        ProteinPreparationError: If pdb4amber fails
    """
    logger.info("Running pdb4amber to normalize PDB for Amber compatibility")
    logger.info(f"  Input:  {input_pdb}")
    logger.info(f"  Output: {output_pdb}")
    
    cmd = [
        "pdb4amber",
        "-i", str(input_pdb.absolute()),
        "-o", str(output_pdb.absolute()),
        "-y",  # Strip hydrogens (tleap will add them back with correct naming)
    ]
    
    logger.debug(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,  # Don't raise on non-zero exit
        )
        
        # Log output
        if result.stdout:
            logger.debug(f"pdb4amber stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"pdb4amber stderr:\n{result.stderr}")
        
        # Check if output was created
        if not output_pdb.exists():
            raise ProteinPreparationError(
                f"pdb4amber failed to generate output file: {output_pdb}\n"
                f"Exit code: {result.returncode}\n\n"
                f"Output:\n{result.stdout}\n\n"
                f"Errors:\n{result.stderr}"
            )
        
        # Check file size
        if output_pdb.stat().st_size == 0:
            raise ProteinPreparationError(
                f"pdb4amber generated empty output file: {output_pdb}\n"
                f"This usually indicates the input PDB is severely malformed."
            )
        
        logger.success(f"pdb4amber completed: {output_pdb}")
        
    except FileNotFoundError:
        raise ProteinPreparationError(
            "pdb4amber not found in PATH.\n\n"
            "pdb4amber is part of AmberTools. Install with:\n"
            "  conda install -c conda-forge ambertools\n"
        )
    except Exception as e:
        if isinstance(e, ProteinPreparationError):
            raise
        raise ProteinPreparationError(f"Unexpected error running pdb4amber: {e}")


def validate_protein_only(pdb_path: Path) -> None:
    """
    Validate that PDB contains only standard amino acid residues.
    
    Rejects:
    - Ligands (non-standard residues)
    - Metals
    - Waters (HOH, WAT)
    - Ions (NA, CL, etc.)
    - Modified residues
    
    Accepts:
    - Standard amino acids (20 common + protonation variants)
    - Terminal variants (NXXX, CXXX)
    
    Args:
        pdb_path: Path to protein PDB file
        
    Raises:
        ProteinPreparationError: If non-protein residues found
    """
    logger.info(f"Validating protein-only structure: {pdb_path}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    
    non_protein_residues = []
    residue_counts = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                
                # Skip hetero flags for waters/ions and histidine variants
                # H_HSE, H_HSD, H_HIP are valid for protonated histidines
                if residue.id[0] not in [' ', 'H_MSE', 'H_M3L', 'H_CAS', 'H_HSE', 'H_HSD', 'H_HIP', 'H_HIE', 'H_HID']:
                    # Check if it's a water
                    if resname in ['HOH', 'WAT', 'TIP', 'TIP3', 'SOL']:
                        non_protein_residues.append(
                            f"Water: {resname} (chain {chain.id}, residue {residue.id[1]})"
                        )
                        continue
                    
                    # Check if it's a common ion
                    if resname in ['NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'FE', 'CU', 'MN']:
                        non_protein_residues.append(
                            f"Ion: {resname} (chain {chain.id}, residue {residue.id[1]})"
                        )
                        continue
                    
                    # Other HETATM - likely ligand or modified residue
                    non_protein_residues.append(
                        f"Non-standard: {resname} (chain {chain.id}, residue {residue.id[1]})"
                    )
                    continue
                
                # Check if residue is a standard amino acid
                if resname not in STANDARD_AMINO_ACIDS:
                    # Check if it's a terminal variant (N-terminal or C-terminal)
                    if not (resname.startswith('N') or resname.startswith('C')):
                        non_protein_residues.append(
                            f"Unknown residue: {resname} (chain {chain.id}, residue {residue.id[1]})"
                        )
                        continue
                
                # Count residues
                residue_counts[resname] = residue_counts.get(resname, 0) + 1
    
    if non_protein_residues:
        error_lines = [
            "PDB file contains non-protein residues:",
            "",
            "Found issues:",
        ]
        
        # Group by type
        waters = [r for r in non_protein_residues if r.startswith("Water:")]
        ions = [r for r in non_protein_residues if r.startswith("Ion:")]
        nonstandard = [r for r in non_protein_residues if r.startswith("Non-standard:")]
        unknown = [r for r in non_protein_residues if r.startswith("Unknown:")]
        
        if waters:
            error_lines.append(f"  Waters ({len(waters)}):")
            for w in waters[:5]:
                error_lines.append(f"    - {w}")
            if len(waters) > 5:
                error_lines.append(f"    ... and {len(waters) - 5} more")
        
        if ions:
            error_lines.append(f"  Ions ({len(ions)}):")
            for i in ions[:5]:
                error_lines.append(f"    - {i}")
            if len(ions) > 5:
                error_lines.append(f"    ... and {len(ions) - 5} more")
        
        if nonstandard:
            error_lines.append(f"  Non-standard residues ({len(nonstandard)}):")
            for n in nonstandard[:5]:
                error_lines.append(f"    - {n}")
            if len(nonstandard) > 5:
                error_lines.append(f"    ... and {len(nonstandard) - 5} more")
        
        if unknown:
            error_lines.append(f"  Unknown residues ({len(unknown)}):")
            for u in unknown:
                error_lines.append(f"    - {u}")
        
        error_lines.extend([
            "",
            "Action required:",
            "  1. Remove ligands, waters, and ions from PDB",
            "  2. Use a structure editor (PyMOL, Chimera) to clean the structure",
            "  3. Save protein-only PDB with command like:",
            "     PyMOL> select protein, polymer.protein",
            "     PyMOL> save protein_only.pdb, protein",
            "",
            "For ligand parameterization, use: rescore parameterize",
        ])
        
        raise ProteinPreparationError("\n".join(error_lines))
    
    logger.success(f"Protein validation passed")
    logger.info(f"  Total residues: {sum(residue_counts.values())}")
    
    # Show protonation variants if present
    protonation_variants = {k: v for k, v in residue_counts.items() 
                           if k in ['HIS', 'HID', 'HIE', 'HIP', 'HSE', 'HSD', 'HSP']}
    if protonation_variants:
        logger.info(f"  Histidine variants: {protonation_variants}")


def prepare_protein(
    pdb_path: Path,
    output_dir: Path,
    force_field: str = PROTEIN_FF,
    skip_pdb4amber: bool = False,
) -> tuple[Path, Path]:
    """
    Prepare protein topology and coordinates using AmberTools tleap.
    
    Workflow:
    1. Run pdb4amber to normalize PDB (unless --skip-pdb4amber)
    2. Validate protein-only structure (no ligands, waters, ions)
    3. Generate tleap input script
    4. Run tleap with ff14SB
    5. Parse output for errors
    6. Return paths to prmtop and inpcrd files
    
    pdb4amber preprocessing (recommended):
    - Converts non-standard atom names (HN → H)
    - Removes alternate conformations
    - Fixes residue naming
    - Handles terminal residues
    - Removes waters and ions
    
    This function does NOT:
    - Infer protonation states (user must provide correct PDB)
    - Add missing heavy atoms (tleap will fail with clear error)
    - Minimize structure
    
    Args:
        pdb_path: Path to protein PDB file
        output_dir: Directory for output files
        force_field: Amber force field (default: ff14SB)
        skip_pdb4amber: If True, skip pdb4amber preprocessing (not recommended)
        
    Returns:
        Tuple of (prmtop_path, inpcrd_path)
        
    Raises:
        ProteinPreparationError: If validation or tleap fails
    """
    pdb_path = Path(pdb_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not pdb_path.exists():
        raise ProteinPreparationError(f"PDB file not found: {pdb_path}")
    
    logger.info(f"Preparing protein: {pdb_path}")
    logger.info(f"  Force field: {force_field}")
    logger.info(f"  Output directory: {output_dir}")
    
    # Step 1: Run pdb4amber preprocessing (unless skipped)
    if not skip_pdb4amber:
        # Generate amber-normalized filename
        amber_pdb = output_dir / f"{pdb_path.stem}_amber.pdb"
        
        run_pdb4amber(pdb_path, amber_pdb)
        
        # Use the amber-processed PDB for subsequent steps
        working_pdb = amber_pdb
        logger.info(f"Using pdb4amber-processed file: {working_pdb}")
    else:
        logger.warning("Skipping pdb4amber preprocessing (--skip-pdb4amber enabled)")
        logger.warning("This may cause issues with non-standard atom naming")
        
        # Copy original to output dir
        import shutil
        working_pdb = output_dir / "protein_input.pdb"
        shutil.copy(pdb_path, working_pdb)
        logger.info(f"Using original PDB file: {working_pdb}")
    
    # Step 2: Validate protein-only structure
    try:
        validate_protein_only(working_pdb)
    except ProteinPreparationError as e:
        if skip_pdb4amber:
            raise ProteinPreparationError(
                f"Validation failed on original PDB:\n\n{e}\n\n"
                f"TIP: Try without --skip-pdb4amber to let pdb4amber clean the structure."
            )
        else:
            raise ProteinPreparationError(
                f"Validation failed even after pdb4amber preprocessing:\n\n{e}\n\n"
                f"This indicates structural issues beyond what pdb4amber can fix.\n"
                f"Manual structure preparation required."
            )
    
    # Step 3: Generate tleap input script
    tleap_in = output_dir / "tleap_protein.in"
    prmtop_out = output_dir / "protein.prmtop"
    inpcrd_out = output_dir / "protein.inpcrd"
    
    tleap_script = f"""# tleap script for protein preparation
# Force field: {force_field}
# Generated by rescore
# Input: {working_pdb.name}

# Load force field
source leaprc.protein.{force_field}

# Load protein PDB
protein = loadpdb {working_pdb.name}

# Check for issues
check protein

# Save topology and coordinates
saveamberparm protein protein.prmtop protein.inpcrd

# Quit
quit
"""
    
    with open(tleap_in, 'w') as f:
        f.write(tleap_script)
    
    logger.debug(f"Generated tleap input: {tleap_in}")
    
    # Step 4: Run tleap
    logger.info(f"Running tleap on {working_pdb.name}...")
    
    tleap_log = output_dir / "tleap_protein.log"
    
    try:
        result = subprocess.run(
            ["tleap", "-f", "tleap_protein.in"],
            capture_output=True,
            text=True,
            cwd=str(output_dir.absolute()),
        )
        
        # Save log
        with open(tleap_log, 'w') as f:
            f.write("=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n\n=== STDERR ===\n")
            f.write(result.stderr)
        
        logger.debug(f"tleap log saved: {tleap_log}")
        
        # Check for errors in output
        output = result.stdout + result.stderr
        
        # Parse for common errors
        errors = []
        warnings = []
        
        for line in output.split('\n'):
            line_lower = line.lower()
            
            if 'error' in line_lower or 'fatal' in line_lower:
                errors.append(line.strip())
            elif 'warning' in line_lower:
                warnings.append(line.strip())
            elif 'missing' in line_lower and 'atom' in line_lower:
                errors.append(line.strip())
            elif 'unknown residue' in line_lower:
                errors.append(line.strip())
        
        # Log warnings
        if warnings:
            logger.warning(f"tleap warnings ({len(warnings)}):")
            for w in warnings[:5]:
                logger.warning(f"  {w}")
            if len(warnings) > 5:
                logger.warning(f"  ... and {len(warnings) - 5} more warnings")
        
        # Check if output files were created
        if not prmtop_out.exists() or not inpcrd_out.exists():
            error_msg = [
                "tleap failed to generate output files",
                "",
            ]
            
            # Note whether pdb4amber was used
            if skip_pdb4amber:
                error_msg.extend([
                    "NOTE: pdb4amber preprocessing was SKIPPED (--skip-pdb4amber)",
                    "This may have caused naming issues.",
                    "",
                ])
            else:
                error_msg.extend([
                    f"NOTE: pdb4amber preprocessing WAS USED on {pdb_path.name}",
                    f"Processed file: {working_pdb.name}",
                    "Remaining errors indicate structural issues beyond naming.",
                    "",
                ])
            
            # Check for specific common errors
            if "HN" in output and "does not have a type" in output:
                if skip_pdb4amber:
                    error_msg.extend([
                        "DETECTED: Non-standard atom naming (HN instead of H)",
                        "",
                        "The PDB file uses 'HN' for backbone amide hydrogens.",
                        "Amber force fields expect 'H' for these atoms.",
                        "",
                        "FIX: Remove --skip-pdb4amber flag to let pdb4amber fix this automatically.",
                        "",
                    ])
                else:
                    error_msg.extend([
                        "UNEXPECTED: Non-standard atom naming still present after pdb4amber",
                        "",
                        "pdb4amber should have fixed 'HN' → 'H' conversion.",
                        "This may indicate a pdb4amber version issue or unusual PDB format.",
                        "",
                        f"Check the processed file: {working_pdb}",
                        "",
                    ])
            elif "missing" in output.lower() or "does not have a type" in output:
                error_msg.extend([
                    "DETECTED: Missing or incorrectly named atoms",
                    "",
                    "Common causes:",
                    "  - Incomplete residues (missing backbone or sidechain atoms)",
                    "  - Non-standard residues that pdb4amber cannot handle",
                    "  - Severely malformed PDB structure",
                    "",
                ])
            
            if errors:
                error_msg.append("tleap errors (first 15):")
                for e in errors[:15]:
                    error_msg.append(f"  {e}")
                if len(errors) > 15:
                    error_msg.append(f"  ... and {len(errors) - 15} more errors")
                error_msg.append("")
            
            error_msg.extend([
                f"Check tleap log for full details: {tleap_log}",
                "",
                "Next steps:",
                "  1. Inspect the processed PDB file manually",
                "  2. Use structure editor (PyMOL, Chimera) to fix issues",
                "  3. Ensure complete residues with standard naming",
            ])
            
            raise ProteinPreparationError("\n".join(error_msg))
        
        # Check if tleap reported errors even if files exist
        if errors:
            logger.warning("tleap completed but reported errors:")
            for e in errors[:10]:
                logger.warning(f"  {e}")
            logger.warning(f"Check log file: {tleap_log}")
        
        logger.success("Protein preparation complete")
        logger.info(f"  Topology:    {prmtop_out}")
        logger.info(f"  Coordinates: {inpcrd_out}")
        
        return prmtop_out, inpcrd_out
        
    except subprocess.CalledProcessError as e:
        raise ProteinPreparationError(
            f"tleap execution failed with exit code {e.returncode}\n"
            f"Check log file: {tleap_log}"
        )
    except Exception as e:
        raise ProteinPreparationError(f"Unexpected error during tleap execution: {e}")
