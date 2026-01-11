"""
MM/GBSA single-frame rescoring calculations.

This module implements MINIMAL single-frame MM/GBSA rescoring using MMPBSA.py.

CRITICAL SCOPE LIMITATIONS:
- Single structure only (no trajectory, no MD)
- Implicit solvent (GB or PB) only
- Optional restrained minimization before scoring
- No entropy calculation
- No per-residue decomposition

This is POST-DOCKING RESCORING, not thermodynamically rigorous binding free energy.
Use for relative ranking only, not absolute ΔG values.
"""

import json
import re
import subprocess
from pathlib import Path
from typing import Dict, Optional

from loguru import logger

from rescore.parameterization.errors import ParameterizationError


class MMGBSACalculationError(ParameterizationError):
    """Raised when MM/GBSA calculation fails."""
    pass


def write_metadata(output_dir: Path, method: str) -> Path:
    """
    Write metadata file describing the calculation parameters.
    
    Args:
        output_dir: Directory where metadata will be written
        method: Solvation method ("gb" or "pb")
        
    Returns:
        Path to metadata file
    """
    metadata = {
        "method": method.lower(),
        "mode": "single-frame",
        "purpose": "post-docking rescoring",
        "entropy": False
    }
    
    prefix = "mmpbsa" if method.lower() == "pb" else "rescore"
    metadata_file = output_dir / f"{prefix}_metadata.json"
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)
    
    logger.debug(f"Metadata written: {metadata_file}")
    return metadata_file


def validate_complex_inputs(complex_dir: Path) -> tuple[Path, Path]:
    """
    Validate that complex topology and coordinates exist.
    
    Args:
        complex_dir: Directory containing complex parameters
        
    Returns:
        Tuple of (complex_prmtop, complex_inpcrd)
        
    Raises:
        MMGBSACalculationError: If required files are missing
    """
    logger.info("Validating complex input files...")
    
    complex_prmtop = complex_dir / "complex.prmtop"
    complex_inpcrd = complex_dir / "complex.inpcrd"
    
    if not complex_prmtop.exists():
        raise MMGBSACalculationError(
            f"Missing complex topology: {complex_prmtop}\n"
            f"Run: rescore assemble --protein <protein_dir> --ligand <ligand_dir> -o {complex_dir}"
        )
    
    if not complex_inpcrd.exists():
        raise MMGBSACalculationError(
            f"Missing complex coordinates: {complex_inpcrd}\n"
            f"Run: rescore assemble --protein <protein_dir> --ligand <ligand_dir> -o {complex_dir}"
        )
    
    logger.success("Complex input files validated")
    return complex_prmtop, complex_inpcrd


def minimize_complex(
    complex_prmtop: Path,
    complex_inpcrd: Path,
    output_dir: Path,
    method: str = "gb",
    restraint_weight: float = 2.0,
    maxcyc: int = 100,
) -> Path:
    """
    Run restrained energy minimization on assembled complex.
    
    This prepares the complex for MM/GBSA by relaxing steric clashes
    and electrostatics while preserving the overall protein structure.
    
    Protocol:
    - Positional restraints on protein backbone (N, CA, C) atoms
    - Ligand and side chains free to relax
    - Steepest descent (50%) + conjugate gradient (50%)
    - Implicit solvent (GB) with reasonable cutoff
    
    IMPORTANT: This is NOT molecular dynamics. This is minimization only
    to remove bad contacts from docking. Still single-frame rescoring.
    
    Args:
        complex_prmtop: Complex topology file
        complex_inpcrd: Complex coordinate file
        output_dir: Directory for minimization outputs
        method: Solvation method ("gb" or "pb") - always uses GB for speed
        restraint_weight: Restraint force constant (kcal/mol/Å²)
        maxcyc: Maximum minimization cycles (default: 100)
        
    Returns:
        Path to minimized coordinates (minimized.inpcrd)
        
    Raises:
        MMGBSACalculationError: If sander fails
    """
    logger.info("Running restrained minimization (Prime-style MM/GBSA preparation)")
    logger.info(f"  Restraining protein backbone (N, CA, C)")
    logger.info(f"  Ligand and side chains free to relax")
    logger.info(f"  Restraint weight: {restraint_weight} kcal/mol/Å²")
    logger.info(f"  Max cycles: {maxcyc}")
    logger.info(f"  Solvation: GB implicit solvent")
    
    if method.lower() == "pb":
        logger.info("  Note: Using GB for minimization speed; MM/GBSA will use PB")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output files
    minimized_inpcrd = output_dir / "minimized.inpcrd"
    min_log = output_dir / "minimization.log"
    min_info = output_dir / "minimization_info.txt"
    min_input = output_dir / "minimization.in"
    
    # Generate sander input for restrained minimization
    ncyc = maxcyc // 2  # First half steepest descent
    
    min_input_content = f"""Restrained minimization for MM/GBSA preparation
&cntrl
  imin=1,
  maxcyc={maxcyc},
  ncyc={ncyc},
  ntb=0,
  igb=5,
  saltcon=0.15,
  cut=12.0,
  ntr=1,
  restraintmask='@N,CA,C',
  restraint_wt={restraint_weight},
  ntpr=50,
/
"""
    
    with open(min_input, "w") as f:
        f.write(min_input_content)
    
    logger.debug(f"Minimization input written: {min_input}")
    
    # Run sander
    logger.info("Running sander minimization...")
    logger.debug(f"Command: sander -O -i {min_input} -p {complex_prmtop} -c {complex_inpcrd} -o {min_info} -r {minimized_inpcrd} -ref {complex_inpcrd}")
    
    # Run without capturing output first to see if there's an immediate error
    result = subprocess.run(
        [
            "sander",
            "-O",
            "-i", str(min_input),
            "-p", str(complex_prmtop),
            "-c", str(complex_inpcrd),
            "-o", str(min_info),
            "-r", str(minimized_inpcrd),
            "-ref", str(complex_inpcrd),
        ],
        capture_output=True,
        text=True,
        timeout=120,  # 2 minute timeout
    )
    
    # Write output to log
    with open(min_log, "w") as f:
        f.write(result.stdout)
        if result.stderr:
            f.write("\n\n=== STDERR ===\n")
            f.write(result.stderr)
    
    if result.returncode != 0 or not minimized_inpcrd.exists():
        with open(min_log, "r") as f:
            log_content = f.read()
        raise MMGBSACalculationError(
            f"sander minimization failed (exit code {result.returncode})\n"
            f"Check log: {min_log}\n"
            f"Last 50 lines:\n{log_content[-2000:]}"
        )
    
    logger.success(f"Minimization complete: {minimized_inpcrd}")
    logger.info(f"  Log: {min_log}")
    logger.info(f"  Details: {min_info}")
    
    # Convert minimized structure to PDB for visualization
    minimized_pdb = output_dir / "minimized_complex.pdb"
    pdb_debug_log = output_dir / "pdb_generation_debug.log"
    logger.info("Converting minimized structure to PDB format...")
    
    try:
        # Write debug info to file (visible even when logger is suppressed)
        with open(pdb_debug_log, "w") as debug:
            debug.write(f"PDB Generation Debug Log\n")
            debug.write(f"========================\n")
            debug.write(f"complex_prmtop: {complex_prmtop}\n")
            debug.write(f"complex_prmtop exists: {complex_prmtop.exists()}\n")
            debug.write(f"minimized_inpcrd: {minimized_inpcrd}\n")
            debug.write(f"minimized_inpcrd exists: {minimized_inpcrd.exists()}\n")
            debug.write(f"output_dir: {output_dir}\n")
            debug.write(f"minimized_pdb target: {minimized_pdb}\n\n")
        
        # ambpdb has issues with spaces in absolute paths, use relative paths
        # We're in output_dir (minimization/), complex_prmtop is in ../../complex/
        rel_prmtop = Path("../../complex") / complex_prmtop.name
        rel_inpcrd = minimized_inpcrd.name  # Already in current dir
        
        pdb_result = subprocess.run(
            ["ambpdb", "-p", str(rel_prmtop), "-c", str(rel_inpcrd)],
            capture_output=True,
            text=True,
            check=False,  # Don't raise exception on non-zero exit
            cwd=str(output_dir.absolute()),  # Run from output directory
        )
        
        # Append result to debug log
        with open(pdb_debug_log, "a") as debug:
            debug.write(f"ambpdb command: ambpdb -p {complex_prmtop} -c {minimized_inpcrd}\n")
            debug.write(f"returncode: {pdb_result.returncode}\n")
            debug.write(f"stdout length: {len(pdb_result.stdout) if pdb_result.stdout else 0}\n")
            debug.write(f"stdout empty: {not pdb_result.stdout}\n")
            debug.write(f"stderr length: {len(pdb_result.stderr) if pdb_result.stderr else 0}\n")
            debug.write(f"stderr: {pdb_result.stderr}\n")
            if pdb_result.stdout:
                debug.write(f"\nFirst 500 chars of stdout:\n{pdb_result.stdout[:500]}\n")
        
        logger.debug(f"ambpdb returncode: {pdb_result.returncode}")
        logger.debug(f"ambpdb stdout length: {len(pdb_result.stdout) if pdb_result.stdout else 0}")
        logger.debug(f"ambpdb stderr: {pdb_result.stderr[:200] if pdb_result.stderr else 'None'}")
        
        if pdb_result.returncode == 0 and pdb_result.stdout:
            with open(minimized_pdb, "w") as f:
                f.write(pdb_result.stdout)
            logger.success(f"Minimized PDB saved: {minimized_pdb} ({len(pdb_result.stdout)} bytes)")
        else:
            logger.warning(f"Could not generate PDB file")
            logger.warning(f"  ambpdb returncode: {pdb_result.returncode}")
            logger.warning(f"  stdout empty: {not pdb_result.stdout}")
            logger.warning(f"  stderr: {pdb_result.stderr[:500] if pdb_result.stderr else 'None'}")
            # Create empty file so user knows it was attempted
            minimized_pdb.touch()
    except Exception as e:
        logger.error(f"Exception during PDB generation: {e}")
        # Log exception to debug file
        with open(pdb_debug_log, "a") as debug:
            debug.write(f"\nEXCEPTION: {e}\n")
            import traceback
            debug.write(traceback.format_exc())
        # Create empty file to indicate attempt was made
        minimized_pdb.touch()
    
    return minimized_inpcrd


def generate_mmpbsa_input(output_dir: Path, method: str = "gb") -> Path:
    """
    Generate MMPBSA.py input file for single-frame calculation.
    
    This creates a minimal input for post-docking rescoring:
    - Single snapshot (startframe=1, endframe=1, interval=1)
    - GB or PB implicit solvent
    - No entropy calculation
    - No decomposition
    
    Args:
        output_dir: Directory where input file will be written
        method: Solvation method - "gb" (default) or "pb"
        
    Returns:
        Path to generated input file
        
    Raises:
        ValueError: If method is not "gb" or "pb"
    """
    method = method.lower()
    if method not in ["gb", "pb"]:
        raise ValueError(f"Invalid method: {method}. Must be 'gb' or 'pb'")
    
    logger.info(f"Generating MMPBSA.py input file (method={method.upper()})...")
    
    if method == "pb":
        logger.warning(
            "PB is slower and more sensitive to structure quality. "
            "In single-frame mode, PB does NOT make results thermodynamically rigorous."
        )
    
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = "mmpbsa" if method == "pb" else "rescore"
    input_file = output_dir / f"{prefix}.in"
    
    # Common header
    header = """# Single-frame MM/GBSA rescoring
# Generated by rescore.calculation.rescore
#
# SCOPE: Post-docking rescoring for relative ranking
# NOT thermodynamically rigorous binding free energy
#
# Method: MM/GBSA with {} implicit solvent
# Frames: Single structure (no MD trajectory)
# Entropy: NONE (too expensive, not reliable for single frame)
""".format("GB (igb=5)" if method == "gb" else "PB")
    
    # General section (same for both methods)
    general_section = """&general
    startframe=1,
    endframe=1,
    interval=1,
    verbose=2,
    keep_files=0,
/
"""
    
    # Method-specific section
    if method == "gb":
        method_section = """&gb
    igb=5,
    saltcon=0.15,
/
"""
    else:  # pb
        method_section = """&pb
    istrng=0.15,
/
"""
    
    input_content = header + "\n" + general_section + "\n" + method_section
    
    with open(input_file, "w") as f:
        f.write(input_content)
    
    logger.success(f"MM/GBSA input generated: {input_file}")
    return input_file


def generate_ligand_topology(
    ligand_mol2: Path,
    ligand_frcmod: Path,
    output_prmtop: Path,
    output_inpcrd: Path,
    output_dir: Path,
) -> None:
    """
    Generate ligand topology using tleap.
    
    Args:
        ligand_mol2: Ligand MOL2 file
        ligand_frcmod: Ligand force field modifications
        output_prmtop: Output topology file
        output_inpcrd: Output coordinates file
        output_dir: Working directory
        
    Raises:
        MMGBSACalculationError: If tleap fails
    """
    # Generate tleap input
    tleap_input = output_dir / "tleap_ligand.in"
    
    # Use relative paths
    rel_frcmod = Path("..") / ligand_frcmod.relative_to(output_dir.parent)
    rel_mol2 = Path("..") / ligand_mol2.relative_to(output_dir.parent)
    
    tleap_script = f"""# Ligand topology for MM/GBSA
source leaprc.gaff2
loadamberparams {rel_frcmod}
ligand = loadmol2 {rel_mol2}
check ligand
saveamberparm ligand ligand.prmtop ligand.inpcrd
quit
"""
    
    with open(tleap_input, "w") as f:
        f.write(tleap_script)
    
    # Run tleap
    result = subprocess.run(
        ["tleap", "-f", "tleap_ligand.in"],
        cwd=output_dir,
        capture_output=True,
        text=True,
    )
    
    if not output_prmtop.exists():
        raise MMGBSACalculationError(
            f"Failed to generate ligand topology\n"
            f"tleap output: {result.stdout}\n"
            f"tleap errors: {result.stderr}"
        )
    
    logger.success(f"Ligand topology generated: {output_prmtop.name}")


def run_mmpbsa_calculation(
    complex_prmtop: Path,
    complex_inpcrd: Path,
    protein_prmtop: Path,
    ligand_prmtop: Path,
    input_file: Path,
    output_dir: Path,
    method: str = "gb",
) -> None:
    """
    Run MMPBSA.py for single-frame MM/GBSA calculation.
    
    This runs AmberTools MMPBSA.py with:
    - Complex, receptor, and ligand topologies
    - Single structure input (no trajectory)
    - GB implicit solvent
    
    Args:
        complex_prmtop: Complex topology file
        complex_inpcrd: Complex coordinates (single frame)
        protein_prmtop: Receptor topology file
        ligand_prmtop: Ligand topology file
        input_file: MMPBSA.py input file
        output_dir: Directory for output files
        
    Raises:
        MMGBSACalculationError: If MMPBSA.py fails
    """
    logger.info("Running MMPBSA.py calculation...")
    logger.info(f"  Complex topology: {complex_prmtop.name}")
    logger.info(f"  Receptor topology: {protein_prmtop.name}")
    logger.info(f"  Ligand topology: {ligand_prmtop.name}")
    logger.info(f"  Complex coords:   {complex_inpcrd.name}")
    logger.info(f"  Input file:       {input_file.name}")
    
    # MMPBSA.py command for single-frame calculation
    # Run from parent directory to avoid path issues with spaces
    # Use relative paths
    
    parent_dir = complex_prmtop.parent.parent  # Go up from complex_params/ to workspace root
    rel_input = input_file.relative_to(parent_dir)
    rel_complex_prmtop = complex_prmtop.relative_to(parent_dir)
    rel_complex_inpcrd = complex_inpcrd.relative_to(parent_dir)
    rel_protein_prmtop = protein_prmtop.relative_to(parent_dir)
    rel_ligand_prmtop = ligand_prmtop.relative_to(parent_dir)
    rel_output = output_dir.relative_to(parent_dir)
    
    prefix = "mmpbsa" if method.lower() == "pb" else "rescore"
    
    cmd = [
        "MMPBSA.py",
        "-O",
        "-i", str(rel_input),
        "-cp", str(rel_complex_prmtop),
        "-rp", str(rel_protein_prmtop),
        "-lp", str(rel_ligand_prmtop),
        "-y", str(rel_complex_inpcrd),
        "-o", str(rel_output / f"{prefix}_output.dat"),
    ]
    
    logger.debug(f"Command: {' '.join(cmd)}")
    logger.debug(f"Working directory: {parent_dir}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=parent_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        
        # Save output
        prefix = "mmpbsa" if method.lower() == "pb" else "rescore"
        log_file = output_dir / f"{prefix}.log"
        with open(log_file, "w") as f:
            f.write("=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)
        
        logger.debug(f"MMPBSA.py log saved: {log_file}")
        
        # Check for output files - use the output file we specified
        output_file = output_dir / f"{prefix}_output.dat"
        
        if not output_file.exists():
            # Parse errors from output
            error_lines = []
            for line in result.stderr.splitlines():
                if "Error" in line or "ERROR" in line or "Fatal" in line:
                    error_lines.append(line.strip())
            
            error_msg = "MMPBSA.py failed to generate results\n"
            error_msg += f"Log file: {log_file}\n\n"
            
            if error_lines:
                error_msg += "Errors:\n"
                for err in error_lines[:10]:
                    error_msg += f"  {err}\n"
            else:
                error_msg += "Check log file for details\n"
            
            raise MMGBSACalculationError(error_msg)
        
        logger.success("MMPBSA.py calculation complete")
        logger.info(f"  Results: {output_file}")
        
    except FileNotFoundError:
        raise MMGBSACalculationError(
            "MMPBSA.py command not found\n"
            "Ensure AmberTools is properly installed and activated"
        )


def parse_rescore_results(output_dir: Path, method: str = "gb") -> Dict[str, float]:
    """
    Parse MMPBSA.py results and extract energy components.
    
    Extracts from FINAL_RESULTS_MMPBSA.dat:
    - ΔG_bind: Total binding free energy
    - ΔH: Enthalpy contribution
    - ΔG_GB: Polar solvation energy
    - ΔG_SA: Non-polar solvation energy
    
    Args:
        output_dir: Directory containing MMPBSA.py output
        method: Solvation method ("gb" or "pb")
        
    Returns:
        Dictionary of energy components (kcal/mol)
        
    Raises:
        MMGBSACalculationError: If results file cannot be parsed
    """
    logger.info("Parsing MM/GBSA results...")
    
    prefix = "mmpbsa" if method.lower() == "pb" else "rescore"
    results_file = output_dir / f"{prefix}_output.dat"
    
    if not results_file.exists():
        raise MMGBSACalculationError(
            f"Results file not found: {results_file}\n"
            "MMPBSA.py may have failed"
        )
    
    # Read results file
    with open(results_file, "r") as f:
        content = f.read()
    
    # Parse energy components using regex
    # Format: "DELTA Energy = value" or "DELTA Energy     value"
    energies = {}
    
    patterns = {
        "DELTA_G": r"DELTA TOTAL\s+(-?\d+\.\d+)",
        "DELTA_H": r"DELTA G gas\s+(-?\d+\.\d+)",
        "DELTA_G_GB": r"DELTA G solv\s+(-?\d+\.\d+)",
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            energies[key] = float(match.group(1))
        else:
            logger.debug(f"Could not parse {key} from results")
    
    if not energies:
        raise MMGBSACalculationError(
            f"Failed to parse energy values from {results_file}\n"
            "Results file may be corrupted or incomplete"
        )
    
    logger.success("Results parsed successfully")
    
    # Log key results
    if "DELTA_G" in energies:
        logger.info(f"  ΔG_bind = {energies['DELTA_G']:.2f} kcal/mol")
    
    return energies


def run_rescore(
    complex_dir: Path,
    output_dir: Path,
    method: str = "gb",
    minimize: bool = True,
) -> Dict[str, float]:
    """
    Run single-frame MM/GBSA rescoring calculation.
    
    This is a MINIMAL implementation for post-docking rescoring:
    - Single structure (no MD trajectory)
    - Optional restrained minimization (removes clashes)
    - GB or PB implicit solvent
    - No entropy calculation
    - No decomposition
    
    SCIENTIFIC POSITIONING:
    This calculates a MM/GBSA score for RELATIVE RANKING of compounds.
    It does NOT provide:
    - Thermodynamically rigorous binding free energies
    - Entropy contributions
    - Conformational sampling
    - Explicit solvent effects
    
    Use for prioritizing docking poses or compound series, not for
    reporting absolute ΔG values in publications.
    
    Args:
        complex_dir: Directory containing complex.prmtop and complex.inpcrd
        output_dir: Directory for calculation outputs
        method: Solvation method - "gb" (default) or "pb"
        minimize: Run restrained minimization before scoring (default: True)
        
    Returns:
        Dictionary of energy components (kcal/mol)
        
    Raises:
        MMGBSACalculationError: If calculation fails
    """
    logger.info("Starting MM/GBSA single-frame rescoring")
    logger.info(f"  Complex: {complex_dir}")
    logger.info(f"  Output:  {output_dir}")
    logger.info(f"  Method:  {method.upper()}")
    logger.info(f"  Minimize: {'YES' if minimize else 'NO'}")
    logger.warning("This is POST-DOCKING RESCORING, not thermodynamic ΔG")
    
    # Step 1: Validate inputs - need complex, protein, and ligand directories
    complex_prmtop, complex_inpcrd = validate_complex_inputs(complex_dir)
    
    # Check for protein and ligand topologies in parent directory
    # Supports two directory naming conventions:
    # 1. receptor_params/, ligand_params/, complex_params/ (from individual commands)
    # 2. protein/, ligand/, complex/ (from integrate command)
    parent_dir = complex_dir.parent
    
    # Try new naming convention first (protein/ligand/complex)
    protein_dir = parent_dir / "protein"
    ligand_dir = parent_dir / "ligand"
    
    # Fall back to old naming convention (receptor_params/ligand_params/complex_params)
    if not protein_dir.exists():
        protein_dir = parent_dir / "receptor_params"
        ligand_dir = parent_dir / "ligand_params"
    
    protein_prmtop = protein_dir / "protein.prmtop"
    ligand_mol2 = ligand_dir / "ligand.mol2"
    ligand_frcmod = ligand_dir / "ligand.frcmod"
    
    if not protein_prmtop.exists():
        raise MMGBSACalculationError(
            f"Missing receptor topology: {protein_prmtop}\n"
            f"MMPBSA.py requires separate receptor and ligand topologies\n"
            f"Expected directory structure:\n"
            f"  {parent_dir}/protein/protein.prmtop (or receptor_params/protein.prmtop)\n"
            f"  {parent_dir}/ligand/ligand.mol2 (or ligand_params/ligand.mol2)\n"
            f"  {parent_dir}/complex/complex.prmtop (or complex_params/complex.prmtop)"
        )
    
    if not ligand_mol2.exists():
        raise MMGBSACalculationError(
            f"Missing ligand structure: {ligand_mol2}\n"
            f"MMPBSA.py requires separate receptor and ligand topologies\n"
            f"Expected directory structure:\n"
            f"  {parent_dir}/protein/protein.prmtop (or receptor_params/protein.prmtop)\n"
            f"  {parent_dir}/ligand/ligand.mol2 (or ligand_params/ligand.mol2)\n"
            f"  {parent_dir}/complex/complex.prmtop (or complex_params/complex.prmtop)"
        )
    
    if not ligand_frcmod.exists():
        raise MMGBSACalculationError(
            f"Missing ligand parameters: {ligand_frcmod}\n"
            f"Run: rescore parameterize <ligand.mol2> -o {ligand_dir}"
        )
    
    logger.success("All required files found")
    
    # Step 2: Run minimization if requested
    if minimize:
        minimized_inpcrd = minimize_complex(
            complex_prmtop,
            complex_inpcrd,
            output_dir / "minimization",
            method=method,
        )
        # Use minimized coordinates for MM/GBSA
        complex_inpcrd = minimized_inpcrd
        logger.info("Using minimized coordinates for MM/GBSA calculation")
    else:
        logger.info("Skipping minimization (using original coordinates)")
    
    # Step 3: Generate MMPBSA.py input
    output_dir.mkdir(parents=True, exist_ok=True)
    input_file = generate_mmpbsa_input(output_dir, method=method)
    
    # Step 3.5: Write metadata
    write_metadata(output_dir, method)
    
    # Step 4: Generate ligand topology using tleap
    logger.info("Generating ligand topology with tleap...")
    ligand_prmtop = output_dir / "ligand.prmtop"
    ligand_inpcrd = output_dir / "ligand.inpcrd"
    
    generate_ligand_topology(ligand_mol2, ligand_frcmod, ligand_prmtop, ligand_inpcrd, output_dir)
    
    # Step 5: Run MMPBSA.py with protein and ligand topologies
    run_mmpbsa_calculation(
        complex_prmtop,
        complex_inpcrd,
        protein_prmtop,
        ligand_prmtop,
        input_file,
        output_dir,
        method=method,
    )
    
    # Step 6: Parse results
    energies = parse_rescore_results(output_dir, method=method)
    
    prefix = "mmpbsa" if method.lower() == "pb" else "rescore"
    logger.success(f"MM/{method.upper()} calculation complete")
    logger.info(f"\nOutput files in {output_dir}:")
    logger.info(f"  - {prefix}_output.dat (energy components)")
    logger.info(f"  - {prefix}.in (input file)")
    logger.info(f"  - {prefix}.log (detailed log)")
    
    return energies
