"""
Batch MM/GBSA rescoring for multiple ligands against a single receptor.

This module provides batch processing capabilities for MM/GBSA rescoring,
operating strictly as an orchestrator without modifying core chemistry or physics.

CRITICAL GUARDRAILS:
- Receptor topology is prepared once and reused (with hash validation)
- Only MOL2 ligands accepted (chemistry validation enforced)
- Results are for RELATIVE RANKING ONLY
- No modifications to parameterization, assembly, or calculation internals
"""

import csv
import hashlib
from pathlib import Path
from typing import List, Dict, Any

from loguru import logger
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.table import Table

from rescore.parameterization import parameterize_ligand
from rescore.parameterization.protein import prepare_protein
from rescore.parameterization.complex import prepare_complex
from rescore.calculation import run_rescore


class BatchError(Exception):
    """Base exception for batch processing errors."""
    pass


console = Console()


def calculate_file_hash(filepath: Path) -> str:
    """
    Calculate SHA256 hash of a file for integrity checking.
    
    Args:
        filepath: Path to file
        
    Returns:
        Hexadecimal hash string
    """
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def validate_receptor_topology_integrity(
    prmtop_path: Path,
    expected_hash: str,
) -> None:
    """
    Verify receptor topology has not been modified during batch run.
    
    This is a CRITICAL GUARDRAIL to ensure all ligands are compared
    against the same receptor structure.
    
    Args:
        prmtop_path: Path to protein.prmtop
        expected_hash: Expected SHA256 hash
        
    Raises:
        BatchError: If hash mismatch detected
    """
    current_hash = calculate_file_hash(prmtop_path)
    if current_hash != expected_hash:
        raise BatchError(
            f"CRITICAL: Receptor topology hash mismatch!\n"
            f"Expected: {expected_hash}\n"
            f"Current:  {current_hash}\n"
            f"The receptor topology has been modified during batch processing.\n"
            f"This would invalidate relative ranking comparisons.\n"
            f"Aborting batch run."
        )


def discover_ligand_files(ligands_input: Path) -> List[Path]:
    """
    Discover MOL2 ligand files from directory or explicit list.
    
    Args:
        ligands_input: Directory containing MOL2 files, or path to single MOL2
        
    Returns:
        List of MOL2 file paths
        
    Raises:
        BatchError: If no valid MOL2 files found or non-MOL2 files detected
    """
    ligand_files = []
    
    if ligands_input.is_dir():
        # Scan directory for MOL2 files
        ligand_files = sorted(ligands_input.glob("*.mol2"))
        
        # Check for non-MOL2 files and warn
        other_files = [
            f for f in ligands_input.iterdir()
            if f.is_file() and f.suffix.lower() not in [".mol2"]
        ]
        if other_files:
            logger.warning(
                f"Found {len(other_files)} non-MOL2 files in directory - ignoring"
            )
    
    elif ligands_input.is_file():
        # Single file - must be MOL2
        if ligands_input.suffix.lower() != ".mol2":
            raise BatchError(
                f"Only MOL2 ligands are accepted for batch processing.\n"
                f"File: {ligands_input.name}\n"
                f"Format: {ligands_input.suffix}\n"
                f"Reason: Batch mode requires consistent chemistry (GAFF2 + AM1-BCC)"
            )
        ligand_files = [ligands_input]
    
    else:
        raise BatchError(
            f"Ligands input not found: {ligands_input}\n"
            f"Expected: Directory of MOL2 files or single MOL2 file"
        )
    
    if not ligand_files:
        raise BatchError(
            f"No MOL2 ligands found in: {ligands_input}\n"
            f"Batch processing requires at least one ligand"
        )
    
    return ligand_files


def process_single_ligand(
    ligand_file: Path,
    receptor_topology_dir: Path,
    output_dir: Path,
    receptor_hash: str,
    method: str = "gb",
) -> Dict[str, Any]:
    """
    Process a single ligand against the pre-prepared receptor.
    
    This is a PURE ORCHESTRATOR - it calls existing functions without
    modifying any chemistry or physics logic.
    
    Args:
        ligand_file: Path to ligand MOL2
        receptor_topology_dir: Directory with protein.prmtop/inpcrd
        output_dir: Base output directory for batch
        receptor_hash: Expected receptor topology hash
        method: Implicit solvent method ("gb" or "pb")
        
    Returns:
        Dictionary with ligand_name and energy components
        
    Raises:
        BatchError: If processing fails or receptor integrity violated
    """
    ligand_name = ligand_file.stem
    method_label = "MM/GBSA" if method == "gb" else "MM/PBSA"
    logger.info(f"Processing ligand: {ligand_name}")
    
    # Verify receptor integrity before processing
    receptor_prmtop = receptor_topology_dir / "protein.prmtop"
    validate_receptor_topology_integrity(receptor_prmtop, receptor_hash)
    
    # Create ligand-specific subdirectory
    ligand_output_dir = output_dir / ligand_name
    ligand_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Parameterize ligand (GAFF2 + AM1-BCC)
        ligand_params_dir = ligand_output_dir / "ligand"
        ligand_params_dir.mkdir(exist_ok=True)
        
        parameterize_ligand(
            ligand_mol2=ligand_file,
            output_dir=ligand_params_dir,
        )
        logger.success(f"  ✓ Ligand parameterized")
        
        # Step 2: Assemble complex (combine receptor + ligand)
        complex_dir = ligand_output_dir / "complex"
        complex_dir.mkdir(exist_ok=True)
        
        prepare_complex(
            protein_dir=receptor_topology_dir,
            ligand_dir=ligand_params_dir,
            output_dir=complex_dir,
        )
        logger.success(f"  ✓ Complex assembled")
        
        # Step 3: Run rescoring
        rescore_dir = ligand_output_dir / "rescore"
        rescore_dir.mkdir(exist_ok=True)
        
        # Copy receptor topology to expected location
        # (needed because run_rescore expects protein/ subdirectory in parent of complex/)
        # Structure: ligand_output_dir/protein/, ligand_output_dir/ligand/, ligand_output_dir/complex/
        parent_protein_dir = ligand_output_dir / "protein"
        parent_protein_dir.mkdir(exist_ok=True)
        
        import shutil
        shutil.copy2(receptor_topology_dir / "protein.prmtop", parent_protein_dir / "protein.prmtop")
        shutil.copy2(receptor_topology_dir / "protein.inpcrd", parent_protein_dir / "protein.inpcrd")
        
        energies = run_rescore(
            complex_dir=complex_dir,
            output_dir=rescore_dir,
            method=method,
        )
        logger.success(
            f"  ✓ {method_label} complete: ΔG = {energies.get('DELTA_G', 'N/A'):.2f} kcal/mol"
        )
        
        # Verify receptor integrity after processing
        validate_receptor_topology_integrity(receptor_prmtop, receptor_hash)
        
        delta_g_gb = energies.get("DELTA_G_GB") if method == "gb" else None
        delta_g_pb = energies.get("DELTA_G_GB") if method == "pb" else None

        return {
            "ligand_name": ligand_name,
            "delta_g_bind": energies.get("DELTA_G"),
            "delta_h": energies.get("DELTA_H"),
            "delta_g_gb": delta_g_gb,
            "delta_g_pb": delta_g_pb,
            "status": "success",
        }
    
    except Exception as e:
        logger.error(f"  ✗ Failed: {e}")
        return {
            "ligand_name": ligand_name,
            "delta_g_bind": None,
            "delta_h": None,
            "delta_g_gb": None,
            "delta_g_pb": None,
            "status": f"failed: {str(e)}",
        }


def run_batch_rescore(
    receptor: Path,
    ligands: Path,
    output_dir: Path,
    skip_pdb4amber: bool = False,
    method: str = "gb",
) -> List[Dict[str, Any]]:
    """
    Run batch MM/GBSA rescoring for multiple ligands against one receptor.
    
    This is a PURE ORCHESTRATOR that:
    - Prepares receptor once
    - Reuses receptor topology for all ligands
    - Enforces chemistry consistency (MOL2 only, GAFF2 + AM1-BCC)
    - Validates receptor integrity throughout
    - Collects results for relative ranking
    
    NO MODIFICATIONS to core chemistry or physics code are made.
    
    Args:
        receptor: Path to receptor PDB
        ligands: Path to directory of MOL2 files or single MOL2
        output_dir: Output directory for batch results
        skip_pdb4amber: Skip pdb4amber processing
        method: Implicit solvent method ("gb" or "pb")
        
    Returns:
        List of result dictionaries for each ligand
        
    Raises:
        BatchError: If batch processing fails
    """
    method = method.lower()
    if method not in ["gb", "pb"]:
        raise BatchError("Invalid method. Expected 'gb' or 'pb'.")
    method_label = "MM/GBSA" if method == "gb" else "MM/PBSA"

    logger.info(f"Starting batch {method_label} rescoring (method={method.upper()})")
    logger.warning("⚠ RELATIVE RANKING ONLY - not thermodynamically rigorous ΔG")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Discover ligand files (MOL2 only)
    console.print("\n[bold cyan]Step 1:[/bold cyan] Discovering ligands...", style="cyan")
    ligand_files = discover_ligand_files(ligands)
    console.print(f"  ✓ Found {len(ligand_files)} MOL2 ligands", style="green")
    for ligand in ligand_files:
        console.print(f"    - {ligand.name}", style="dim")
    
    # Step 2: Prepare receptor once (reuse for all ligands)
    console.print("\n[bold cyan]Step 2:[/bold cyan] Preparing receptor (once)...", style="cyan")
    receptor_dir = output_dir / "receptor"
    receptor_dir.mkdir(exist_ok=True)
    
    prepare_protein(
        pdb_path=receptor,
        output_dir=receptor_dir,
        skip_pdb4amber=skip_pdb4amber,
    )
    console.print("  ✓ Receptor prepared", style="green")
    console.print(f"    {receptor_dir}/protein.prmtop", style="dim")
    console.print(f"    {receptor_dir}/protein.inpcrd", style="dim")
    
    # Convert receptor to PDB once (needed for complex assembly)
    import subprocess
    receptor_prmtop = receptor_dir / "protein.prmtop"
    receptor_inpcrd = receptor_dir / "protein.inpcrd"
    receptor_pdb = receptor_dir / "protein.pdb"
    
    with open(receptor_pdb, "w") as pdb_out:
        result = subprocess.run(
            ["ambpdb", "-p", "protein.prmtop", "-c", "protein.inpcrd"],
            stdout=pdb_out,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(receptor_dir),  # Run from receptor directory
        )
    
    if result.returncode != 0 or not receptor_pdb.exists():
        raise BatchError(
            f"Failed to convert receptor to PDB format:\n{result.stderr}"
        )
    
    console.print(f"    {receptor_dir}/protein.pdb", style="dim")
    
    # Calculate receptor topology hash (critical guardrail)
    receptor_hash = calculate_file_hash(receptor_prmtop)
    logger.info(f"Receptor topology hash: {receptor_hash[:16]}...")
    
    # Step 3: Process each ligand
    console.print("\n[bold cyan]Step 3:[/bold cyan] Processing ligands...", style="cyan")
    results = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"Processing {len(ligand_files)} ligands...",
            total=len(ligand_files)
        )
        
        for ligand_file in ligand_files:
            console.print(f"\n  [cyan]→ {ligand_file.name}[/cyan]")
            
            result = process_single_ligand(
                ligand_file=ligand_file,
                receptor_topology_dir=receptor_dir,
                output_dir=output_dir / "ligands",
                receptor_hash=receptor_hash,
                method=method,
            )
            results.append(result)
            
            progress.update(task, advance=1)
    
    # Step 4: Write results CSV
    console.print("\n[bold cyan]Step 4:[/bold cyan] Writing results...", style="cyan")
    csv_path = output_dir / "rescore_batch_results.csv"

    delta_field = "delta_g_gb" if method == "gb" else "delta_g_pb"
    with open(csv_path, "w", newline="") as csvfile:
        fieldnames = ["ligand_name", "delta_g_bind", "delta_h", delta_field, "status"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in results:
            writer.writerow({
                "ligand_name": result["ligand_name"],
                "delta_g_bind": result["delta_g_bind"],
                "delta_h": result["delta_h"],
                delta_field: result.get(delta_field),
                "status": result["status"],
            })
    
    console.print(f"  ✓ Results written to: {csv_path}", style="green")
    
    # Display summary table
    console.print("\n" + "=" * 60)
    console.print(f"[bold]BATCH {method_label} RESULTS[/bold]")
    console.print("=" * 60 + "\n")
    
    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Ligand", style="cyan")
    table.add_column("ΔG_bind", justify="right")
    table.add_column("ΔH", justify="right")
    table.add_column("ΔG_GB" if method == "gb" else "ΔG_PB", justify="right")
    table.add_column("Status", style="dim")
    
    # Sort by ΔG_bind (most favorable first)
    successful_results = [r for r in results if r["status"] == "success"]
    failed_results = [r for r in results if r["status"] != "success"]
    sorted_results = sorted(
        successful_results,
        key=lambda x: x["delta_g_bind"] if x["delta_g_bind"] is not None else float('inf')
    )
    
    for result in sorted_results + failed_results:
        dg = f"{result['delta_g_bind']:.2f}" if result['delta_g_bind'] is not None else "N/A"
        dh = f"{result['delta_h']:.2f}" if result['delta_h'] is not None else "N/A"
        dgpol = result.get(delta_field)
        dgb = f"{dgpol:.2f}" if dgpol is not None else "N/A"
        table.add_row(
            result["ligand_name"],
            dg,
            dh,
            dgb,
            result["status"]
        )
    
    console.print(table)
    
    # Summary statistics
    n_success = len(successful_results)
    n_failed = len(failed_results)
    console.print(f"\n[bold]Summary:[/bold]")
    console.print(f"  Total ligands:  {len(results)}")
    console.print(f"  Successful:     {n_success}", style="green")
    if n_failed > 0:
        console.print(f"  Failed:         {n_failed}", style="red")
    
    console.print("\n[yellow]⚠ CRITICAL REMINDER:[/yellow]")
    console.print("  These scores are for RELATIVE RANKING only")
    console.print("  NOT thermodynamically rigorous binding free energies")
    console.print("  Use for compound prioritization, not absolute ΔG values\n")
    
    logger.info(f"Batch processing complete: {n_success}/{len(results)} successful")
    
    return results
