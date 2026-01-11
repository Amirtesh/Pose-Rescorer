"""Command-line interface for physics-based post-docking rescoring.

This module provides a Typer-based CLI for running single-frame MM/GBSA / PBSA
rescoring (AmberTools MMPBSA.py) on protein-ligand complexes.
"""

from pathlib import Path
from typing import Optional

import typer
from Bio.PDB import PDBParser
from loguru import logger
from rich.console import Console
from rich.panel import Panel

from rescore import __version__
from rescore.validation import validate_pdb_complex, MMGBSAValidationError
from rescore.parameterization import (
    parameterize_ligand,
    ParameterizationError,
    ChemistryValidationError,
)
from rescore.parameterization.protein import (
    prepare_protein,
    ProteinPreparationError,
)
from rescore.parameterization.complex import (
    prepare_complex,
    ComplexAssemblyError,
)
from rescore.calculation import (
    run_rescore,
    MMGBSACalculationError,
)
from rescore.batch import (
    run_batch_rescore,
    BatchError,
)

app = typer.Typer(
    name="rescore",
    help="Physics-based post-docking rescoring tool (MM/GBSA under the hood)",
    add_completion=False,
)
console = Console()


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        console.print(f"rescore version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(                                                                                                         
        None,
        "--version",
        "-v",
        help="Show version and exit",                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        callback=version_callback,
        is_eager=True,
    ),
) -> None:
    """
    Physics-based post-docking rescoring tool.

    Uses MM/GBSA calculations (AmberTools) to rescore protein-ligand complexes:
    - Protein force field: ff14SB
    - Ligand force field: GAFF2
    - Charges: AM1-BCC
    - Solvation: GB/PB implicit solvent

    Designed for:
    - Post-docking relative ranking
    - Single-frame rescoring (fast)
    - Compound prioritization

    Limitations:
    - Relative ranking only (not absolute binding free energies)
    - No entropy calculations
    - No conformational sampling
    - No support for metals, covalent ligands, or multiple ligands
    """
    pass


@app.command()
def validate(
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="PDB file containing protein-ligand complex",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Validate a PDB structure for rescoring compatibility.

    Automatically detects ligands and checks for:
    - Exactly one ligand molecule (no zero, no multiple)
    - Missing backbone atoms (N, CA, C, O) in standard residues
    - Presence of metal atoms
    - Alternate locations (ALTLOC) or partial occupancy

    Note: Side-chain atoms are NOT validated in this version.

    Ligand detection:
    - Standard amino acids → protein
    - Water/ions (HOH, NA, CL, etc.) → solvent
    - Everything else → ligand

    Example:
        rescore validate complex.pdb
    """
    # Configure logging
    if verbose:
        logger.remove()
        logger.add(lambda msg: console.print(msg, end=""), level="DEBUG")
    else:
        logger.remove()
        logger.add(lambda msg: console.print(msg, end=""), level="INFO")

    console.print(
        Panel.fit(
            f"[bold]Structure Validation[/bold]\n"
            f"PDB: {pdb_file.name}",
            border_style="blue",
        )
    )

    try:
        ligand_resname = validate_pdb_complex(pdb_file)
        console.print("\n[bold green]✓ Validation passed[/bold green]")
        console.print(f"\nDetected ligand: [cyan]{ligand_resname}[/cyan]")
        console.print("Structure is compatible with rescoring workflow.")

    except MMGBSAValidationError as e:
        console.print(f"\n[bold red]✗ Validation failed[/bold red]")
        console.print(f"\n[red]{e}[/red]")
        raise typer.Exit(code=1)

    except Exception as e:
        console.print(f"\n[bold red]✗ Unexpected error[/bold red]")
        console.print(f"\n[red]{e}[/red]")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1)


@app.command()
@app.command(name="pdb-rescore", hidden=True)
def pdb_rescore(
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="PDB file containing protein-ligand complex",
    ),
    output_dir: Path = typer.Option(
        Path("rescore_output"),
        "--output",
        "-o",
        help="Output directory for results",
    ),
) -> None:
    """
    Run MM/GBSA rescoring on a protein-ligand complex.

    [NOT YET IMPLEMENTED]

    This command will:
    1. Automatically detect and validate the ligand
    2. Prepare protein with tleap (ff14SB)
    3. Parameterize ligand with antechamber (GAFF2, AM1-BCC)
    4. Run MM/GBSA calculation with MMPBSA.py
    5. Parse and report results

    Example:
        rescore pdb-rescore complex.pdb --output results/
    """
    console.print("[bold yellow]NOT YET IMPLEMENTED[/bold yellow]")
    console.print("\nUse 'rescore validate' to check structure compatibility.")
    raise typer.Exit(code=1)


@app.command()
def parameterize(
    ligand_mol2: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Ligand MOL2 file with complete chemistry (bond orders, hydrogens)",
    ),
    output_dir: Path = typer.Option(
        Path("ligand_params"),
        "--output-dir",
        "-o",
        help="Output directory for parameterization files",
    ),
    docked_pdb: Optional[Path] = typer.Option(
        None,
        "--pdb",
        "-p",
        exists=True,
        help="Docked complex PDB for coordinate transfer (optional)",
    ),
    ligand_resname: Optional[str] = typer.Option(
        None,
        "--resname",
        "-r",
        help="Ligand residue name in docked PDB (required with --pdb)",
    ),
    net_charge: int = typer.Option(
        0,
        "--charge",
        "-c",
        help="Net charge of ligand (default: 0)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Parameterize ligand with GAFF2 force field and AM1-BCC charges.

    CRITICAL: Requires MOL2 or SDF input with complete bond information.
    PDB files alone cannot be used because they lack bond orders.

    Workflow:
    1. If --pdb provided: Transfer docked coordinates to MOL2
    2. Run antechamber (GAFF2 + AM1-BCC charges)
    3. Run parmchk2 (generate frcmod for missing parameters)

    Why MOL2/SDF is required:
    - Antechamber needs explicit bond orders for atom typing
    - PDB format does not encode single/double/aromatic bonds
    - Inferring bonds from geometry is unreliable

    Prepare MOL2 files using:
    - PyMOL: save ligand.mol2
    - Avogadro: Export as MOL2
    - OpenBabel: obabel ligand.sdf -O ligand.mol2

    Examples:
        # Basic parameterization (MOL2 already has correct coordinates)
        rescore parameterize ligand.mol2 -o params/

        # With docked coordinates from PDB
        rescore parameterize ligand.mol2 --pdb complex.pdb --resname UNL -o params/

        # Charged ligand
        rescore parameterize ligand.mol2 --charge -1 -o params/
    """
    # Validate arguments
    if docked_pdb and not ligand_resname:
        console.print("[bold red]Error:[/bold red] --resname is required when --pdb is provided")
        raise typer.Exit(code=1)
    
    if not verbose:
        logger.remove()
        logger.add(
            lambda msg: None,
            filter=lambda record: record["level"].no >= logger.level("SUCCESS").no
        )

    # Display header
    console.print()
    console.print(
        Panel.fit(
            f"[bold]Ligand Parameterization[/bold]\n"
            f"MOL2: {ligand_mol2.name}",
            border_style="blue",
        )
    )

    try:
        # Parameterize ligand
        console.print("Running ligand parameterization...", style="yellow")
        
        input_mol2, output_mol2, frcmod = parameterize_ligand(
            ligand_mol2=ligand_mol2,
            output_dir=output_dir,
            net_charge=net_charge,
            docked_pdb=docked_pdb,
            ligand_resname=ligand_resname,
        )

        # Success message
        console.print()
        console.print("✓ Parameterization complete", style="green bold")
        console.print()
        console.print("[bold]Output files:[/bold]")
        console.print(f"  - Input:  {input_mol2}")
        console.print(f"  - Output: {output_mol2}")
        console.print(f"  - Frcmod: {frcmod}")

    except ParameterizationError as e:
        console.print()
        console.print("[bold red]✗ Parameterization failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during parameterization")
        raise typer.Exit(code=1)


@app.command()
def prep_protein(
    protein_pdb: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Protein-only PDB file (no ligands, waters, or ions)",
    ),
    output_dir: Path = typer.Option(
        Path("protein_params"),
        "--output-dir",
        "-o",
        help="Output directory for protein topology files",
    ),
    skip_pdb4amber: bool = typer.Option(
        False,
        "--skip-pdb4amber",
        help="Skip pdb4amber preprocessing (not recommended)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Prepare protein topology using AmberTools tleap with ff14SB.

    DEFAULT WORKFLOW (recommended):
    1. Run pdb4amber to normalize PDB (fixes atom names, removes waters)
    2. Validate protein-only structure
    3. Run tleap with ff14SB to generate topology

    pdb4amber preprocessing:
    - Converts non-standard atom names (HN → H)
    - Removes waters and ions
    - Fixes residue naming
    - Handles alternate conformations

    Supported residues:
    - Standard amino acids (ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, etc.)
    - Histidine variants (HIS, HID, HIE, HIP, HSE, HSD, HSP)
    - Terminal variants (automatically handled by tleap)

    The command will FAIL if:
    - Non-protein residues present (ligands, waters, ions) after pdb4amber
    - Missing atoms in structure
    - Non-standard residue names
    - Incomplete residues

    Prepare protein-only PDB using:
    - PyMOL: select protein, polymer.protein; save protein.pdb, protein
    - Chimera: select protein; save protein.pdb

    Examples:
        # Basic protein preparation (with pdb4amber)
        rescore prep-protein receptor.pdb -o protein_params/

        # Skip pdb4amber (not recommended)
        rescore prep-protein receptor.pdb -o protein_params/ --skip-pdb4amber

        # With verbose output (see pdb4amber and tleap logs)
        rescore prep-protein receptor.pdb -o protein_params/ -v
    """
    if not verbose:
        logger.remove()
        logger.add(
            lambda msg: None,
            filter=lambda record: record["level"].no >= logger.level("SUCCESS").no
        )

    # Display header
    console.print()
    console.print(
        Panel.fit(
            f"[bold]Protein Preparation[/bold]\n"
            f"PDB: {protein_pdb.name}\n"
            f"Force Field: ff14SB\n"
            f"pdb4amber: {'DISABLED' if skip_pdb4amber else 'ENABLED'}",
            border_style="blue",
        )
    )

    try:
        console.print("Running protein preparation...", style="yellow")
        
        prmtop, inpcrd = prepare_protein(
            pdb_path=protein_pdb,
            output_dir=output_dir,
            skip_pdb4amber=skip_pdb4amber,
        )

        # Success message
        console.print()
        console.print("✓ Protein preparation complete", style="green bold")
        console.print()
        console.print("[bold]Output files:[/bold]")
        console.print(f"  - Topology:    {prmtop}")
        console.print(f"  - Coordinates: {inpcrd}")
        console.print()
        console.print("[dim]Next steps:[/dim]")
        console.print("[dim]  1. Parameterize ligand: rescore parameterize ligand.mol2[/dim]")
        console.print("[dim]  2. Assemble complex: rescore assemble ... [/dim]")
        console.print("[dim]  3. Run rescoring: rescore run ... (not yet implemented)[/dim]")

    except ProteinPreparationError as e:
        console.print()
        console.print("[bold red]✗ Protein preparation failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during protein preparation")
        raise typer.Exit(code=1)


@app.command()
def assemble(
    protein_dir: Path = typer.Option(
        ...,
        "--protein",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Directory containing protein.prmtop and protein.inpcrd",
    ),
    ligand_dir: Path = typer.Option(
        ...,
        "--ligand",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Directory containing ligand.mol2 and ligand.frcmod",
    ),
    output_dir: Path = typer.Option(
        Path("complex_params"),
        "--output",
        "-o",
        help="Output directory for complex topology",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Assemble protein-ligand complex for rescoring calculation.

    This command combines:
    - Protein topology (from prep-protein)
    - Ligand topology (from parameterize)

    into a single-frame complex suitable for MM/GBSA rescoring.

    WORKFLOW:
    1. Validate all input files exist
    2. Combine protein + ligand using tleap

    CRITICAL NOTES:
    - This does NOT perform solvation, minimization, or MD
    - Protein coordinates are taken from protein.inpcrd
    - Ligand coordinates are taken from ligand.mol2
    - Chemistry (atom types, bonds, charges) is preserved from parameters

    REQUIREMENTS:
    - Run prep-protein first to generate protein parameters
    - Run parameterize first to generate ligand parameters

    Examples:
        # Basic complex assembly
        rescore assemble \
          --protein receptor_params/ \
          --ligand ligand_params/ \
          -o complex_params/

        # With verbose output
        rescore assemble \
          --protein receptor_params/ \
          --ligand ligand_params/ \
          -o complex_params/ -v
    """
    if not verbose:
        logger.remove()
        logger.add(
            lambda msg: None,
            filter=lambda record: record["level"].no >= logger.level("SUCCESS").no
        )

    # Display header
    console.print()
    console.print(
        Panel.fit(
            f"[bold]Complex Assembly[/bold]\n"
            f"Protein: {protein_dir.name}/\n"
            f"Ligand: {ligand_dir.name}/",
            border_style="blue",
        )
    )

    try:
        console.print("Assembling complex...", style="yellow")
        
        prepare_complex(
            protein_dir=protein_dir,
            ligand_dir=ligand_dir,
            output_dir=output_dir,
        )

        # Success message
        console.print()
        console.print("✓ Complex assembly complete", style="green bold")
        console.print()
        console.print("[bold]Output files:[/bold]")
        console.print(f"  - Topology:    {output_dir}/complex.prmtop")
        console.print(f"  - Coordinates: {output_dir}/complex.inpcrd")
        console.print()
        console.print("[dim]Next steps:[/dim]")
        console.print(f"[dim]  1. Run rescoring: rescore run {output_dir}/ -o results/[/dim]")

    except ComplexAssemblyError as e:
        console.print()
        console.print("[bold red]✗ Complex assembly failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during complex assembly")
        raise typer.Exit(code=1)


@app.command()
def run(
    complex_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Directory containing complex.prmtop and complex.inpcrd",
    ),
    output_dir: Path = typer.Option(
        Path("rescore_results"),
        "--output",
        "-o",
        help="Output directory for MM/GBSA results",
    ),
    method: str = typer.Option(
        "gb",
        "--method",
        "-m",
        help="Solvation method: 'gb' (default, fast) or 'pb' (slower, not more rigorous)",
        case_sensitive=False,
    ),
    minimize: bool = typer.Option(
        True,
        "--minimize/--no-minimize",
        help="Run restrained minimization before scoring (default: ON)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Run single-frame MM/GBSA rescoring.

    This command performs MM/GBSA rescoring for post-docking analysis:
    - Single structure (no MD trajectory)
    - Optional restrained minimization (removes clashes)
    - GB or PB implicit solvent
    - No entropy calculation
    - No decomposition

    CRITICAL: This is POST-DOCKING RESCORING, not thermodynamic ΔG
    
    Solvation Methods:
    - GB (default): Fast, recommended for most cases
    - PB: Slower, more sensitive to structure quality, NOT more rigorous in single-frame mode
    
    Use for:
    - Relative ranking of docked poses
    - Compound prioritization
    - Initial screening
    
    Do NOT use for:
    - Absolute binding free energies
    - Publication-quality ΔG values
    - Entropy-dominated systems
    
    Method: MM/GBSA with Generalized Born implicit solvent
    Output: ΔG_bind and energy components (kcal/mol)
    
    Requirements:
    - Run assemble first to generate complex topology
    - AmberTools MMPBSA.py must be available

    Examples:
        # Basic rescoring
        rescore run complex_params/ -o results/

        # With verbose output
        rescore run complex_params/ -o results/ --verbose
    """
    if not verbose:
        logger.remove()
        logger.add(
            lambda msg: None,
            filter=lambda record: record["level"].no >= logger.level("SUCCESS").no
        )

    # Display header
    # Validate method parameter
    method = method.lower()
    if method not in ["gb", "pb"]:
        console.print()
        console.print(f"[bold red]✗ Invalid method:[/bold red] {method}")
        console.print()
        console.print("Method must be either 'gb' or 'pb'")
        console.print("  - gb: Generalized Born (default, fast)")
        console.print("  - pb: Poisson-Boltzmann (slower, not more rigorous)")
        raise typer.Exit(code=1)
    method_label = "MM/GBSA" if method == "gb" else "MM/PBSA"
    method_label = "MM/GBSA" if method == "gb" else "MM/PBSA"
    
    console.print()
    console.print(
        Panel.fit(
            f"[bold]{method_label} Rescoring[/bold]\n"
            f"Complex: {complex_dir.name}/\n"
            f"Method: Single-frame {method.upper()}\n"
            f"[yellow]⚠ POST-DOCKING RESCORING[/yellow]\n"
            f"[dim]Not thermodynamic ΔG[/dim]",
            border_style="blue",
        )
    )

    try:
        console.print(f"Running {method_label} calculation ({method.upper()})...", style="yellow")
        
        energies = run_rescore(
            complex_dir=complex_dir,
            output_dir=output_dir,
            method=method,
            minimize=minimize,
        )

        # Success message with results
        console.print()
        console.print(f"✓ {method_label} calculation complete", style="green bold")
        console.print()
        console.print("[bold]Energy Components (kcal/mol):[/bold]")
        
        if "DELTA_G" in energies:
            console.print(f"  ΔG_bind:  {energies['DELTA_G']:>8.2f}")
        if "DELTA_H" in energies:
            console.print(f"  ΔH:       {energies['DELTA_H']:>8.2f}")
        delta_g_label = "ΔG_GB" if method == "gb" else "ΔG_PB"
        if "DELTA_G_GB" in energies:
            console.print(f"  {delta_g_label}:    {energies['DELTA_G_GB']:>8.2f}")
        if "DELTA_G_SA" in energies:
            console.print(f"  ΔG_SA:    {energies['DELTA_G_SA']:>8.2f}")
        
        console.print()
        console.print("[bold]Output files:[/bold]")
        console.print(f"  - Results: {output_dir}/rescore_output.dat")
        console.print(f"  - Input:   {output_dir}/rescore.in")
        console.print(f"  - Log:     {output_dir}/rescore.log")
        console.print()
        console.print("[yellow]⚠ IMPORTANT:[/yellow]")
        console.print("[dim]  This is RESCORING for relative ranking only[/dim]")
        console.print("[dim]  Not suitable for absolute ΔG values[/dim]")

    except MMGBSACalculationError as e:
        console.print()
        console.print("[bold red]✗ MM/GBSA calculation failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during MM/GBSA calculation")
        raise typer.Exit(code=1)


@app.command()
def integrate(
    receptor: Path = typer.Option(
        ...,
        "--receptor",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Receptor PDB file (protein structure)",
    ),
    ligand: Path = typer.Option(
        ...,
        "--ligand",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Ligand MOL2 or SDF file (must include bond connectivity)",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output directory for complete workflow",
    ),
    method: str = typer.Option(
        "gb",
        "--method",
        "-m",
        help="Solvation method: 'gb' (default, fast) or 'pb' (slower, not more rigorous)",
        case_sensitive=False,
    ),
    skip_validation: bool = typer.Option(
        False,
        "--skip-validation",
        help="Skip structure validation step",
    ),
    skip_pdb4amber: bool = typer.Option(
        False,
        "--skip-pdb4amber",
        help="Skip pdb4amber preprocessing (not recommended)",
    ),
    minimize: bool = typer.Option(
        True,
        "--minimize/--no-minimize",
        help="Run restrained minimization before scoring (default: ON)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging",
    ),
) -> None:
    """
    Integrated pipeline: receptor + ligand → rescoring in one command.

    This is the PRIMARY USER-FACING COMMAND for the package. It executes
    the complete MM/GBSA workflow from raw inputs to results:

    WORKFLOW (executed sequentially):
    1. Structure validation (optional)
    2. Protein preparation (pdb4amber + tleap)
    3. Ligand parameterization (GAFF2 + AM1-BCC)
    4. Complex assembly (combine topologies)
    5. Restrained minimization (optional, removes clashes)
    6. MM/GBSA single-frame rescoring (GB or PB)

    OUTPUT DIRECTORY STRUCTURE:
    output/
    ├── protein/      (protein.prmtop, protein.inpcrd)
    ├── ligand/       (ligand.mol2, ligand.frcmod)
    ├── complex/      (complex.prmtop, complex.inpcrd)
    └── rescore/       (rescore_output.dat, rescore.log, minimization/)

    REQUIREMENTS:
    - Receptor: PDB file (may contain non-standard naming)
    - Ligand: MOL2 or SDF with full bond connectivity (NOT PDB)

    CRITICAL NOTES:
    - This is POST-DOCKING RESCORING for relative ranking
    - Not thermodynamically rigorous binding free energy
    - Optional minimization relaxes clashes (NOT MD)
    - Fail-fast on any error (no silent fixes)

    Examples:
        # Basic workflow
        rescore integrate \
          --receptor receptor.pdb \
          --ligand ligand.mol2 \
          -o results/

        # Skip validation (faster, but risky)
        rescore integrate \
          --receptor receptor.pdb \
          --ligand ligand.mol2 \
          -o results/ \
          --skip-validation

        # With verbose output
        rescore integrate \
          --receptor receptor.pdb \
          --ligand ligand.mol2 \
          -o results/ -v
    """
    if not verbose:
        logger.remove()
        logger.add(
            lambda msg: None,
            filter=lambda record: record["level"].no >= logger.level("SUCCESS").no
        )
    
    # Validate method parameter
    method = method.lower()
    if method not in ["gb", "pb"]:
        console.print()
        console.print(f"[bold red]✗ Invalid method:[/bold red] {method}")
        console.print()
        console.print("Method must be either 'gb' or 'pb'")
        console.print("  - gb: Generalized Born (default, fast)")
        console.print("  - pb: Poisson-Boltzmann (slower, not more rigorous)")
        raise typer.Exit(code=1)
    method_label = "MM/GBSA" if method == "gb" else "MM/PBSA"

    # Display header
    console.print()
    console.print(
        Panel.fit(
            f"[bold]{method_label} Integrated Pipeline[/bold]\n"
            f"Receptor: {receptor.name}\n"
            f"Ligand: {ligand.name}\n"
            f"Output: {output_dir}/\n"
            f"Method: {method.upper()}\n"
            f"[yellow]⚠ POST-DOCKING RESCORING[/yellow]",
            border_style="blue",
        )
    )
    console.print()

    # Setup directory structure
    protein_dir = output_dir / "protein"
    ligand_dir = output_dir / "ligand"
    complex_dir = output_dir / "complex"
    rescore_dir = output_dir / "rescore"

    try:
        # ============================================================
        # STEP 1: Protein Preparation
        # ============================================================
        console.print("[bold cyan]Step 1/5:[/bold cyan] Preparing protein...", style="cyan")
        
        prepare_protein(
            pdb_path=receptor,
            output_dir=protein_dir,
            skip_pdb4amber=skip_pdb4amber,
        )
        console.print("  ✓ Protein prepared", style="green")
        console.print(f"    {protein_dir}/protein.prmtop", style="dim")
        console.print(f"    {protein_dir}/protein.inpcrd", style="dim")

        # ============================================================
        # STEP 2: Ligand Parameterization
        # ============================================================
        console.print()
        console.print("[bold cyan]Step 2/5:[/bold cyan] Parameterizing ligand...", style="cyan")
        
        # Validate ligand format
        ligand_suffix = ligand.suffix.lower()
        if ligand_suffix not in [".mol2", ".sdf"]:
            console.print()
            console.print("[bold red]✗ Invalid ligand format[/bold red]")
            console.print()
            console.print(f"Ligand must be MOL2 or SDF format, not {ligand_suffix}")
            console.print("PDB ligands are NOT supported (no bond information)")
            raise typer.Exit(code=1)
        
        parameterize_ligand(
            ligand_mol2=ligand,
            output_dir=ligand_dir,
        )
        console.print("  ✓ Ligand parameterized", style="green")
        console.print(f"    {ligand_dir}/ligand.mol2", style="dim")
        console.print(f"    {ligand_dir}/ligand.frcmod", style="dim")

        # ============================================================
        # STEP 3: Complex Assembly
        # ============================================================
        console.print()
        console.print("[bold cyan]Step 3/5:[/bold cyan] Assembling complex...", style="cyan")
        
        prepare_complex(
            protein_dir=protein_dir,
            ligand_dir=ligand_dir,
            output_dir=complex_dir,
        )
        console.print("  ✓ Complex assembled", style="green")
        console.print(f"    {complex_dir}/complex.prmtop", style="dim")
        console.print(f"    {complex_dir}/complex.inpcrd", style="dim")

        # ============================================================
        # STEP 4: Structure Validation (Optional)
        # ============================================================
        if not skip_validation:
            console.print()
            console.print("[bold cyan]Step 4/5:[/bold cyan] Validating assembled complex...", style="cyan")
            
            # Convert complex to PDB for validation
            import subprocess
            complex_pdb = complex_dir / "complex.pdb"
            result = subprocess.run(
                ["ambpdb", "-p", str(complex_dir / "complex.prmtop"), 
                 "-c", str(complex_dir / "complex.inpcrd")],
                stdout=open(complex_pdb, "w"),
                stderr=subprocess.PIPE,
                text=True,
            )
            
            if result.returncode == 0 and complex_pdb.exists():
                try:
                    validate_pdb_complex(complex_pdb)
                    console.print("  ✓ Validation passed", style="green")
                except MMGBSAValidationError as e:
                    console.print()
                    console.print("[bold red]✗ Validation failed[/bold red]")
                    console.print()
                    console.print(str(e))
                    console.print()
                    console.print("[dim]Hint: Use --skip-validation to bypass (not recommended)[/dim]")
                    raise typer.Exit(code=1)
            else:
                console.print("  [yellow]⚠ Could not generate PDB for validation[/yellow]", style="yellow")
        else:
            console.print()
            console.print("[bold cyan]Step 4/5:[/bold cyan] Validation skipped", style="yellow")

        # ============================================================
        # STEP 5: MM/GBSA Calculation
        # ============================================================
        console.print()
        console.print(f"[bold cyan]Step 5/5:[/bold cyan] Running {method_label} rescoring...", style="cyan")
        
        energies = run_rescore(
            complex_dir=complex_dir,
            output_dir=rescore_dir,
            method=method,
            minimize=minimize,
        )
        console.print(f"  ✓ {method_label} complete", style="green")

        # ============================================================
        # FINAL RESULTS
        # ============================================================
        console.print()
        console.print("=" * 60)
        console.print()
        console.print(f"[bold green]✓ {method_label} PIPELINE COMPLETE[/bold green]")
        console.print()
        console.print("[bold]Energy Components (kcal/mol):[/bold]")
        
        if "DELTA_G" in energies:
            console.print(f"  ΔG_bind:  {energies['DELTA_G']:>8.2f}")
        if "DELTA_H" in energies:
            console.print(f"  ΔH:       {energies['DELTA_H']:>8.2f}")
        delta_g_label = "ΔG_GB" if method == "gb" else "ΔG_PB"
        if "DELTA_G_GB" in energies:
            console.print(f"  {delta_g_label}:    {energies['DELTA_G_GB']:>8.2f}")
        
        console.print()
        console.print("[bold]Output Structure:[/bold]")
        console.print(f"  {output_dir}/")
        console.print(f"  ├── protein/    (receptor parameters)")
        console.print(f"  ├── ligand/     (ligand parameters)")
        console.print(f"  ├── complex/    (combined topology)")
        console.print(f"  └── rescore/     (rescoring results)")
        console.print()
        console.print("[yellow]⚠ CRITICAL REMINDER:[/yellow]")
        console.print("[dim]  This is POST-DOCKING RESCORING for relative ranking[/dim]")
        console.print("[dim]  NOT thermodynamically rigorous binding free energy[/dim]")
        console.print("[dim]  Use for compound prioritization, not absolute ΔG values[/dim]")
        console.print()

    except ProteinPreparationError as e:
        console.print()
        console.print("[bold red]✗ Protein preparation failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except (ParameterizationError, ChemistryValidationError) as e:
        console.print()
        console.print("[bold red]✗ Ligand parameterization failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except ComplexAssemblyError as e:
        console.print()
        console.print("[bold red]✗ Complex assembly failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except MMGBSACalculationError as e:
        console.print()
        console.print(f"[bold red]✗ {method_label} calculation failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)

    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during integrated pipeline")
        raise typer.Exit(code=1)


@app.command()
def batch(
    receptor: Path = typer.Option(
        ...,
        "--receptor",
        "-r",
        help="Receptor PDB file",
        exists=True,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    ),
    ligands: Path = typer.Option(
        ...,
        "--ligands",
        "-l",
        help="Directory containing MOL2 files or single MOL2 file",
        exists=True,
        resolve_path=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output directory for batch results",
        resolve_path=True,
    ),
    skip_pdb4amber: bool = typer.Option(
        False,
        "--skip-pdb4amber",
        help="Skip pdb4amber processing (use if receptor already prepared)",
    ),
    method: str = typer.Option(
        "gb",
        "--method",
        "-m",
        help="Implicit solvent method: gb (Generalized Born) or pb (Poisson-Boltzmann)",
        case_sensitive=False,
    ),
    minimize: bool = typer.Option(
        True,
        "--minimize/--no-minimize",
        help="Run restrained minimization before scoring (default: ON)",
    ),
) -> None:
    """
    Batch rescoring for multiple ligands against one receptor.
    
    This command is a PURE ORCHESTRATOR for high-throughput rescoring:
    
    \b
    - Prepares receptor once (topology reused for all ligands)
    - Processes each ligand: parameterize → assemble → minimize (optional) → rescore
    - Validates receptor integrity throughout (hash checking)
    - Outputs CSV with all results for relative ranking
    
    \b
    Solvation Methods:
    - GB (default): Faster, sufficient for relative ranking
    - PB: Slower, more sensitive to structure quality
      WARNING: In single-frame mode, PB does NOT make results more rigorous!
    
    \b
    CRITICAL GUARDRAILS:
    - Only MOL2 ligands accepted (consistent chemistry)
    - Receptor topology must not change during batch
    - Results are for RELATIVE RANKING ONLY
    
    \b
    Example:
        rescore batch -r protein.pdb -l ligands/ -o batch_results/
        rescore batch -r protein.pdb -l ligands/ -o batch_results/ --method pb
    
    \b
    Output structure:
        batch_results/
        ├── receptor/              (prepared once)
        ├── ligands/               (per-ligand subdirs)
        │   ├── ligand1/
        │   │   ├── ligand/
        │   │   ├── complex/
        │   │   └── rescore/
        │   └── ligand2/
        │       └── ...
        └── rescore_batch_results.csv
    """
    logger.remove()  # Remove default handler
    logger.add(
        lambda msg: None,  # Suppress loguru output (use Rich for user messages)
        level="INFO"
    )
    
    # Validate method
    method = method.lower()
    if method not in ["gb", "pb"]:
        console.print("[bold red]✗ Error:[/bold red] Method must be 'gb' or 'pb'")
        raise typer.Exit(code=1)
    
    console.print()
    console.print(
        Panel.fit(
            f"[bold]MM/GBSA Batch Rescoring[/bold]\n"
            f"Receptor: {receptor.name}\n"
            f"Ligands: {ligands}\n"
            f"Method: {method.upper()}\n"
            f"Output: {output_dir}/\n"
            f"[yellow]⚠ RELATIVE RANKING ONLY[/yellow]",
            border_style="blue",
        )
    )
    console.print()
    
    try:
        run_batch_rescore(
            receptor=receptor,
            ligands=ligands,
            output_dir=output_dir,
            skip_pdb4amber=skip_pdb4amber,
            method=method,
            minimize=minimize,
        )
        
        console.print("[bold green]✓ Batch processing complete![/bold green]\n")
    
    except BatchError as e:
        console.print()
        console.print("[bold red]✗ Batch processing failed[/bold red]")
        console.print()
        console.print(str(e))
        console.print()
        raise typer.Exit(code=1)
    
    except Exception as e:
        console.print()
        console.print(f"[bold red]✗ Unexpected error:[/bold red] {e}")
        logger.exception("Unexpected error during batch processing")
        raise typer.Exit(code=1)
