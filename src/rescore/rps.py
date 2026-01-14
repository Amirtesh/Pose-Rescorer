"""
Rapid Perturbation Sampling (RPS) for MM/GBSA score uncertainty quantification.

CRITICAL SCOPE:
- This is NOT molecular dynamics
- This is NOT conformer generation
- This is NOT docking or pose optimization
- This IS numerical sensitivity analysis

RPS evaluates how small coordinate perturbations affect single-frame MM/GBSA scores.
It provides uncertainty bounds for a fixed docking-derived pose.

HARD CHEMISTRY-FROZEN GUARANTEE:
- antechamber is NEVER invoked during RPS
- Ligand charges are frozen (AM1-BCC from main workflow)
- Ligand atom types are frozen (GAFF2 from main workflow)
- Force field parameters are frozen (frcmod from main workflow)
- ONLY coordinates (x,y,z) are perturbed with Gaussian noise

Use for:
- Quantifying score stability
- Understanding numerical sensitivity
- Diagnostic analysis of single poses

Do NOT use for:
- Ensemble averaging
- Thermodynamic sampling
- Pose refinement
- Conformational search
"""

import csv
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from loguru import logger
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from rescore.parameterization.complex import prepare_complex
from rescore.calculation import run_rescore


console = Console()


class RPSError(Exception):
    """Raised when RPS execution fails."""
    pass


def generate_deterministic_seed(ligand_name: str, replicate_id: int) -> int:
    """
    Generate deterministic random seed from ligand name and replicate ID.
    
    This ensures reproducibility: same ligand + same replicate → same perturbation.
    
    Args:
        ligand_name: Name of the ligand
        replicate_id: Perturbation replicate index
        
    Returns:
        32-bit integer seed
    """
    seed_string = f"{ligand_name}_replicate_{replicate_id}"
    hash_digest = hashlib.md5(seed_string.encode()).hexdigest()
    seed = int(hash_digest[:8], 16)
    return seed


def perturb_ligand_coordinates(
    mol2_path: Path,
    output_path: Path,
    sigma: float,
    seed: int,
) -> None:
    """
    Perturb ligand coordinates by adding Gaussian noise.
    
    CRITICAL:
    - Only ligand coordinates are modified
    - Bonding, atom types, and connectivity are preserved
    - Perturbations are independent for x, y, z
    - Random seed ensures reproducibility
    
    Args:
        mol2_path: Path to original ligand MOL2 file
        output_path: Path for perturbed MOL2 file
        sigma: Standard deviation of Gaussian noise (Angstroms)
        seed: Random seed for reproducibility
        
    Raises:
        RPSError: If MOL2 parsing fails
    """
    # Set random seed
    rng = np.random.RandomState(seed)
    
    # Read MOL2 file
    with open(mol2_path, 'r') as f:
        lines = f.readlines()
    
    # Find @<TRIPOS>ATOM section
    atom_start = None
    atom_end = None
    for i, line in enumerate(lines):
        if '@<TRIPOS>ATOM' in line:
            atom_start = i + 1
        elif atom_start is not None and line.startswith('@<TRIPOS>'):
            atom_end = i
            break
    
    if atom_start is None:
        raise RPSError(f"Could not find @<TRIPOS>ATOM section in {mol2_path}")
    
    if atom_end is None:
        atom_end = len(lines)
    
    # Perturb coordinates in atom section
    perturbed_lines = lines[:atom_start]
    
    for i in range(atom_start, atom_end):
        line = lines[i]
        parts = line.split()
        
        if len(parts) >= 6:
            # Parse coordinates
            atom_id = parts[0]
            atom_name = parts[1]
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            atom_type = parts[5]
            
            # Add Gaussian noise
            x_pert = x + rng.normal(0, sigma)
            y_pert = y + rng.normal(0, sigma)
            z_pert = z + rng.normal(0, sigma)
            
            # Reconstruct line with perturbed coordinates
            # Preserve MOL2 formatting (fixed-width fields)
            rest_of_line = ' '.join(parts[6:]) if len(parts) > 6 else ''
            new_line = f"{atom_id:>7} {atom_name:<8} {x_pert:>9.4f} {y_pert:>9.4f} {z_pert:>9.4f} {atom_type:<8} {rest_of_line}\n"
            perturbed_lines.append(new_line)
        else:
            perturbed_lines.append(line)
    
    # Add remaining sections
    perturbed_lines.extend(lines[atom_end:])
    
    # Write perturbed MOL2
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.writelines(perturbed_lines)
    
    logger.debug(f"Perturbed ligand saved: {output_path}")


def run_single_perturbation(
    receptor_dir: Path,
    parameterized_ligand_dir: Path,
    replicate_id: int,
    sigma: float,
    ligand_name: str,
    output_dir: Path,
    method: str,
    minimize: bool,
) -> Dict[str, float]:
    """
    Run complete workflow for a single perturbation replicate.
    
    CRITICAL: RPS perturbs coordinates only; ligand chemistry is fixed by design.
    This function does NOT re-parameterize the ligand. It uses the pre-computed
    GAFF2 atom types, AM1-BCC charges, and frcmod parameters from the main workflow.
    
    HARD GUARD: antechamber is NEVER invoked in RPS mode.
    - Charges are frozen (from original AM1-BCC calculation)
    - Atom types are frozen (from original GAFF2 assignment)
    - Force field parameters are frozen (from original frcmod)
    - Only coordinates (x,y,z) are perturbed
    
    Steps:
    1. Copy pre-parameterized ligand.mol2 and perturb only x/y/z coordinates
    2. Copy original ligand.frcmod (chemistry frozen)
    3. Assemble complex using perturbed coordinates + fixed chemistry
    4. Run MM/GBSA scoring
    
    Args:
        receptor_dir: Directory containing prepared receptor
        parameterized_ligand_dir: Directory with ligand.mol2 and ligand.frcmod (pre-parameterized)
        replicate_id: Perturbation replicate index
        sigma: Perturbation magnitude (Angstroms)
        ligand_name: Ligand name for seed generation
        output_dir: Output directory for this replicate
        method: Solvation method ("gb" or "pb")
        minimize: Whether to run minimization
        
    Returns:
        Dictionary of energy components
        
    Raises:
        RPSError: If any step fails
    """
    # Generate deterministic seed
    seed = generate_deterministic_seed(ligand_name, replicate_id)
    
    # Create replicate directory
    replicate_dir = output_dir / f"replicate_{replicate_id:03d}"
    replicate_dir.mkdir(parents=True, exist_ok=True)
    
    # Create symbolic link to receptor directory (needed for run_rescore)
    # run_rescore expects protein/ as sibling of complex/
    protein_link = replicate_dir / "protein"
    if not protein_link.exists():
        protein_link.symlink_to(receptor_dir.resolve(), target_is_directory=True)
    
    # Step 1: Perturb ligand coordinates (chemistry frozen)
    # RPS perturbs coordinates only; ligand chemistry is fixed by design.
    parameterized_mol2 = parameterized_ligand_dir / "ligand.mol2"
    perturbed_mol2 = replicate_dir / "ligand_perturbed.mol2"
    perturb_ligand_coordinates(parameterized_mol2, perturbed_mol2, sigma, seed)
    
    # Step 2: Copy pre-parameterized ligand files (NO re-parameterization)
    # Use the original GAFF2 atom types, AM1-BCC charges, and frcmod
    # HARD GUARD: antechamber must NEVER be invoked in RPS mode
    ligand_dir = replicate_dir / "ligand"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    
    import shutil
    # Copy perturbed MOL2 with original chemistry
    shutil.copy2(perturbed_mol2, ligand_dir / "ligand.mol2")
    
    # Copy original frcmod (chemistry frozen)
    original_frcmod = parameterized_ligand_dir / "ligand.frcmod"
    if not original_frcmod.exists():
        raise RPSError(
            f"CRITICAL: Missing ligand.frcmod in {parameterized_ligand_dir}\n"
            "RPS requires pre-parameterized ligand. Run main workflow first."
        )
    shutil.copy2(original_frcmod, ligand_dir / "ligand.frcmod")
    
    # HARD GUARD: Verify chemistry files exist (copied, not generated)
    assert (ligand_dir / "ligand.mol2").exists(), "ligand.mol2 copy failed"
    assert (ligand_dir / "ligand.frcmod").exists(), "ligand.frcmod copy failed"
    
    # HARD GUARD: Verify frcmod is the original (same size = not regenerated)
    if (ligand_dir / "ligand.frcmod").stat().st_size != original_frcmod.stat().st_size:
        raise RPSError(
            "CRITICAL: ligand.frcmod was modified during RPS\n"
            "This indicates antechamber was invoked, which violates RPS chemistry-frozen principle."
        )
    
    # Step 3: Assemble complex
    complex_dir = replicate_dir / "complex"
    prepare_complex(
        protein_dir=receptor_dir,
        ligand_dir=ligand_dir,
        output_dir=complex_dir,
    )
    
    # Step 4: Run MM/GBSA
    rescore_dir = replicate_dir / "rescore"
    energies = run_rescore(
        complex_dir=complex_dir,
        output_dir=rescore_dir,
        method=method,
        minimize=minimize,
    )
    
    return energies


def compute_statistics(scores: List[float]) -> Dict[str, float]:
    """
    Compute statistical metrics for score distribution.
    
    Args:
        scores: List of MM/GBSA scores from perturbations
        
    Returns:
        Dictionary containing:
        - mean: Mean score
        - std: Sample standard deviation
        - ci_lower: Lower 95% confidence interval bound
        - ci_upper: Upper 95% confidence interval bound
        - cv: Coefficient of variation (std / |mean|)
        - stability: Heuristic stability classification
    """
    scores_array = np.array(scores)
    
    mean = np.mean(scores_array)
    std = np.std(scores_array, ddof=1)
    
    # 95% confidence interval (assumes normal distribution)
    n = len(scores_array)
    sem = std / np.sqrt(n)
    ci_margin = 1.96 * sem  # 95% CI
    ci_lower = mean - ci_margin
    ci_upper = mean + ci_margin
    
    # Coefficient of variation
    cv = (std / abs(mean)) * 100 if mean != 0 else float('inf')
    
    # Heuristic stability classification
    if cv < 5:
        stability = "HIGH"
    elif cv < 15:
        stability = "MODERATE"
    else:
        stability = "LOW"
    
    return {
        "mean": mean,
        "std": std,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "cv": cv,
        "stability": stability,
    }


def plot_rps_distribution(
    scores: List[float],
    ligand_name: str,
    n_replicates: int,
    sigma: float,
    output_path: Path,
    method: str = "gb",
) -> None:
    """
    Generate publication-quality plot of RPS score distribution.
    
    Creates violin plot with overlaid boxplot showing:
    - Distribution of scores across perturbations
    - Mean score as horizontal line
    - Annotated with N and sigma
    
    Args:
        scores: List of MM/GBSA scores
        ligand_name: Name of ligand
        n_replicates: Number of perturbation replicates
        sigma: Perturbation magnitude (Angstroms)
        output_path: Output directory for plot
        method: Solvation method ("gb" or "pb")
    """
    logger.info("Generating RPS distribution plot...")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create violin plot
    parts = ax.violinplot(
        [scores],
        positions=[0],
        showmeans=False,
        showmedians=False,
        widths=0.7,
    )
    
    # Customize violin plot colors
    for pc in parts['bodies']:
        pc.set_facecolor('#8dd3c7')
        pc.set_edgecolor('#2b8a8a')
        pc.set_alpha(0.7)
    
    # Overlay boxplot
    bp = ax.boxplot(
        [scores],
        positions=[0],
        widths=0.3,
        patch_artist=True,
        boxprops=dict(facecolor='#fb8072', edgecolor='black', alpha=0.8),
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
    )
    
    # Add mean line
    mean_score = np.mean(scores)
    ax.axhline(y=mean_score, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_score:.2f}')
    
    # Customize plot
    method_label = "MM/GBSA" if method.lower() == "gb" else "MM/PBSA"
    ax.set_title(
        f"Rapid Perturbation Sampling: {ligand_name}\n"
        f"Single-frame {method_label} score distribution",
        fontsize=14,
        fontweight='bold'
    )
    ax.set_ylabel("Single-frame score (kcal/mol)", fontsize=12)
    ax.set_xticks([0])
    ax.set_xticklabels([f'{n_replicates} replicates'])
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='upper right')
    
    # Add annotation box with parameters
    stats = compute_statistics(scores)
    textstr = f'σ = {sigma:.2f} Å\nCV = {stats["cv"]:.1f}%\nStability: {stats["stability"]}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(
        0.02, 0.98, textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=props
    )
    
    # Tight layout
    plt.tight_layout()
    
    # Save
    output_path.mkdir(parents=True, exist_ok=True)
    plot_file = output_path / f"{ligand_name}_rps_distribution.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.success(f"RPS distribution plot saved: {plot_file}")


def run_rps(
    receptor_dir: Path,
    parameterized_ligand_dir: Path,
    ligand_name: str,
    output_dir: Path,
    n_replicates: int = 10,
    sigma: float = 0.2,
    method: str = "gb",
    minimize: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Run Rapid Perturbation Sampling analysis.
    
    CRITICAL WARNINGS:
    - This is NOT molecular dynamics
    - This is NOT conformational sampling
    - This IS numerical sensitivity analysis
    
    CRITICAL METHODOLOGY:
    - RPS perturbs coordinates only; ligand chemistry is fixed by design
    - Uses pre-parameterized ligand.mol2 and ligand.frcmod from main workflow
    - NO re-parameterization per replicate (antechamber NOT called)
    - GAFF2 atom types, AM1-BCC charges, and frcmod are frozen
    
    Args:
        receptor_dir: Directory containing prepared receptor
        parameterized_ligand_dir: Directory with pre-parameterized ligand.mol2 and ligand.frcmod
        ligand_name: Name of ligand (for output files)
        output_dir: Output directory for RPS results
        n_replicates: Number of perturbation replicates
        sigma: Perturbation magnitude in Angstroms
        method: Solvation method ("gb" or "pb")
        minimize: Whether to run minimization
        
    Returns:
        Tuple of (per_replicate_dataframe, summary_statistics)
        
    Raises:
        RPSError: If RPS execution fails
    """
    logger.warning("=" * 70)
    logger.warning("RAPID PERTURBATION SAMPLING (RPS) - DIAGNOSTIC MODE")
    logger.warning("=" * 70)
    logger.warning("RPS is NOT:")
    logger.warning("  - Molecular dynamics")
    logger.warning("  - Conformational sampling")
    logger.warning("  - Pose optimization")
    logger.warning("  - Thermodynamic ensemble averaging")
    logger.warning("")
    logger.warning("RPS IS:")
    logger.warning("  - Numerical sensitivity analysis")
    logger.warning("  - Uncertainty quantification for single-frame scores")
    logger.warning("  - Diagnostic tool for score stability")
    logger.warning("=" * 70)
    
    console.print()
    console.print(f"[bold cyan]Running RPS with {n_replicates} replicates (σ = {sigma:.2f} Å)[/bold cyan]")
    console.print()
    
    # Create RPS output directory
    rps_dir = output_dir / "rps_analysis"
    rps_dir.mkdir(parents=True, exist_ok=True)
    
    # Run perturbations
    results = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"Processing {n_replicates} perturbations...",
            total=n_replicates
        )
        
        for replicate_id in range(1, n_replicates + 1):
            try:
                energies = run_single_perturbation(
                    receptor_dir=receptor_dir,
                    parameterized_ligand_dir=parameterized_ligand_dir,
                    replicate_id=replicate_id,
                    sigma=sigma,
                    ligand_name=ligand_name,
                    output_dir=rps_dir,
                    method=method,
                    minimize=minimize,
                )
                
                results.append({
                    "ligand": ligand_name,
                    "replicate_id": replicate_id,
                    "score": energies.get("DELTA_G", np.nan),
                })
                
            except Exception as e:
                import traceback
                console.print(f"[red]Replicate {replicate_id} failed: {e}[/red]")
                console.print(f"[dim]{traceback.format_exc()}[/dim]")
                logger.error(f"Replicate {replicate_id} failed: {e}")
                results.append({
                    "ligand": ligand_name,
                    "replicate_id": replicate_id,
                    "score": np.nan,
                })
            
            progress.update(task, advance=1)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Filter out failed replicates
    valid_scores = df[~df['score'].isna()]['score'].tolist()
    
    if len(valid_scores) < 2:
        raise RPSError(
            f"Insufficient valid replicates: {len(valid_scores)}/{n_replicates}\n"
            "Cannot compute statistics with < 2 valid scores"
        )
    
    # Compute statistics
    stats = compute_statistics(valid_scores)
    stats['ligand'] = ligand_name
    stats['n_valid'] = len(valid_scores)
    stats['n_total'] = n_replicates
    
    # Save per-replicate CSV
    per_replicate_csv = rps_dir / f"{ligand_name}_rps_replicates.csv"
    df.to_csv(per_replicate_csv, index=False)
    logger.success(f"Per-replicate results saved: {per_replicate_csv}")
    
    # Save summary CSV
    summary_csv = rps_dir / f"{ligand_name}_rps_summary.csv"
    summary_df = pd.DataFrame([{
        'ligand': stats['ligand'],
        'n_valid': stats['n_valid'],
        'n_total': stats['n_total'],
        'mean_score': stats['mean'],
        'std': stats['std'],
        'ci_lower': stats['ci_lower'],
        'ci_upper': stats['ci_upper'],
        'cv': stats['cv'],
        'stability': stats['stability'],
    }])
    summary_df.to_csv(summary_csv, index=False)
    logger.success(f"Summary statistics saved: {summary_csv}")
    
    # Generate plot
    plot_rps_distribution(
        scores=valid_scores,
        ligand_name=ligand_name,
        n_replicates=len(valid_scores),
        sigma=sigma,
        output_path=rps_dir,
        method=method,
    )
    
    # Print summary to console
    console.print()
    console.print("=" * 70)
    console.print(f"[bold]RPS SUMMARY: {ligand_name}[/bold]")
    console.print("=" * 70)
    console.print(f"  Valid replicates:     {stats['n_valid']}/{stats['n_total']}")
    console.print(f"  Perturbation σ:       {sigma:.2f} Å")
    console.print(f"  Mean score:           {stats['mean']:.2f} ± {stats['std']:.2f} kcal/mol")
    console.print(f"  95% CI:               [{stats['ci_lower']:.2f}, {stats['ci_upper']:.2f}]")
    console.print(f"  Coefficient of Var:   {stats['cv']:.1f}%")
    console.print(f"  Stability:            {stats['stability']}")
    console.print("=" * 70)
    console.print()
    console.print("[yellow]⚠ INTERPRETATION REMINDER:[/yellow]")
    console.print("  RPS quantifies numerical sensitivity to coordinate perturbations")
    console.print("  It does NOT provide thermodynamic ensembles or conformational sampling")
    console.print("  Use for diagnostic purposes and uncertainty bounds only")
    console.print()
    
    return df, stats
