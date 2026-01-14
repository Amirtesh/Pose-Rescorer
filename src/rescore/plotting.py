"""
Plotting utilities for MM/GBSA and MM/PBSA score visualization.

CRITICAL TERMINOLOGY RULE:
- Use "score" everywhere, NEVER "energy"
- Plots are descriptive tools for result inspection
- They do not alter score interpretation

This module provides optional plotting functionality for:
1. Single ligand score components (rescore integrate)
2. Multi-ligand score comparison (rescore batch)
"""

from pathlib import Path
from typing import Dict, List, Any

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless systems
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger


def plot_integrate_scores(
    energies: Dict[str, float],
    ligand_name: str,
    output_path: Path,
    method: str = "gb",
) -> None:
    """
    Generate bar chart of score components for a single ligand.
    
    Creates a 3-bar chart showing:
    - MM score (ΔH)
    - Solvation score (ΔG_GB or ΔG_PB)
    - Total score (ΔG_bind)
    
    Args:
        energies: Dictionary of energy components from parse_rescore_results()
        ligand_name: Name of the ligand for filename
        output_path: Output directory for the plot
        method: Solvation method ("gb" or "pb")
        
    Saves:
        <ligand_name>_score_components.png in output_path
    """
    logger.info(f"Generating score component plot for {ligand_name}...")
    
    # Extract values with fallbacks
    mm_score = energies.get("DELTA_H", 0.0)
    solvation_score = energies.get("DELTA_G_GB", 0.0)
    total_score = energies.get("DELTA_G", 0.0)
    
    # Data for plotting
    components = ["MM score", "Solvation score", "Total score"]
    scores = [mm_score, solvation_score, total_score]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Generate distinct colors using seaborn palette
    colors = sns.color_palette("Set2", n_colors=3)
    
    # Create bar chart
    bars = ax.bar(components, scores, color=colors, edgecolor='black', linewidth=1.2)
    
    # Customize appearance
    method_label = "MM/GBSA" if method.lower() == "gb" else "MM/PBSA"
    ax.set_title(f"Single-frame {method_label} score components (descriptive)", fontsize=14, fontweight='bold')
    ax.set_xlabel("Score component", fontsize=12)
    ax.set_ylabel("Relative score (single-frame, rank-oriented)", fontsize=12)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, score in zip(bars, scores):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.,
            height,
            f'{score:.2f}',
            ha='center',
            va='bottom' if height >= 0 else 'top',
            fontsize=10,
            fontweight='bold'
        )
    
    # Tight layout
    plt.tight_layout()
    
    # Save to file
    output_path.mkdir(parents=True, exist_ok=True)
    plot_file = output_path / f"{ligand_name}_score_components.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.success(f"Score component plot saved: {plot_file}")


def plot_batch_scores(
    results: List[Dict[str, Any]],
    output_path: Path,
    method: str = "gb",
) -> None:
    """
    Generate bar chart comparing total scores across multiple ligands.
    
    Creates a multi-bar chart showing only total MM/GBSA or MM/PBSA scores,
    sorted from best (most negative) to worst.
    
    Args:
        results: List of result dictionaries from run_batch_rescore()
        output_path: Output directory for the plot
        method: Solvation method ("gb" or "pb")
        
    Saves:
        batch_relative_scores.png in output_path
    """
    logger.info("Generating batch score comparison plot...")
    
    # Filter successful results only
    successful = [r for r in results if r["status"] == "success" and r["delta_g_bind"] is not None]
    
    if not successful:
        logger.warning("No successful results to plot")
        return
    
    # Sort by total score (most favorable first)
    sorted_results = sorted(successful, key=lambda x: x["delta_g_bind"])
    
    # Extract data
    ligand_names = [r["ligand_name"] for r in sorted_results]
    total_scores = [r["delta_g_bind"] for r in sorted_results]
    
    # Create figure with dynamic width based on number of ligands
    fig_width = max(8, len(ligand_names) * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    
    # Generate distinct colors using seaborn palette
    n_colors = len(ligand_names)
    if n_colors <= 10:
        colors = sns.color_palette("Set3", n_colors=n_colors)
    else:
        colors = sns.color_palette("husl", n_colors=n_colors)
    
    # Create bar chart
    bars = ax.bar(range(len(ligand_names)), total_scores, color=colors, edgecolor='black', linewidth=1.2)
    
    # Customize appearance
    method_label = "MM/GBSA" if method.lower() == "gb" else "MM/PBSA"
    ax.set_title(f"Relative single-frame {method_label} scores", fontsize=14, fontweight='bold')
    ax.set_xlabel("Ligand", fontsize=12)
    ax.set_ylabel("Relative score (single-frame, rank-oriented)", fontsize=12)
    ax.set_xticks(range(len(ligand_names)))
    ax.set_xticklabels(ligand_names, rotation=45, ha='right')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, score in zip(bars, total_scores):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.,
            height,
            f'{score:.2f}',
            ha='center',
            va='bottom' if height >= 0 else 'top',
            fontsize=9,
            fontweight='bold',
            rotation=0
        )
    
    # Tight layout
    plt.tight_layout()
    
    # Save to file
    output_path.mkdir(parents=True, exist_ok=True)
    plot_file = output_path / "batch_relative_scores.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.success(f"Batch score comparison plot saved: {plot_file}")
