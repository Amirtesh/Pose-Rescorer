"""
Configuration constants for MM/GBSA calculations.

This module defines force field parameters, metal element lists,
and other constants used throughout the package.
"""

from typing import Final, Set

# Force field parameters (non-negotiable)
PROTEIN_FF: Final[str] = "ff14SB"
LIGAND_FF: Final[str] = "gaff2"
CHARGE_METHOD: Final[str] = "bcc"  # AM1-BCC

# Metal elements that are not supported
# Covers common metals in PDB structures (transition metals, alkali, alkaline earth)
METAL_ELEMENTS: Final[Set[str]] = {
    # Transition metals
    "FE", "CU", "ZN", "MN", "CO", "NI", "MO", "W", "V", "CR",
    "PT", "PD", "AU", "AG", "CD", "HG", "RU", "RH", "IR", "OS",
    # Alkali and alkaline earth
    "NA", "K", "CA", "MG", "LI", "RB", "CS", "SR", "BA",
    # Lanthanides and actinides (occasionally seen)
    "LA", "CE", "GD", "YB", "EU", "TB", "U", "TH",
}

# Standard amino acid residue names
STANDARD_AMINO_ACIDS: Final[Set[str]] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    # Common protonation/modification variants
    "HIE", "HID", "HIP",  # Histidine (epsilon, delta, doubly protonated)
    "HSE", "HSD", "HSP",  # Histidine alternative naming
    "CYX",  # Disulfide-bonded cysteine
}

# Common solvent and ion residue names to ignore during validation
SOLVENT_RESIDUES: Final[Set[str]] = {
    "HOH", "WAT", "H2O", "TIP3", "TIP4", "SPC", "SOL",
    "CL", "NA", "K", "MG", "CA",  # Simple ions
}

# File naming conventions
DEFAULT_PROTEIN_NAME: Final[str] = "protein"
DEFAULT_LIGAND_NAME: Final[str] = "ligand"
DEFAULT_COMPLEX_NAME: Final[str] = "complex"

# MM/GBSA calculation defaults
DEFAULT_GB_MODEL: Final[int] = 2  # GB model (1=GB-HCT, 2=GB-OBC1, 5=GB-OBC2)
DEFAULT_SALT_CONCENTRATION: Final[float] = 0.15  # M (physiological)
