"""
Input validation module for MM/GBSA calculations.

This module provides functions to validate protein-ligand complexes
before attempting MM/GBSA calculations. All checks fail fast with
explicit error messages.
"""

from mmgbsa.validation.errors import (
    MMGBSAValidationError,
    MissingBackboneAtomsError,
    NoLigandDetectedError,
    MultipleLigandsError,
    MetalDetectedError,
    AltlocOccupancyError,
    UnsupportedSystemError,
)
from mmgbsa.validation.pdb_checks import (
    validate_pdb_complex,
    check_missing_backbone_atoms,
    check_for_metals,
    check_altloc_and_occupancy,
    detect_ligands,
)

__all__ = [
    "MMGBSAValidationError",
    "MissingBackboneAtomsError",
    "NoLigandDetectedError",
    "MultipleLigandsError",
    "MetalDetectedError",
    "AltlocOccupancyError",
    "UnsupportedSystemError",
    "validate_pdb_complex",
    "check_missing_backbone_atoms",
    "check_for_metals",
    "check_altloc_and_occupancy",
    "detect_ligands",
]
