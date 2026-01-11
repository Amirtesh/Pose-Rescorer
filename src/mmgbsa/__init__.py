"""
MM/GBSA rescoring package for protein-ligand complexes using AmberTools.

This package provides a minimal, defensible workflow for MM/GBSA calculations
using only AmberTools executables (tleap, antechamber, parmchk2, sander, MMPBSA.py).

Key limitations:
- Single-structure or small ensemble only
- No entropy calculations
- Relative ranking only (not absolute binding free energies)
- No support for metals, covalent ligands, or multiple ligands
"""

__version__ = "0.1.0"
__author__ = "AmberTools MM/GBSA Team"

from mmgbsa.validation.errors import (
    MMGBSAValidationError,
    MissingBackboneAtomsError,
    NoLigandDetectedError,
    MultipleLigandsError,
    MetalDetectedError,
    AltlocOccupancyError,
    UnsupportedSystemError,
)

__all__ = [
    "__version__",
    "__author__",
    "MMGBSAValidationError",
    "MissingBackboneAtomsError",
    "NoLigandDetectedError",
    "MultipleLigandsError",
    "MetalDetectedError",
    "AltlocOccupancyError",
    "UnsupportedSystemError",
]
