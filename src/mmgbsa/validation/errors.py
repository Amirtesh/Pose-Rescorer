"""
Custom exceptions for MM/GBSA validation.

All validation errors inherit from MMGBSAValidationError to allow
catch-all error handling while maintaining specific error types
for different failure modes.
"""


class MMGBSAValidationError(Exception):
    """Base exception for all MM/GBSA validation failures."""
    pass


class MissingBackboneAtomsError(MMGBSAValidationError):
    """Raised when PDB structure has missing backbone atoms (N, CA, C, O) in standard residues."""
    pass


class NoLigandDetectedError(MMGBSAValidationError):
    """Raised when no ligand molecules are detected in the structure."""
    pass


class MultipleLigandsError(MMGBSAValidationError):
    """Raised when PDB structure contains multiple ligand molecules."""
    pass


class MetalDetectedError(MMGBSAValidationError):
    """Raised when PDB structure contains metal atoms (not supported by GAFF2)."""
    pass


class AltlocOccupancyError(MMGBSAValidationError):
    """Raised when PDB structure contains alternate locations or partial occupancy."""
    pass


class UnsupportedSystemError(MMGBSAValidationError):
    """Raised for other unsupported system types (e.g., covalent ligands, DNA/RNA)."""
    pass
