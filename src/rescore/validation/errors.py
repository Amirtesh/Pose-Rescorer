"""Validation exception classes for rescoring workflow."""

class MMGBSAValidationError(Exception):
    """Base class for validation errors."""


class MissingBackboneAtomsError(MMGBSAValidationError):
    """Backbone atoms are missing in the protein structure."""


class NoLigandDetectedError(MMGBSAValidationError):
    """No ligand detected in the structure when one was expected."""


class MultipleLigandsError(MMGBSAValidationError):
    """Multiple ligands detected when only one is supported."""


class MetalDetectedError(MMGBSAValidationError):
    """Metal ions detected (not supported in this workflow)."""


class AltlocOccupancyError(MMGBSAValidationError):
    """Alternate locations or partial occupancies detected."""


class UnsupportedSystemError(MMGBSAValidationError):
    """Unsupported system composition detected."""
