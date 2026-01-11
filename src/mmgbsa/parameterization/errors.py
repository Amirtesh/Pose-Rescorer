"""
Custom exceptions for parameterization failures.
"""


class ParameterizationError(Exception):
    """Base exception for all parameterization errors."""
    pass


class AntechamberError(ParameterizationError):
    """Raised when antechamber fails."""
    pass


class ParmchkError(ParameterizationError):
    """Raised when parmchk2 fails."""
    pass


class LigandExtractionError(ParameterizationError):
    """Raised when ligand extraction fails."""
    pass


class ChemistryValidationError(ParameterizationError):
    """Raised when ligand chemistry validation fails (invalid electron count, missing H, etc.)."""
    pass
