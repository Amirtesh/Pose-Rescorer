"""Parameterization exception classes for rescoring workflow."""

class ParameterizationError(Exception):
    """Base class for parameterization-related errors."""


class AntechamberError(ParameterizationError):
    """Errors raised when antechamber fails."""


class ParmchkError(ParameterizationError):
    """Errors raised when parmchk2 fails."""


class LigandExtractionError(ParameterizationError):
    """Errors raised when ligand extraction or file handling fails."""


class ChemistryValidationError(ParameterizationError):
    """Errors raised when ligand chemistry validation fails."""
