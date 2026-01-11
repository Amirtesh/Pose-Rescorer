"""
Tests for PDB validation module.

These tests verify that validation checks properly detect:
- Missing backbone atoms (N, CA, C, O only)
- Metal atoms
- Zero ligands
- Exactly one ligand
- Multiple ligands
"""

import pytest
from pathlib import Path
from io import StringIO
from Bio.PDB import PDBParser

from mmgbsa.validation import (
    check_missing_backbone_atoms,
    check_for_metals,
    check_altloc_and_occupancy,
    detect_ligands,
    MissingBackboneAtomsError,
    NoLigandDetectedError,
    MetalDetectedError,
    MultipleLigandsError,
    AltlocOccupancyError,
)


def create_mock_structure(pdb_string: str) -> "Structure":
    """Helper to create BioPython Structure from PDB string."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure("test", StringIO(pdb_string))


def test_check_for_metals_detects_zinc() -> None:
    """Test that zinc atoms are properly detected."""
    pdb_with_zinc = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5 ZN   ZN  B   1       5.000   5.000   5.000  1.00  0.00          ZN
END
"""
    structure = create_mock_structure(pdb_with_zinc)

    with pytest.raises(MetalDetectedError) as exc_info:
        check_for_metals(structure)

    assert "ZN" in str(exc_info.value)
    assert "Metal coordination is not supported" in str(exc_info.value)


def test_check_for_metals_passes_clean_structure() -> None:
    """Test that structures without metals pass validation."""
    clean_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
END
"""
    structure = create_mock_structure(clean_pdb)
    # Should not raise
    check_for_metals(structure)


def test_detect_ligands_raises_on_zero_ligands() -> None:
    """Test that structures with no ligands raise NoLigandDetectedError."""
    pdb_protein_only = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       3.321   1.525   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       4.100   2.750   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       5.600   2.500   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.100   1.400   0.000  1.00  0.00           O
HETATM    9  O   HOH A 300       8.000   8.000   8.000  1.00  0.00           O
END
"""
    structure = create_mock_structure(pdb_protein_only)

    with pytest.raises(NoLigandDetectedError) as exc_info:
        detect_ligands(structure)

    assert "No ligand molecules detected" in str(exc_info.value)


def test_detect_ligands_detects_single_ligand() -> None:
    """Test that single ligand is correctly detected and returned."""
    pdb_with_one_ligand = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    6  C2  LIG A 200       6.000   5.000   5.000  1.00  0.00           C
HETATM    7  O   HOH A 300       8.000   8.000   8.000  1.00  0.00           O
END
"""
    structure = create_mock_structure(pdb_with_one_ligand)
    ligand_resname = detect_ligands(structure)
    assert ligand_resname == "LIG"


def test_detect_ligands_raises_on_multiple_ligands() -> None:
    """Test that multiple ligand molecules raise MultipleLigandsError."""
    pdb_with_two_ligands = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    6  C2  LIG A 200       6.000   5.000   5.000  1.00  0.00           C
HETATM    7  C1  LIG B 200      10.000  10.000  10.000  1.00  0.00           C
HETATM    8  C2  LIG B 200      11.000  10.000  10.000  1.00  0.00           C
END
"""
    structure = create_mock_structure(pdb_with_two_ligands)

    with pytest.raises(MultipleLigandsError) as exc_info:
        detect_ligands(structure)

    assert "Detected 2 ligand molecule(s)" in str(exc_info.value)
    assert "Only single-ligand systems are supported" in str(exc_info.value)


def test_detect_ligands_raises_on_different_ligand_types() -> None:
    """Test that multiple different ligand types raise MultipleLigandsError."""
    pdb_with_mixed_ligands = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    6  C1  INH B 200      10.000  10.000  10.000  1.00  0.00           C
END
"""
    structure = create_mock_structure(pdb_with_mixed_ligands)

    with pytest.raises(MultipleLigandsError) as exc_info:
        detect_ligands(structure)

    assert "Detected 2 ligand molecule(s)" in str(exc_info.value)
    assert "LIG" in str(exc_info.value)
    assert "INH" in str(exc_info.value)


def test_detect_ligands_ignores_water() -> None:
    """Test that water molecules are correctly ignored."""
    pdb_with_ligand_and_water = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    6  O   HOH A 301       8.000   8.000   8.000  1.00  0.00           O
HETATM    7  O   HOH A 302       9.000   9.000   9.000  1.00  0.00           O
HETATM    8  O   WAT A 303      10.000  10.000  10.000  1.00  0.00           O
END
"""
    structure = create_mock_structure(pdb_with_ligand_and_water)
    ligand_resname = detect_ligands(structure)
    assert ligand_resname == "LIG"


def test_detect_ligands_ignores_ions() -> None:
    """Test that common ions are correctly ignored."""
    pdb_with_ligand_and_ions = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
HETATM    5  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    6 NA   NA  B 401       8.000   8.000   8.000  1.00  0.00          NA
HETATM    7 CL   CL  B 402       9.000   9.000   9.000  1.00  0.00          CL
END
"""
    structure = create_mock_structure(pdb_with_ligand_and_ions)
    ligand_resname = detect_ligands(structure)
    assert ligand_resname == "LIG"


def test_check_missing_backbone_atoms_detects_incomplete_backbone() -> None:
    """Test that missing backbone atoms (N, CA, C, O) are detected."""
    # ALA missing CA atom
    incomplete_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  O   ALA A   1       1.251   2.391   0.000  1.00  0.00           O
END
"""
    structure = create_mock_structure(incomplete_pdb)

    with pytest.raises(MissingBackboneAtomsError) as exc_info:
        check_missing_backbone_atoms(structure)

    assert "missing backbone atoms" in str(exc_info.value).lower()
    assert "ALA" in str(exc_info.value)
    assert "Side-chain atoms are not validated" in str(exc_info.value)


def test_check_backbone_allows_c_terminal_with_ot_atoms() -> None:
    """Test that C-terminal residues with OT1/OT2 instead of O are accepted."""
    # Complete structure with C-terminal using OT1/OT2
    c_terminal_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  1.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
ATOM      5  CB  ALA A   1       2.000  -0.700   1.250  1.00  1.00           C
ATOM      6  N   GLY A   2       3.321   1.525   0.000  1.00  1.00           N
ATOM      7  CA  GLY A   2       4.100   2.750   0.000  1.00  1.00           C
ATOM      8  C   GLY A   2       5.600   2.500   0.000  1.00  1.00           C
ATOM      9  OT1 GLY A   2       6.100   1.400   0.000  1.00  1.00           O
ATOM     10  OT2 GLY A   2       6.300   3.500   0.000  1.00  1.00           O
END
"""
    structure = create_mock_structure(c_terminal_pdb)
    # Should not raise
    check_missing_backbone_atoms(structure)


def test_check_backbone_rejects_c_terminal_missing_both_o_types() -> None:
    """Test that C-terminal missing both O and OT1/OT2 is rejected."""
    # C-terminal missing all oxygen atoms
    bad_c_terminal_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  1.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
ATOM      5  N   GLY A   2       3.321   1.525   0.000  1.00  1.00           N
ATOM      6  CA  GLY A   2       4.100   2.750   0.000  1.00  1.00           C
ATOM      7  C   GLY A   2       5.600   2.500   0.000  1.00  1.00           C
END
"""
    structure = create_mock_structure(bad_c_terminal_pdb)

    with pytest.raises(MissingBackboneAtomsError) as exc_info:
        check_missing_backbone_atoms(structure)

    assert "C-terminal" in str(exc_info.value)
    assert "GLY" in str(exc_info.value)


def test_check_backbone_rejects_c_terminal_with_only_ot1() -> None:
    """Test that C-terminal with only OT1 (missing OT2) is rejected."""
    # C-terminal with incomplete terminal oxygens
    incomplete_ot_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  1.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
ATOM      5  N   GLY A   2       3.321   1.525   0.000  1.00  1.00           N
ATOM      6  CA  GLY A   2       4.100   2.750   0.000  1.00  1.00           C
ATOM      7  C   GLY A   2       5.600   2.500   0.000  1.00  1.00           C
ATOM      8  OT1 GLY A   2       6.100   1.400   0.000  1.00  1.00           O
END
"""
    structure = create_mock_structure(incomplete_ot_pdb)

    with pytest.raises(MissingBackboneAtomsError) as exc_info:
        check_missing_backbone_atoms(structure)

    assert "C-terminal" in str(exc_info.value)
    assert "OT2" in str(exc_info.value) or "missing O" in str(exc_info.value).lower()


def test_check_backbone_requires_o_for_internal_residues() -> None:
    """Test that internal residues must have O (OT1/OT2 not accepted)."""
    # Internal residue with OT1/OT2 instead of O (invalid)
    internal_ot_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  1.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  OT1 ALA A   1       1.251   2.391   0.000  1.00  1.00           O
ATOM      5  OT2 ALA A   1       2.500   2.500   0.000  1.00  1.00           O
ATOM      6  N   GLY A   2       3.321   1.525   0.000  1.00  1.00           N
ATOM      7  CA  GLY A   2       4.100   2.750   0.000  1.00  1.00           C
ATOM      8  C   GLY A   2       5.600   2.500   0.000  1.00  1.00           C
ATOM      9  O   GLY A   2       6.100   1.400   0.000  1.00  1.00           O
END
"""
    structure = create_mock_structure(internal_ot_pdb)

    with pytest.raises(MissingBackboneAtomsError) as exc_info:
        check_missing_backbone_atoms(structure)

    assert "internal" in str(exc_info.value).lower()
    assert "ALA" in str(exc_info.value)


def test_check_altloc_detects_alternate_locations() -> None:
    """Test that alternate locations are properly detected."""
    # PDB with alternate locations (A and B conformers)
    pdb_with_altloc = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.50  1.00           C
ATOM      3  CA BALA A   1       1.500   0.100   0.000  0.50  1.00           C
ATOM      4  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      5  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
HETATM    6  C1  LIG B 100      10.000  10.000  10.000  1.00  1.00           C
END
"""
    structure = create_mock_structure(pdb_with_altloc)

    with pytest.raises(AltlocOccupancyError) as exc_info:
        check_altloc_and_occupancy(structure)

    assert "alternate locations" in str(exc_info.value).lower()
    assert "altloc 'A'" in str(exc_info.value) or "altloc 'B'" in str(exc_info.value)
    assert "AmberTools requires single-conformation" in str(exc_info.value)


def test_check_altloc_detects_partial_occupancy() -> None:
    """Test that partial occupancy is properly detected."""
    # PDB with partial occupancy
    pdb_with_partial_occ = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  0.60  0.60           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
HETATM    5  C1  LIG B 100      10.000  10.000  10.000  1.00  0.75           C
END
"""
    structure = create_mock_structure(pdb_with_partial_occ)

    with pytest.raises(AltlocOccupancyError) as exc_info:
        check_altloc_and_occupancy(structure)

    assert "partial occupancy" in str(exc_info.value).lower()
    assert "occupancy 0.60" in str(exc_info.value) or "occupancy 0.75" in str(exc_info.value)
    assert "full occupancy" in str(exc_info.value).lower()


def test_check_altloc_passes_clean_structure() -> None:
    """Test that clean structures with no altloc or partial occupancy pass."""
    clean_pdb = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  1.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      4  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
HETATM    5  C1  LIG B 100      10.000  10.000  10.000  1.00  1.00           C
END
"""
    structure = create_mock_structure(clean_pdb)
    # Should not raise
    check_altloc_and_occupancy(structure)


def test_check_altloc_detects_both_issues() -> None:
    """Test that both altloc and partial occupancy are reported together."""
    pdb_with_both = """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.80           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.50  0.50           C
ATOM      3  CA BALA A   1       1.500   0.100   0.000  0.50  0.50           C
ATOM      4  C   ALA A   1       2.009   1.420   0.000  1.00  1.00           C
ATOM      5  O   ALA A   1       1.251   2.391   0.000  1.00  1.00           O
HETATM    6  C1  LIG B 100      10.000  10.000  10.000  1.00  1.00           C
END
"""
    structure = create_mock_structure(pdb_with_both)

    with pytest.raises(AltlocOccupancyError) as exc_info:
        check_altloc_and_occupancy(structure)

    error_msg = str(exc_info.value)
    assert "alternate locations" in error_msg.lower()
    assert "partial occupancy" in error_msg.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
