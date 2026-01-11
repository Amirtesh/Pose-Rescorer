"""
Ligand extraction and GAFF2 parameterization using AmberTools.

This module handles:
1. Extracting ligand molecules from validated PDB complexes
2. Running antechamber for AM1-BCC charge calculation and GAFF2 typing
3. Running parmchk2 to generate missing force field parameters
"""

import subprocess
from pathlib import Path
from typing import Optional

from Bio.PDB import PDBIO, Structure
from loguru import logger

from rescore.config import LIGAND_FF, CHARGE_METHOD
from rescore.parameterization.errors import (
    AntechamberError,
    ParmchkError,
    LigandExtractionError,
    ChemistryValidationError,
)
from rescore.parameterization.connectivity import (
    extract_ligand_simple,
    extract_ligand_with_connectivity,
)


# Valence electrons for common elements
VALENCE_ELECTRONS = {
    "H": 1, "C": 4, "N": 5, "O": 6, "F": 7,
    "P": 5, "S": 6, "CL": 7, "BR": 7, "I": 7,
}


def extract_ligand(
    structure: Structure,
    ligand_resname: str,
    output_path: Path,
    source_pdb_path: Optional[Path] = None,
    preserve_connectivity: bool = True,
) -> None:
    """
    Extract ligand residue(s) from a PDB structure and write to a standalone PDB file.

    Extracts all ATOM/HETATM records matching the specified ligand residue name.
    Preserves atom names, coordinates, chain IDs, and residue numbering.

    NOTE: Due to BioPython PDB parser limitations with certain PDB formats,
    this function falls back to simple text extraction when the source PDB path
    is provided.

    Args:
        structure: BioPython Structure object containing the ligand
        ligand_resname: Three-letter residue name of the ligand to extract
        output_path: Path where the ligand PDB file will be written
        source_pdb_path: Optional path to source PDB file (enables connectivity preservation)
        preserve_connectivity: If True and source_pdb_path provided, extract with CONECT records

    Raises:
        LigandExtractionError: If no ligand residues found or file write fails
    """
    logger.info(f"Extracting ligand '{ligand_resname}' from structure")

    # If source PDB provided and connectivity preservation requested, use text-based extraction
    if source_pdb_path and preserve_connectivity:
        logger.debug("Attempting connectivity-preserving extraction")
        try:
            extract_ligand_with_connectivity(source_pdb_path, ligand_resname, output_path)
            return
        except LigandExtractionError as e:
            # CONECT records missing - fall back to simple text extraction
            logger.warning(f"Connectivity extraction failed: {e}")
            logger.info("Falling back to simple text extraction (no CONECT records)")
            try:
                extract_ligand_simple(source_pdb_path, ligand_resname, output_path)
                return
            except LigandExtractionError as e2:
                logger.error(f"Simple text extraction also failed: {e2}")
                logger.info("Falling back to BioPython extraction as last resort")

    # BioPython-based extraction (no connectivity preservation)
    # This is the last resort fallback and may have parsing issues
    logger.debug("Using BioPython-based extraction")
    ligand_atom_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_resname:
                    ligand_atom_count += len(list(residue.get_atoms()))

    if ligand_atom_count == 0:
        raise LigandExtractionError(
            f"No residues with name '{ligand_resname}' found in structure. "
            "Ensure the ligand residue name matches the validated ligand."
        )

    logger.debug(f"Found ligand with {ligand_atom_count} atoms")

    class LigandSelect:
        def __init__(self, target_resname: str):
            self.target_resname = target_resname

        def accept_model(self, model):
            return True

        def accept_chain(self, chain):
            return True

        def accept_residue(self, residue):
            return residue.get_resname() == self.target_resname

        def accept_atom(self, atom):
            return True

    try:
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), select=LigandSelect(ligand_resname))
        
        # Verify output
        if not output_path.exists():
            raise LigandExtractionError(f"Failed to create output file: {output_path}")
        
        with open(output_path) as f:
            atom_lines = sum(1 for line in f if line.startswith(("ATOM", "HETATM")))
        
        if atom_lines == 0:
            raise LigandExtractionError(
                f"Ligand PDB file is empty. BioPython failed to extract ligand atoms."
            )
        
        logger.success(f"Ligand extracted to: {output_path} ({atom_lines} atoms)")
        
    except Exception as e:
        raise LigandExtractionError(
            f"Failed to write ligand PDB file to {output_path}: {e}"
        )


def validate_ligand_chemistry(
    ligand_pdb: Path,
    net_charge: Optional[int] = None,
) -> int:
    """
    Validate ligand chemistry before antechamber: check electron count and hydrogens.

    Performs strict validation without attempting to fix chemistry:
    1. Counts atoms by element from PDB file
    2. Estimates total valence electrons
    3. Checks if electron count is valid for given charge
    4. Warns about missing hydrogens (common issue)

    Args:
        ligand_pdb: Path to ligand PDB file
        net_charge: Expected net charge (if None, tries to infer valid charge)

    Returns:
        Recommended net charge if validation passes

    Raises:
        ChemistryValidationError: If structure has invalid electron count or other issues
    """
    logger.info("Validating ligand chemistry before parameterization")

    # Parse PDB file and count atoms by element
    element_counts = {}
    total_atoms = 0
    hydrogen_count = 0

    try:
        with open(ligand_pdb) as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue

                total_atoms += 1
                # Element is in columns 77-78 (right-justified, may have leading space)
                element = line[76:78].strip().upper()
                
                if not element:
                    # Fallback: try to infer from atom name (columns 13-16)
                    atom_name = line[12:16].strip()
                    element = atom_name[0].upper()
                    logger.debug(f"No element column, inferred '{element}' from atom name '{atom_name}'")

                element_counts[element] = element_counts.get(element, 0) + 1
                
                if element == "H":
                    hydrogen_count += 1

    except Exception as e:
        raise ChemistryValidationError(
            f"Failed to parse ligand PDB file: {ligand_pdb}\nError: {e}"
        )

    if total_atoms == 0:
        raise ChemistryValidationError(
            f"No atoms found in ligand PDB file: {ligand_pdb}"
        )

    logger.debug(f"Atom counts: {element_counts}")
    logger.debug(f"Total atoms: {total_atoms}, Hydrogens: {hydrogen_count}")

    # Calculate total valence electrons
    total_valence_electrons = 0
    unsupported_elements = []

    for element, count in element_counts.items():
        if element in VALENCE_ELECTRONS:
            total_valence_electrons += VALENCE_ELECTRONS[element] * count
        else:
            unsupported_elements.append(element)

    if unsupported_elements:
        logger.warning(
            f"Unsupported elements (not counted in electron calculation): {unsupported_elements}"
        )

    logger.info(f"Total valence electrons: {total_valence_electrons}")

    # Check electron parity with charge
    if net_charge is None:
        # Try to infer charge - electrons should be even for closed-shell neutral molecules
        if total_valence_electrons % 2 == 0:
            inferred_charge = 0
            logger.info("Electron count is even - inferring neutral charge (0)")
        else:
            # Odd electrons - molecule is either charged or a radical
            # Common cases: +1 (cation, removes 1 electron) or -1 (anion, adds 1 electron)
            raise ChemistryValidationError(
                f"Ligand has odd number of valence electrons ({total_valence_electrons}).\n\n"
                "This indicates either:\n"
                "  1. Missing or extra atoms (most common)\n"
                "  2. Charged species (requires --charge flag)\n"
                "  3. Radical species (not supported)\n\n"
                f"Atom composition: {dict(element_counts)}\n\n"
                "Action required:\n"
                "  - If ligand is charged, specify --charge flag (e.g., --charge -1 for anion)\n"
                "  - If neutral, check for missing/extra hydrogens or atoms\n"
                "  - Use a structure preparation tool (e.g., PyMOL, Avogadro, OpenBabel)\n\n"
                f"Ligand file: {ligand_pdb}"
            )
        net_charge = inferred_charge
    else:
        # Validate provided charge
        electrons_after_charge = total_valence_electrons - net_charge
        
        if electrons_after_charge % 2 != 0:
            raise ChemistryValidationError(
                f"Invalid charge specification!\n\n"
                f"Total valence electrons: {total_valence_electrons}\n"
                f"Specified net charge: {net_charge}\n"
                f"Electrons after charge: {electrons_after_charge} (ODD - invalid!)\n\n"
                "For a closed-shell molecule, (electrons - charge) must be even.\n\n"
                f"Atom composition: {dict(element_counts)}\n\n"
                "Possible issues:\n"
                "  - Wrong charge value specified\n"
                "  - Missing or extra atoms in structure\n"
                "  - Incorrect protonation state\n\n"
                "Action required:\n"
                "  - Verify ligand structure in molecular viewer\n"
                "  - Check protonation state at physiological pH\n"
                "  - Use structure preparation tools to add/remove hydrogens\n\n"
                f"Ligand file: {ligand_pdb}"
            )

    # Warn about missing hydrogens (very common issue)
    heavy_atom_count = total_atoms - hydrogen_count
    
    if hydrogen_count == 0:
        raise ChemistryValidationError(
            f"No hydrogen atoms found in ligand structure!\n\n"
            f"Heavy atoms: {heavy_atom_count}\n"
            f"Hydrogens: {hydrogen_count}\n\n"
            "Antechamber requires explicit hydrogens for:\n"
            "  - Correct charge calculation (AM1-BCC)\n"
            "  - Proper atom typing (GAFF2)\n"
            "  - Bond order assignment\n\n"
            "Action required:\n"
            "  1. Add hydrogens using a structure preparation tool:\n"
            "     - PyMOL: h_add command\n"
            "     - Avogadro: Build → Add Hydrogens\n"
            "     - OpenBabel: obabel -h flag\n"
            "     - RDKit or other cheminformatics tools\n"
            "  2. Consider protonation state at target pH\n"
            "  3. Energy minimize after adding hydrogens\n\n"
            f"Ligand file: {ligand_pdb}"
        )

    # Heuristic check: typical H:heavy atom ratio is 1:1 to 3:1 for organic molecules
    if heavy_atom_count > 0:
        h_ratio = hydrogen_count / heavy_atom_count
        if h_ratio < 0.5:
            logger.warning(
                f"Low hydrogen:heavy-atom ratio ({h_ratio:.2f}). "
                f"Typical ratio is 1.0-2.0 for organic molecules. "
                "Structure may be missing hydrogens."
            )

    logger.success(
        f"Chemistry validation passed: {total_atoms} atoms, "
        f"{hydrogen_count} H, charge={net_charge}, "
        f"electrons={total_valence_electrons - net_charge} (even)"
    )

    return net_charge


def run_antechamber(
    input_file: Path,
    output_mol2: Path,
    charge_method: str = CHARGE_METHOD,
    force_field: str = LIGAND_FF,
    net_charge: Optional[int] = None,
) -> None:
    """
    Run antechamber to parameterize a ligand with GAFF2 and calculate AM1-BCC charges.

    Args:
        input_file: Path to ligand file (PDB, MOL2, or SDF)
        output_mol2: Path where output MOL2 file will be written
        charge_method: Charge calculation method (default: 'bcc' for AM1-BCC)
        force_field: Force field to use (default: 'gaff2')
        net_charge: Net charge of the ligand (if None, antechamber will guess)

    Raises:
        AntechamberError: If antechamber command fails
    """
    # Auto-detect input format
    input_suffix = input_file.suffix.lower()
    if input_suffix == ".pdb":
        input_format = "pdb"
    elif input_suffix == ".mol2":
        input_format = "mol2"
    elif input_suffix in [".sdf", ".sd"]:
        input_format = "sdf"
    else:
        input_format = "pdb"  # fallback
    
    logger.info(f"Running antechamber: {input_file} -> {output_mol2}")
    logger.debug(f"Input format: {input_format}, Charge method: {charge_method}, Force field: {force_field}")

    # Build antechamber command
    cmd = [
        "antechamber",
        "-i", str(input_file.absolute()),
        "-fi", input_format,
        "-o", str(output_mol2.absolute()),
        "-fo", "mol2",
        "-c", charge_method,
        "-at", force_field,
        "-pf", "y",  # Remove intermediate files
    ]

    if net_charge is not None:
        cmd.extend(["-nc", str(net_charge)])
        logger.debug(f"Net charge specified: {net_charge}")
    else:
        logger.warning(
            "Net charge not specified - antechamber will guess from atom types"
        )

    # Execute antechamber
    try:
        logger.debug(f"Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )

        # Log output
        if result.stdout:
            logger.debug(f"Antechamber stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"Antechamber stderr:\n{result.stderr}")

        # Verify output file was created
        if not output_mol2.exists():
            raise AntechamberError(
                f"Antechamber completed but output file not found: {output_mol2}"
            )

        logger.success(f"Antechamber completed: {output_mol2}")

    except subprocess.CalledProcessError as e:
        error_msg = (
            f"Antechamber failed with exit code {e.returncode}\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{e.stdout}\n\n"
            f"STDERR:\n{e.stderr}\n\n"
            "Inspect the ligand PDB file for structural issues.\n"
            "Common causes:\n"
            "  - Unrecognized atom types\n"
            "  - Invalid bond connectivity\n"
            "  - Missing hydrogens (antechamber requires all hydrogens)\n"
            "  - Incorrect geometry\n\n"
            f"Input file: {input_pdb}\n"
        )
        raise AntechamberError(error_msg)

    except FileNotFoundError:
        raise AntechamberError(
            "antechamber command not found. Ensure AmberTools is installed "
            "and the conda environment is activated."
        )


def run_parmchk2(
    mol2_file: Path,
    frcmod_file: Path,
    force_field: str = LIGAND_FF,
) -> None:
    """
    Run parmchk2 to generate missing force field parameters.

    Parmchk2 creates an frcmod file containing parameters not available
    in the standard GAFF2 force field.

    Args:
        mol2_file: Path to ligand MOL2 file from antechamber
        frcmod_file: Path where frcmod file will be written
        force_field: Force field to use (default: 'gaff2')

    Raises:
        ParmchkError: If parmchk2 command fails
    """
    logger.info(f"Running parmchk2: {mol2_file} -> {frcmod_file}")

    # Build parmchk2 command
    cmd = [
        "parmchk2",
        "-i", str(mol2_file.absolute()),
        "-f", "mol2",
        "-o", str(frcmod_file.absolute()),
        "-s", force_field,
    ]

    # Execute parmchk2
    try:
        logger.debug(f"Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )

        # Log output
        if result.stdout:
            logger.debug(f"Parmchk2 stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"Parmchk2 stderr:\n{result.stderr}")

        # Verify output file was created
        if not frcmod_file.exists():
            raise ParmchkError(
                f"Parmchk2 completed but output file not found: {frcmod_file}"
            )

        logger.success(f"Parmchk2 completed: {frcmod_file}")

    except subprocess.CalledProcessError as e:
        error_msg = (
            f"Parmchk2 failed with exit code {e.returncode}\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{e.stdout}\n\n"
            f"STDERR:\n{e.stderr}\n\n"
            "Inspect the MOL2 file for issues.\n"
            f"Input file: {mol2_file}\n"
        )
        raise ParmchkError(error_msg)

    except FileNotFoundError:
        raise ParmchkError(
            "parmchk2 command not found. Ensure AmberTools is installed "
            "and the conda environment is activated."
        )


def parameterize_ligand(
    ligand_mol2: Path,
    output_dir: Path,
    net_charge: Optional[int] = None,
    docked_pdb: Optional[Path] = None,
    ligand_resname: Optional[str] = None,
) -> tuple[Path, Path, Path]:
    """
    Complete ligand parameterization workflow for MM/GBSA.

    CRITICAL: This function requires a MOL2 or SDF file with complete bond information.
    PDB files alone are NOT sufficient because they lack bond orders required by antechamber.

    Workflow:
    1. If docked_pdb provided: transfer coordinates from PDB to MOL2 (preserves chemistry)
    2. Run antechamber for GAFF2 typing and AM1-BCC charges  
    3. Run parmchk2 to generate frcmod file

    Why MOL2/SDF is required:
    - Antechamber needs explicit bond orders for GAFF2 atom typing
    - PDB format does not encode bond orders (single/double/aromatic)
    - Inferring bonds from PDB geometry is unreliable
    - MOL2/SDF from ligand preparation tools (PyMOL, Avogadro, OpenBabel) have correct chemistry

    Args:
        ligand_mol2: Path to ligand MOL2 file with complete chemistry
        output_dir: Directory where output files will be written
        net_charge: Net charge of the ligand (default: 0)
        docked_pdb: Optional docked complex PDB for coordinate transfer
        ligand_resname: Required if docked_pdb provided - ligand residue name

    Returns:
        Tuple of (input_mol2_path, output_mol2_path, frcmod_path)

    Raises:
        LigandExtractionError: If coordinate transfer fails
        ChemistryValidationError: If chemistry validation fails
        AntechamberError: If antechamber fails
        ParmchkError: If parmchk2 fails
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ligand_mol2 = Path(ligand_mol2)
    if not ligand_mol2.exists():
        raise LigandExtractionError(f"Ligand MOL2 file not found: {ligand_mol2}")

    logger.info(f"Starting ligand parameterization")
    logger.info(f"  Input MOL2: {ligand_mol2}")
    logger.info(f"  Output directory: {output_dir}")

    # Define output file paths
    working_mol2 = output_dir / "ligand_input.mol2"
    output_mol2 = output_dir / "ligand.mol2"
    ligand_frcmod = output_dir / "ligand.frcmod"

    # Step 1: Coordinate transfer if docked PDB provided
    if docked_pdb:
        if not ligand_resname:
            raise LigandExtractionError(
                "ligand_resname is required when docked_pdb is provided"
            )
        
        logger.info(f"Transferring docked coordinates from PDB to MOL2")
        from rescore.parameterization.coordinates import transfer_coordinates
        
        transfer_coordinates(
            pdb_path=Path(docked_pdb),
            ligand_resname=ligand_resname,
            mol2_path=ligand_mol2,
            output_mol2=working_mol2,
        )
    else:
        # Use MOL2 as-is
        logger.info("Using MOL2 coordinates (no coordinate transfer)")
        import shutil
        shutil.copy(ligand_mol2, working_mol2)

    # Step 2: Set charge
    if net_charge is None:
        net_charge = 0
        logger.info(f"Using default net charge: {net_charge}")

    # Step 3: Run antechamber with specified charge
    run_antechamber(working_mol2, output_mol2, net_charge=net_charge)

    # Step 4: Run parmchk2
    run_parmchk2(output_mol2, ligand_frcmod)

    logger.success("Ligand parameterization complete")
    logger.info(f"  - Input MOL2:  {working_mol2}")
    logger.info(f"  - Output MOL2: {output_mol2}")
    logger.info(f"  - FRCMOD:      {ligand_frcmod}")

    return working_mol2, output_mol2, ligand_frcmod
