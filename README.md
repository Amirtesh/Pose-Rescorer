# Physics-Based Post-Docking Rescoring Tool

**Physics-based post-docking rescoring tool (MM/GBSA under the hood).**

A streamlined Python package for **post-docking, single-frame rescoring** of protein-ligand complexes using **AmberTools MMPBSA.py** (MM/GBSA with optional PBSA).

## 🎯 Scientific Positioning

**This tool is designed for POST-DOCKING RESCORING and RELATIVE RANKING only:**
- ✅ Rank docked poses or compound series
- ✅ Prioritize compounds for synthesis/testing
- ✅ Compare binding modes qualitatively
- ❌ NOT for absolute binding free energies
- ❌ NOT for publication-quality ΔG values
- ❌ NOT a replacement for rigorous alchemical methods

The MM/GBSA scores provide **relative binding strength** for comparison within a series, not thermodynamically rigorous free energies.

## ✨ Features

### Core Capabilities
- **Single-frame MM/GBSA** rescoring using AmberTools MMPBSA.py
- **Restrained minimization** before scoring 
- **Integrated pipeline** from PDB → parameters → complex → scores
- **Batch processing** for multiple ligands against one receptor
- **Automated parameterization** with GAFF2 + AM1-BCC charges
- **Protein preparation** with pdb4amber + ff14SB force field
- **Structure validation** with chemistry checks
- **Optional plotting** for score visualization
- **Rapid Perturbation Sampling (RPS)** for uncertainty quantification

### Force Fields & Methods
- **Protein**: AMBER ff14SB
- **Ligand**: GAFF2 (General Amber Force Field 2)
- **Charges**: AM1-BCC (semi-empirical quantum mechanics)
- **Solvation**: Generalized Born (GB) implicit solvent (igb=5, OBC II)
- **Minimization**: Restrained (backbone only), 100 cycles, GB implicit solvent
- **Salt**: 0.15 M (physiological)

## 🚫 Limitations

- ❌ No absolute binding free energies (relative ranking only)
- ❌ No entropy calculations (single-frame, no conformational sampling)
- ❌ No explicit solvent (GB implicit solvent only)
- ❌ No metal cofactors support
- ❌ No covalent ligands
- ❌ No multi-ligand complexes
- ❌ Ligands must be complete (no missing heavy atoms)
- ❌ Requires pre-docked complex or separate ligand MOL2

## 📋 Requirements

- **AmberTools** 24.8+ (provides MMPBSA.py, antechamber, parmchk2, tleap, pdb4amber)
- **Python** 3.10+
- **BioPython** 1.86+
- **Typer** 0.21.1+ (CLI framework)
- **Rich** 13.7.1+ (colored output)
- **Loguru** 0.7.3+ (logging)
- **Matplotlib** 3.10+ (plotting)
- **Seaborn** 0.13+ (plotting)
- **Pandas** 2.3+ (RPS data handling)
- **NumPy** 2.2+ (numerical operations)

## 🔧 Installation

```bash
# Create conda environment with AmberTools
conda env create -f environment.yml
conda activate mmgbsa-dev

# Install package in development mode
pip install -e .

# Verify installation
rescore --version
```

## 📖 Usage

### ⚡ Minimization Feature

By default, **restrained energy minimization** is performed before MM/GBSA scoring to remove clashes from docking (similar to Schrödinger Prime MM/GBSA):

- **Protocol**: 100 cycles (50 steepest descent + 50 conjugate gradient)
- **Restraints**: Protein backbone atoms (N, CA, C) held fixed with 2.0 kcal/mol/Å² restraint
- **Free to move**: Ligand and protein side chains
- **Solvation**: GB implicit solvent (igb=5) with 12 Å cutoff
- **Output**: `minimized_complex.pdb` saved for visualization

**Control minimization:**
```bash
rescore integrate ... --minimize      # Default (ON)
rescore integrate ... --no-minimize   # Skip minimization
rescore batch ... --no-minimize       # Also works with batch
```

**Why minimize?** Docking produces poses with steric clashes. Minimization relaxes these while preserving the binding mode, improving MM/GBSA accuracy.

---

### Command Overview

```bash
rescore --help
```

Available commands:
- `validate` - Validate protein-ligand complex structure
- `prep-protein` - Prepare receptor with pdb4amber + tleap
- `parameterize` - Parameterize ligand with GAFF2 + AM1-BCC
- `assemble` - Assemble protein-ligand complex topology
- `run` - Run MM/GBSA rescoring on assembled complex
- `integrate` - End-to-end pipeline (all steps in one command)
- `batch` - Batch rescoring for multiple ligands vs one receptor

---

### 1️⃣ Validate Structure

Check if protein-ligand complex is suitable for rescoring:

```bash
rescore validate complex.pdb --ligand LIG
```

**Checks performed:**
- Structure parsing (valid PDB format)
- Ligand presence (at least one ligand molecule)
- Protein presence (at least one protein chain)
- Chemistry validation (ligand has bond information if from MOL2)

---

### 2️⃣ Prepare Protein

Process receptor PDB with pdb4amber and generate topology:

```bash
rescore prep-protein receptor.pdb -o protein_params/
```

**Output:**
- `protein_params/protein.prmtop` - AMBER topology file
- `protein_params/protein.inpcrd` - Coordinates
- `protein_params/receptor_processed.pdb` - Cleaned PDB

**Options:**
- `--skip-pdb4amber` - Skip pdb4amber cleaning (use if already prepared)

---

### 3️⃣ Parameterize Ligand

Generate GAFF2 parameters and AM1-BCC charges:

```bash
rescore parameterize ligand.mol2 -o ligand_params/
```

**Requirements:**
- Ligand must be in **MOL2 or SDF format** (with bond information)
- All heavy atoms must be present
- 3D coordinates required

**Output:**
- `ligand_params/ligand.mol2` - Parameterized ligand
- `ligand_params/ligand.frcmod` - GAFF2 force field modifications
- `ligand_params/ligand.pdb` - Converted structure

**Options:**
- `--net-charge INT` - Specify net charge (auto-detected if not provided)

---

### 4️⃣ Assemble Complex

Combine protein and ligand topologies:

```bash
rescore assemble \
  --protein protein_params/ \
  --ligand ligand_params/ \
  -o complex_params/
```

**Output:**
- `complex_params/complex.prmtop` - Combined topology
- `complex_params/complex.inpcrd` - Combined coordinates
- `complex_params/protein.pdb` - Protein structure (for tleap)
- `complex_params/tleap_complex.in` - Tleap input script

---

### 5️⃣ Run Rescoring

Execute single-frame MM/GBSA or MM/PBSA rescoring:

```bash
rescore run complex_params/ -o rescore_results/ --method pb  # use --method gb for GB (default)
```

**Output (method-dependent filenames):**
- GB: `rescore_results/rescore_output.dat`, `rescore_results/rescore.in`, `rescore_results/rescore.log`
- PB: `rescore_results/mmpbsa_output.dat`, `rescore_results/mmpbsa.in`, `rescore_results/mmpbsa.log`
- `rescore_results/ligand.prmtop` - Ligand topology (generated)

**Energy components:**
- **ΔG_bind** - Total binding energy (DELTA TOTAL)
- **ΔH** - Enthalpy contribution (DELTA G gas)
- **ΔG_GB** or **ΔG_PB** - Polar solvation energy (DELTA G solv)

---

### 6️⃣ Integrated Pipeline

Run all steps in one command:

```bash
rescore integrate \
  --receptor receptor.pdb \
  --ligand ligand.mol2 \
  -o pipeline_output/
```

**Workflow:**
1. Prepare protein (pdb4amber + tleap)
2. Parameterize ligand (GAFF2 + AM1-BCC)
3. Assemble complex (tleap combine)
4. Validate assembled complex (optional)
5. Restrained minimization (removes clashes, default ON)
6. Run MM/GBSA or MM/PBSA rescoring

**Output structure:**
```
pipeline_output/
├── protein/         # Receptor parameters
├── ligand/          # Ligand parameters
├── complex/         # Combined topology
├── rescore/         # Rescoring results (same folder name for GB/PB)
│   ├── minimization/        # Minimization outputs (if --minimize)
│   │   ├── minimized_complex.pdb      # Minimized structure for visualization
│   │   ├── minimized.inpcrd           # Minimized coordinates
│   │   ├── minimization_info.txt      # Detailed sander output
│   │   └── minimization.in            # Sander input file
│   ├── rescore_output.dat   # Energy results (GB)
│   ├── rescore.in           # MMPBSA.py input
│   ├── rescore.log          # Detailed calculation log
│   └── ligand_score_components.png    # Score plot (if --plot)
└── rps_analysis/    # RPS results (if --rps-replicates used)
    ├── replicate_001/
    ├── replicate_002/
    ├── ligand_rps_replicates.csv
    ├── ligand_rps_summary.csv
    └── ligand_rps_distribution.png
```

**Options:**
- `--skip-validation` - Skip structure validation step
- `--skip-pdb4amber` - Skip pdb4amber processing
- `--no-minimize` - Skip restrained minimization (default: ON)
- `--method gb|pb` - Solvation method (default: gb)
- `--plot` - Generate bar chart of score components
- `--rps-replicates N` - Enable Rapid Perturbation Sampling with N replicates
- `--rps-sigma FLOAT` - RPS perturbation magnitude in Angstroms (default: 0.2)

**Rapid Perturbation Sampling (RPS):**
Optional diagnostic feature for uncertainty quantification. Generates N perturbed ligand structures (Gaussian noise, σ = 0.2 Å default), rescores each, and reports statistics (mean, std, CV, 95% CI).

**Chemistry-Frozen Guarantee:**
- Ligand is parameterized **only once** in the main workflow (GAFF2 + AM1-BCC)
- Each RPS replicate **uses the same charges, atom types, and force field parameters**
- **antechamber is NOT run per replicate** - only coordinates are perturbed
- This ensures chemistry is frozen and only geometry varies

```bash
# With RPS uncertainty analysis
rescore integrate \
  --receptor receptor.pdb \
  --ligand ligand.mol2 \
  -o pipeline_output/ \
  --rps-replicates 20
```

**Critical:** RPS is NOT molecular dynamics, conformational sampling, or pose optimization. It is numerical sensitivity analysis for uncertainty quantification only.

---

### 7️⃣ Batch Processing

Rescore multiple ligands against one receptor:

```bash
rescore batch \
  --receptor receptor.pdb \
  --ligands ligands_directory/ \
  -o batch_results/
```

**Features:**
- Prepares receptor **once** (topology reused for all ligands)
- Processes ligands in parallel-safe workflow
- Validates receptor integrity throughout (hash checking)
- Outputs CSV with all results for easy analysis

**Output:**
```
batch_results/
├── receptor/                    # Prepared once
├── ligands/                     # Per-ligand subdirectories
│   ├── compound_1/
│   │   ├── ligand/
│   │   ├── complex/
│   │   └── rescore/             # Same directory name for any method
│   ├── compound_2/
│   └── ...
└── rescore_batch_results.csv    # Summary table (name does not change)
```

**CSV format (method-specific column):**
```csv
ligand_name,delta_g_bind,delta_h,delta_g_gb,status   # GB runs
ligand_name,delta_g_bind,delta_h,delta_g_pb,status   # PB runs
```

**Requirements:**
- All ligands must be in **MOL2 format**
- Ligands should be pre-docked or have coordinates
- Chemistry consistency enforced (same force field for all)

**Options:**
- `--method gb|pb` - Solvation method (default: gb)
- `--minimize / --no-minimize` - Minimization control (default: ON)
- `--plot` - Generate bar chart comparing ligand scores
- `--skip-pdb4amber` - Skip receptor preprocessing

**Note:** RPS (Rapid Perturbation Sampling) is NOT available in batch mode. Use `rescore integrate` for ligand-specific diagnostic analysis.

---

## 📊 Interpreting Results

### Energy Components

From `rescore_output.dat` (GB) or `mmpbsa_output.dat` (PB):

```
DELTA TOTAL    -6.79    # ΔG_bind: Total binding score
DELTA G gas   -49.75    # ΔH: Gas-phase enthalpy
DELTA G solv   42.96    # ΔG_GB or ΔG_PB: Solvation penalty
```

### Interpretation Guidelines

✅ **What you CAN do:**
- Rank compounds within a series (more negative = more favorable)
- Compare docking poses for the same ligand
- Identify top candidates for experimental validation
- Detect outliers or problematic structures

❌ **What you CANNOT do:**
- Report absolute ΔG values in publications
- Compare results across different receptors
- Predict experimental Ki/Kd values directly
- Use for QSAR model training without validation

### Typical Values

- **Strong binders:** -10 to -15 kcal/mol
- **Moderate binders:** -5 to -10 kcal/mol
- **Weak binders:** 0 to -5 kcal/mol
- **Non-binders:** > 0 kcal/mol

**Note:** These are RELATIVE scores, not experimental affinities!

---

## 🔬 Scientific Background

### MM/GBSA Method

MM/GBSA (Molecular Mechanics - Generalized Born Surface Area) estimates binding free energies by:

```
ΔG_bind = ΔH - TΔS
        ≈ ΔE_MM + ΔG_solv - TΔS
```

Where:
- **ΔE_MM** = Gas-phase molecular mechanics energy (AMBER ff14SB + GAFF2)
- **ΔG_solv** = Solvation free energy (GB implicit solvent)
- **TΔS** = Entropy (NOT calculated in single-frame mode)

### Why Single-Frame?

This tool uses **single-frame rescoring** (one structure per complex):
- ✅ Fast (seconds per complex)
- ✅ Suitable for post-docking ranking
- ✅ Consistent with docking pose selection
- ❌ No conformational sampling
- ❌ No entropy calculation
- ❌ High sensitivity to input structure quality

### Comparison to Alternatives

| Method | Speed | Accuracy | Use Case |
|--------|-------|----------|----------|
| **MM/GBSA (this tool)** | Fast | Low-Medium | Post-docking ranking |
| **MM/PBSA** | Medium | Medium | Better solvation model |
| **FEP/TI** | Slow | High | Rigorous ΔG calculations |
| **Docking scores** | Very Fast | Low | Initial screening |

---

## ⚠️ Important Warnings

### Critical Reminders

1. **NOT thermodynamically rigorous** - These are approximations for ranking, not absolute free energies

2. **Single-frame limitations** - No conformational sampling means results depend heavily on input structure quality

3. **No entropy** - Binding entropy is NOT calculated, introducing systematic error

4. **Implicit solvent** - GB model approximates solvation, missing explicit water effects

5. **Force field dependence** - Results depend on AMBER ff14SB + GAFF2 parameterization quality

### Common Pitfalls

- ❌ Using different receptors/conditions between compounds
- ❌ Comparing results to experimental ΔG directly
- ❌ Over-interpreting small differences (< 1 kcal/mol)
- ❌ Ignoring structure validation failures
- ❌ Using PDB ligands without bond information

### Best Practices

- ✅ Always validate structures before rescoring
- ✅ Use consistent docking protocol for all ligands
- ✅ Check for clashes or unusual geometries
- ✅ Compare within series, not across receptors
- ✅ Treat scores as relative rankings, not predictions
- ✅ Validate top hits experimentally

---

## 🐛 Troubleshooting

### Common Issues

**Issue:** `ambpdb failed to convert protein coordinates`
- **Cause:** Spaces in directory paths
- **Solution:** Avoid spaces in file/directory names, or use relative paths

**Issue:** `No ligand molecules detected in structure`
- **Cause:** Validation expects ligand in complex PDB
- **Solution:** Use `--skip-validation` for integrate/batch commands, or provide complex PDB

**Issue:** `Only MOL2 ligands are accepted for batch processing`
- **Cause:** Batch mode requires consistent chemistry (GAFF2)
- **Solution:** Convert ligands to MOL2 format with bond information

**Issue:** `MMPBSA.py requires separate receptor and ligand topologies`
- **Cause:** Directory structure doesn't match expected layout
- **Solution:** Follow output structure conventions (protein/, ligand/, complex/ subdirs)

**Issue:** `ligand.mol2 missing or invalid`
- **Cause:** Ligand parameterization failed or file corrupted
- **Solution:** Check antechamber logs, ensure ligand has complete chemistry

---

## 📁 Output File Descriptions

### Topology Files (.prmtop)
AMBER topology files containing atom types, charges, bonds, angles, dihedrals

### Coordinate Files (.inpcrd)
AMBER coordinate files with 3D positions (Å)

### Force Field Modifications (.frcmod)
GAFF2 parameter additions for ligand-specific atom types

### MM/GBSA Results (.dat)
Energy component table with ΔG_bind, ΔH, ΔG_GB values

### Tleap Scripts (.in)
Input scripts for tleap topology generation (for reproducibility)

### Logs (.log)
Detailed execution logs from AmberTools programs

---

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=rescore --cov-report=html

# Run specific test
pytest tests/test_validation.py
```

---

## 📚 References

### AmberTools
- **AMBER Force Fields:** Case et al., J. Comput. Chem. (2005)
- **GAFF:** Wang et al., J. Comput. Chem. (2004)
- **AM1-BCC:** Jakalian et al., J. Comput. Chem. (2002)
- **MM/PBSA:** Kollman et al., Acc. Chem. Res. (2000)

### Recommended Reading
- Genheden & Ryde, *Expert Opin. Drug Discov.* (2015) - MM/PBSA review
- Wang et al., *Front. Mol. Biosci.* (2019) - End-point free energy methods
- Miller et al., *J. Chem. Theory Comput.* (2012) - MMPBSA.py tool

---

## 📄 License

This project uses AmberTools, which is distributed under GNU GPL v3.

---

## Acknowledgments

- **AmberTools** developers for MMPBSA.py and force field infrastructure
- **BioPython** team for PDB parsing tools
- **Typer/Rich** for modern CLI development

---

## 📮 Contact

For bug reports, feature requests, or questions about usage, please open an issue on GitHub.

**Remember:** This tool is for **relative ranking only**. Always validate top hits experimentally!
