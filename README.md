# MM/GBSA Rescoring Tool

A minimal Python package for MM/GBSA rescoring of protein-ligand complexes using **AmberTools only**.

## Features

- **Protein force field**: ff14SB
- **Ligand force field**: GAFF2
- **Charges**: AM1-BCC
- **MM/GBSA**: Single-structure or small ensemble, no entropy

## Limitations

This tool is designed for **relative ranking only**:
- ❌ No absolute binding free energies
- ❌ No support for metals
- ❌ No support for covalent ligands
- ❌ No support for multiple ligands
- ❌ No support for missing heavy atoms

## Installation

```bash
# Activate the mmgbsa-dev environment
conda activate mmgbsa-dev

# Install in development mode
pip install -e .
```

## Usage

### Validate a structure

```bash
mmgbsa validate complex.pdb --ligand LIG
```

### Run MM/GBSA rescoring (not yet implemented)

```bash
mmgbsa rescore complex.pdb --ligand LIG --output results/
```

## Development

```bash
# Run tests
pytest tests/

# Format code
black src/ tests/
isort src/ tests/

# Type checking
mypy src/
```

## Requirements

- Python 3.10+
- AmberTools 24+ (must be available in PATH)
- Conda environment: `mmgbsa-dev`

## Project Structure

```
mmgbsa/
├── pyproject.toml          # Modern Python packaging config
├── README.md
├── src/mmgbsa/
│   ├── __init__.py         # Package initialization
│   ├── cli.py              # Typer-based CLI
│   ├── config.py           # Configuration constants
│   └── validation/
│       ├── __init__.py
│       ├── pdb_checks.py   # PDB validation logic
│       └── errors.py       # Custom exceptions
└── tests/
    └── test_validation.py  # Validation tests
```

## License

MIT
