# FMAPIq

Calculate small-angle X-ray scattering intensity profiles I(q) from protein configurations in atomistic simulations.

## Quick Start

### Installation
```bash
# Using uv (recommended)
uv sync

# Or using pip
pip install ../fmapb2-pkg ../maeppi-pkg ../sim2iq-pkg
```

### Usage

**GUI Mode (recommended):**
```bash
uv run main.py
```

**Command Line:**
```bash
# Run complete pipeline
uv run main.py protein.pqr

# Load existing work directory
uv run main.py --load /path/to/workdir

# Run with custom parameters
uv run main.py protein.pqr --n-proteins 100 --n-reps 5
```

## What It Does

FMAPIq runs a 3-stage pipeline:

1. **Energy Maps** - Calculates protein interaction potentials from PQR structure
2. **Simulation** - Runs molecular dynamics to sample protein configurations
3. **Scattering** - Computes I(q) intensity profiles from configurations

## Requirements

- Python ≥ 3.8
- Dependencies listed in `pyproject.toml`
- For GUI: `gradio`
- For plotting: `matplotlib`

## Output

- `Iqs.txt` - Raw I(q) data for all configurations
- `Iqs.stat.txt` - Statistical summary (mean ± standard error)
- `iq_plot.png` - Visualization of results
- Packed archive for sharing results

## Options

Common parameters:
- `--n-proteins` - Number of proteins in simulation (default: 200)
- `--n-reps` - Number of replicate simulations (default: 3)
- `--ion` - Ionic strength in M (default: 0.016)
- `--temp` - Temperature in °C (default: 25)
- `--threads` - CPU threads to use

See `uv run main.py --help` for all options.

## Building Executable

```bash
chmod +x build.sh
./build.sh
```