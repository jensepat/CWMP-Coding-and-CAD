
# CWMP TRIGA Reactor OpenMC Models

## Overview

This directory contains OpenMC computational models for the Oregon State TRIGA Reactor (OSTR) with three different target configurations. All models are developed in Python using the OpenMC nuclear simulation package.

## Repository Structure

The models are organized into three main directories:

### **Model 1: Slug Model** (`Model_1/`)
- **Target Material**: Solid Neptunium Dioxide (NpO₂) slug
- **Configuration**: Cylindrical target geometry
- **Key Files**:
  - `Slug_Model.py` - Main geometry and reactor model setup
  - `Slug_Materials.py` - Material definitions including NpO₂ target

### **Model 2: Pellet Model** (`Model_2/`)
- **Target Material**: HFIR-style Neptunium Cermet pellets (NpO₂-Al composite)
- **Configuration**: Stacked pellet arrangement in target locations
- **Key Files**:
  - `Pellet_Model.py` - Geometry setup for pellet stack configuration
  - `Pellet_Materials.py` - Material definitions including NpO₂-Al cermet composition

### **Model 3: Aqueous Model** (`Model_3/`)
- **Target Material**: Aqueous Neptunium solution (1.25 M Np-237 in water)
- **Configuration**: Liquid target in the lazy susan target location
- **Key Files**:
  - `Aqueous_Model.py` - Geometry setup for aqueous solution target
  - `Aqueous_Materials.py` - Material definitions including aqueous Np solution

## Dependencies

All models require the following Python packages:

```
numpy
os
argparse
openmc (v0.15.0)
openmc.deplete
```

## Running the Models

Each model can be run independently. The command-line interface is consistent across all models:

### Command-Line Arguments

```
-n              = (int)         Number of particle histories to track per cycle
-i              = (int)         Number of inactive cycles (for settling)
-a              = (int)         Number of active cycles (statistics tracked)
-m              = (bool)        Full power (1 MW) temperatures if included
--run           = (bool)        Execute OpenMC after generating input files
--plot          = (bool)        Generate geometry plots after completion
--rods          = (float) x4    Control rod heights: <TR> <SA> <SH> <REG> (0-100%)
--num_threads   = (int)         Number of OpenMP threads for parallel execution
--dir           = (str)         Output directory for input files and results
```

### Example Command

```bash
python Model_1/Slug_Model.py -n 25000 -i 30 -a 470 -m --run --plot --rods 30 40.1 50.2 60.3 --num_threads 16 --dir ./output_slug
```

This command runs Model 1 with:
- **Transient Rod**: 30% withdrawn
- **Safety Rod**: 40.1% withdrawn
- **Shim Rod**: 50.2% withdrawn
- **Regulating Rod**: 60.3% withdrawn
- **Temperatures**: 1 MW (full power) conditions
- **Simulation**: 25,000 particles/cycle, 30 inactive + 470 active cycles
- **Execution**: 16 parallel threads
- **Output**: Files saved to `./output_slug/`

### Default Configuration

```bash
python Model_1/Slug_Model.py -n 25000 -i 30 -a 470 -m --run --plot --rods 70 70 70 70 --num_threads 16 --dir ./output_slug
```

This is the standard baseline configuration with all control rods at 70% withdrawn.

## Model Outputs

Upon execution, each model generates:
1. **OpenMC XML input files** - Geometry, materials, and settings definitions
2. **Plot files** - Visual representations of the reactor geometry (if `--plot` enabled)
3. **Simulation results** - Keff and other neutronics metrics (if `--run` enabled)

## File Organization

```
CWMP Code/
├── README.md
├── Model_1/
│   ├── Slug_Model.py
│   └── Slug_Materials.py
├── Model_2/
│   ├── Pellet_Model.py
│   └── Pellet_Materials.py
└── Model_3/
    ├── Aqueous_Model.py
    └── Aqueous_Materials.py
```

## Key Specifications

### Reactor Configuration (All Models)
- **Reactor Type**: TRIGA (Training, Research, Isotope Production, General Atomics)
- **Core Composition**: TRIGA fuel elements with control rods
- **Coolant**: Light water at adjustable temperatures
- **Operating Power**: Up to 1 MW

### Material Library
All models include comprehensive material definitions:
- TRIGA fuel (U-235/U-238 in ZrH)
- Structural materials (aluminum, stainless steel, titanium)
- Moderators (graphite, boron carbide)
- Thermal scattering data (S(α,β) tables for hydrogen and thermal scattering)

## Notes

- Ensure both the model file and its corresponding materials file are in the same directory or Python path
- The OpenMC chain file path may need to be updated based on your system configuration
- Temperature settings affect cross-section interpolation and physics accuracy

## Citation

Courtesy of Paul Sprague for initial code framework and reactor model design.
