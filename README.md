# Hydrothermal Circulation Simulation

Numerical simulation of hydrothermal circulation around a cooling pluton in a geological setting, implemented in Python using FiPy.

## Overview

This project implements a 2D numerical simulation of hydrothermal circulation around a granite gneiss pluton. The simulation models coupled heat transfer and fluid flow in a porous geological medium, considering temperature-dependent fluid properties and heterogeneous geological units.

## Features

- 2D simulation of coupled heat transfer and fluid flow
- Temperature-dependent fluid density and viscosity
- Heterogeneous geological domain with pluton, host rock, and sedimentary layers
- Adaptive time stepping for numerical stability and efficiency
- Comprehensive visualization of simulation results

## Requirements

- Python 3.x
- NumPy, SciPy, Matplotlib
- FiPy (finite volume PDE solver)
- Gmsh and meshio (mesh creation and handling)
- H5py (HDF5 file format support)
- PyVista (for visualization)

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/hydrothermal_simulation.git
   cd hydrothermal_simulation
   ```

2. Create a virtual environment and activate it:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install required packages:
   ```
   pip install -r requirements.txt
   ```

## Usage

Run the simulation with default parameters:
```
python hydrothermal/src/main.py
```

With custom parameters:
```
python hydrothermal/src/main.py --max-time 10000 --mesh-size 20 --pluton-temp 650
```

Run visualization only on existing outputs:
```
python hydrothermal/src/main.py --visualize-only
```

Check all available options:
```
python hydrothermal/src/main.py --help
```

## Configuration

All simulation parameters can be adjusted in `hydrothermal/src/config.py`. See the [Configuration Guide](hydrothermal/CONFIG_GUIDE.md) for details.

## Project Structure

```
hydrothermal/
├── src/                  # Source code
│   ├── main.py           # Main entry point
│   ├── config.py         # Centralized configuration
│   ├── domain.py         # Domain and mesh generation
│   ├── properties.py     # Physical properties
│   ├── equations.py      # Governing equations
│   ├── solver.py         # Numerical solver
│   └── visualization.py  # Visualization tools
├── data/                 # Input data and meshes
├── output/               # Simulation outputs
└── visualization/        # Generated visualizations
```

## Physical Background

The simulation models a granite gneiss pluton cooling in the Earth's crust, driving hydrothermal circulation through the surrounding rock. This process is important in the formation of ore deposits, geothermal systems, and understanding crustal heat transfer.

Key physical processes included:
- Heat conduction and advection
- Buoyancy-driven fluid flow
- Permeability contrasts between geological units
- Temperature-dependent fluid properties

## License

[MIT License](LICENSE)

## Acknowledgments

This project uses the FiPy finite volume PDE solver developed at NIST.

## Contact

Your Name - youremail@example.com