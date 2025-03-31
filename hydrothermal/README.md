# Hydrothermal Circulation Simulation

This project implements a 2D numerical simulation of hydrothermal circulation around a granite gneiss pluton using Python. The simulation models coupled heat transfer and fluid flow in a porous geological medium.

## Overview

The simulation solves coupled partial differential equations for heat transport and fluid flow in a heterogeneous geological domain. It includes an elliptical pluton (heat source) at the base of the model, with overlying host rock and sedimentary layers.

## Features

- 2D simulation with coupled heat transfer and fluid flow
- Temperature-dependent fluid properties (density, viscosity)
- Heterogeneous geological domain (pluton, host rock, sediments)
- Adaptive time stepping for efficient computation
- Comprehensive visualization of results
- Parametric study capabilities

## Requirements

- Python 3.x
- NumPy, SciPy, Matplotlib
- FiPy (finite volume PDE solver)
- Gmsh and meshio (mesh creation and handling)
- H5py (HDF5 file format support)
- PyVista (3D visualization)
- Optional: ParaView (advanced scientific visualization)

## Project Structure

```
hydrothermal/
├── src/
│   ├── main.py          # Main simulation driver
│   ├── domain.py        # Domain and mesh generation
│   ├── properties.py    # Physical properties
│   ├── equations.py     # Governing equations
│   ├── solver.py        # Numerical solver
│   └── visualization.py # Visualization tools
├── data/                # Input data and meshes
├── output/              # Simulation outputs
└── visualization/       # Generated visualizations
```

## Running the Simulation

To run the simulation with default parameters:

```bash
python src/main.py
```

### Command-line Options

- `--max-time`: Maximum simulation time in years (default: 100000)
- `--output-interval`: Output interval in years (default: 1000)
- `--mesh-size`: Mesh size near pluton in meters (default: 25)
- `--pluton-temp`: Initial pluton temperature in °C (default: 600)
- `--visualize-only`: Only run visualization on existing outputs
- `--show-gui`: Show Gmsh GUI during mesh generation
- `--parametric-study`: Run parametric study varying key properties

### Examples

Run a shorter simulation with higher resolution:
```bash
python src/main.py --max-time 10000 --mesh-size 15
```

Visualize existing results:
```bash
python src/main.py --visualize-only
```

Run a parametric study:
```bash
python src/main.py --parametric-study
```

## Physical Parameters

### Thermal Properties
- Thermal conductivity: pluton (3.0 W/m·K), host rock (2.5 W/m·K), sediments (1.8 W/m·K)
- Heat capacity: pluton (850 J/kg·K), host rock (800 J/kg·K), sediments (1000 J/kg·K)
- Density: pluton (2700 kg/m³), host rock (2500 kg/m³), sediments (2200 kg/m³)
- Initial temperature: 20°C at surface with 25°C/km geothermal gradient
- Pluton temperature: 600°C (cooling)

### Fluid Properties
- Density: temperature-dependent (ρ = 1000 - 0.375·(T-4)² kg/m³)
- Viscosity: temperature-dependent (μ = 2.414×10⁻⁵·10^(247.8/(T+133)) Pa·s)
- Heat capacity: 4200 J/kg·K

### Hydraulic Properties
- Permeability: pluton (10⁻¹⁸ m²), host rock (10⁻¹⁶ m²), sediments (10⁻¹⁴ m²)
- Porosity: pluton (0.05), host rock (0.1), sediments (0.25)

## Governing Equations

1. Heat transport equation with conduction and advection:
   ```
   ρᵣcᵣ(1-φ)∂T/∂t + ρₗcₗφvᵢ∂T/∂xᵢ = ∂/∂xᵢ(λ∂T/∂xᵢ)
   ```

2. Darcy's law for fluid flow:
   ```
   vᵢ = -(k/μ)(∂P/∂xᵢ + ρₗg)
   ```

3. Continuity equation:
   ```
   ∂/∂t(ρₗφ) + ∂/∂xᵢ(ρₗvᵢ) = 0
   ```

## Output and Visualization

The simulation generates HDF5 output files containing:
- Temperature distribution
- Velocity field
- Pressure field
- Geological unit information
- Material properties

Visualization tools provide:
- Temperature maps at different time steps
- Velocity vectors superimposed on temperature fields
- Heat flux distribution
- Temperature profiles along vertical and horizontal transects
- Time series of heat extraction rate and maximum velocity

## License

[MIT License](LICENSE)