# Hydrothermal Simulation Usage Guide

This document provides detailed instructions on how to use the hydrothermal simulation software.

## Setup and Installation

1. Ensure you have the required Python packages installed:
   ```
   source ../activate_env.sh  # Activate the environment
   ```

2. Run the setup script:
   ```
   ./setup.sh
   ```
   This will create the directory structure and test that all dependencies are properly installed.

## Getting Started

1. Initialize the simulation:
   ```
   python initialize.py
   ```
   This creates the initial mesh for the simulation domain.

2. Run a test to verify the implementation:
   ```
   python src/test_implementation.py
   ```
   This will test the properties and domain modules to ensure they are working correctly.

## Running Simulations

### Basic Simulation

To run a simulation with default parameters:

```
python src/main.py
```

This will:
- Create a 10km × 5km domain with an elliptical pluton
- Run the simulation for 100,000 years
- Save outputs every 1,000 years
- Generate visualizations in the "visualization" directory

### Command-line Options

The simulation can be customized with several command-line options:

```
python src/main.py [options]
```

Options:
- `--max-time VALUE`: Maximum simulation time in years (default: 100000)
- `--output-interval VALUE`: Output interval in years (default: 1000)
- `--mesh-size VALUE`: Mesh size near pluton in meters (default: 25)
- `--pluton-temp VALUE`: Initial pluton temperature in °C (default: 600)
- `--visualize-only`: Only run visualization on existing outputs
- `--show-gui`: Show Gmsh GUI during mesh generation
- `--parametric-study`: Run parametric study varying key properties

### Examples

1. Run a shorter simulation with higher resolution:
   ```
   python src/main.py --max-time 10000 --mesh-size 15
   ```

2. Run a simulation with a hotter pluton:
   ```
   python src/main.py --pluton-temp 700
   ```

3. Generate visualizations from previous simulation results:
   ```
   python src/main.py --visualize-only
   ```

4. Run a parametric study to explore different property values:
   ```
   python src/main.py --parametric-study
   ```

## Working with Results

### Output Files

Simulation results are saved in the "output" directory as HDF5 files:
- Each file corresponds to a specific time step
- File naming format: `hydrothermal_{time}yrs.h5`
- Tracking data is saved in `tracking_data.h5`

### Visualization

The simulation automatically generates visualizations in the "visualization" directory:
- Temperature distributions
- Velocity fields
- Geological unit maps
- Permeability maps
- Temperature profiles
- Time series of key parameters

### Custom Visualization

You can create custom visualizations by importing the Visualization class:

```python
from src.visualization import Visualization

# Create visualization object
vis = Visualization(output_dir="output", vis_dir="my_custom_vis")

# List output files
files = vis.list_output_files()

# Load a specific output file
data = vis.load_output_file(files[0])

# Create temperature plot
vis.plot_temperature(data, output_file="my_temperature_plot.png")

# Create temperature with velocity plot
vis.plot_temperature_with_velocity(data, output_file="my_temp_vel_plot.png")
```

## Advanced Usage

### Customizing Physical Properties

You can modify physical properties by editing the `properties.py` file:

- Thermal properties (conductivity, heat capacity, density)
- Hydraulic properties (permeability, porosity)
- Fluid properties (density and viscosity functions)
- Initial and boundary conditions

### Customizing the Domain

To customize the simulation domain, modify the `domain.py` file:

- Domain dimensions
- Pluton size and location
- Mesh resolution
- Geological unit boundaries

### Adding Fracture Zones

To add fracture zones with enhanced permeability, modify the `setup_fracture_zones` method in the `domain.py` file.

### Using ParaView for Advanced Visualization

For more advanced 3D visualization, you can load the HDF5 output files in ParaView:

1. Open ParaView
2. Choose File > Open
3. Navigate to the output directory and select an HDF5 file
4. Apply the "Table To Points" filter to create a 3D representation
5. Add filters and color mappings as needed

## Troubleshooting

### Common Issues

1. **"Module not found" errors**:
   Make sure you've activated the virtual environment:
   ```
   source ../activate_env.sh
   ```

2. **Mesh generation errors**:
   Check that Gmsh is properly installed:
   ```
   python -c "import gmsh; gmsh.initialize(); gmsh.finalize()"
   ```

3. **Slow simulation**:
   - Increase mesh size (e.g., `--mesh-size 50`)
   - Reduce maximum simulation time (e.g., `--max-time 10000`)

4. **Convergence problems**:
   Modify the relaxation factor and convergence tolerance in `solver.py`

### Getting Help

If you encounter problems, check:
1. The README.md file for general information
2. The code comments for specific implementation details
3. The test_implementation.py file for examples of how to use the modules