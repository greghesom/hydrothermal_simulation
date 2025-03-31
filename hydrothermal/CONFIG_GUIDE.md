# Configuration Guide for Hydrothermal Simulation

This document provides guidance on how to configure the hydrothermal simulation by modifying the central configuration file.

## Overview

All main simulation parameters are centralized in the `src/config.py` file. This makes it easy to adjust the simulation behavior without modifying the core code.

## Key Configuration Sections

### Domain Configuration

These parameters control the physical dimensions of the simulation domain:

```python
# DOMAIN CONFIGURATION
DOMAIN_WIDTH = 10000         # Width of domain in meters (x-direction)
DOMAIN_HEIGHT = 5000         # Height of domain in meters (y-direction)
PLUTON_WIDTH = 3000          # Width of pluton in meters
PLUTON_HEIGHT = 1500         # Height of pluton in meters
```

### Mesh Parameters

These parameters control the mesh resolution:

```python
# Mesh parameters
MESH_SIZE_PLUTON = 25        # Mesh size near pluton in meters
MESH_SIZE_FAR = 200          # Mesh size far from pluton in meters
```

Decreasing these values will create a finer mesh with more detail but will increase computation time.

### Material Properties

These parameters control the physical properties of the geological units:

```python
# Thermal properties
THERMAL_CONDUCTIVITY = {
    'PLUTON': 3.0,           # W/m·K
    'HOST_ROCK': 2.5,        # W/m·K
    'SEDIMENTS': 1.8         # W/m·K
}

# Hydraulic properties
PERMEABILITY = {
    'PLUTON': 1e-18,         # m²
    'HOST_ROCK': 1e-16,      # m²
    'SEDIMENTS': 1e-14       # m²
}
```

### Boundary Conditions

These parameters control the temperature and pressure conditions at the boundaries:

```python
# BOUNDARY CONDITIONS
SURFACE_TEMPERATURE = 20     # °C
GEOTHERMAL_GRADIENT = 25     # °C/km
PLUTON_TEMPERATURE = 600     # °C
BOTTOM_HEAT_FLUX = 0.06      # W/m² (60 mW/m²)
```

### Simulation Control

These parameters control how the simulation runs:

```python
# SIMULATION CONTROL
MAX_SIMULATION_TIME = 100000   # Maximum simulation time in years
OUTPUT_INTERVAL = 1000         # Output interval in years
MIN_TIME_STEP = 1              # Minimum time step in days
MAX_TIME_STEP = 1000           # Maximum time step in years
INITIAL_TIME_STEP = 1          # Initial time step in years
```

### Solver Parameters

These parameters control the numerical solver behavior:

```python
# Solver parameters
MAX_ITERATIONS = 100           # Maximum iterations per time step
TEMPERATURE_TOLERANCE = 1e-3   # Convergence tolerance for temperature (°C)
PRESSURE_TOLERANCE = 10.0      # Convergence tolerance for pressure (Pa)
RELAXATION_FACTOR = 0.7        # Relaxation factor for iterations
```

### Debug Settings

For testing and debugging:

```python
# Debug settings
DEBUG_MAX_STEPS = 10           # Maximum number of steps when debugging
```

Increase this value to run more time steps during testing.

### Fracture Zones

Configure fracture zones in the model:

```python
# Default fracture zones
DEFAULT_FRACTURE_ZONES = [
    {"x": -PLUTON_WIDTH/4, "width": 200, "permeability_factor": 10},
    {"x": PLUTON_WIDTH/4, "width": 200, "permeability_factor": 10}
]
```

### Visualization

Control the appearance of plots:

```python
# VISUALIZATION
TEMPERATURE_COLORMAP = 'viridis'
VELOCITY_COLORMAP = 'Blues'
QUIVER_SCALE = 0.02
QUIVER_DENSITY = 20
```

## Overriding Configuration via Command Line

You can override some configuration parameters via command line arguments:

```bash
python src/main.py --max-time 10000 --pluton-temp 650 --mesh-size 20
```

This will run the simulation with a maximum time of 10,000 years, a pluton temperature of 650°C, and a mesh size of 20m near the pluton.

## Examples

### Running a Short Test Simulation

Edit `config.py` and set:
```python
DEBUG_MAX_STEPS = 5     # Just 5 time steps
MAX_SIMULATION_TIME = 10  # Only simulate 10 years
```

### Increasing Resolution for Final Run

Edit `config.py` and set:
```python
MESH_SIZE_PLUTON = 15   # Finer mesh near pluton
MESH_SIZE_FAR = 150     # Finer mesh everywhere
DEBUG_MAX_STEPS = 1000  # Allow many time steps
```

### Testing Different Pluton Properties

Edit `config.py` and set:
```python
PLUTON_TEMPERATURE = 700  # Hotter pluton
THERMAL_CONDUCTIVITY['PLUTON'] = 3.5  # Higher thermal conductivity
```