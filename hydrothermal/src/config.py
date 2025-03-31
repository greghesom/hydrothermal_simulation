"""
Configuration file for hydrothermal simulation.

This file contains all adjustable parameters for the simulation.
Edit this file to change the behavior of the simulation.
"""

# DOMAIN CONFIGURATION
# -------------------
# Physical dimensions
DOMAIN_WIDTH = 10000         # Width of domain in meters (x-direction)
DOMAIN_HEIGHT = 5000         # Height of domain in meters (y-direction)
PLUTON_WIDTH = 3000          # Width of pluton in meters
PLUTON_HEIGHT = 1500         # Height of pluton in meters

# Mesh parameters
MESH_SIZE_PLUTON = 25        # Mesh size near pluton in meters
MESH_SIZE_FAR = 200          # Mesh size far from pluton in meters

# Simulation regions
SEDIMENT_LAYER_DEPTH = 0.3   # Fraction of domain height (0.3 = 30% of domain height)

# MATERIAL PROPERTIES
# ------------------
# Thermal properties
THERMAL_CONDUCTIVITY = {
    'PLUTON': 3.0,           # W/m·K
    'HOST_ROCK': 2.5,        # W/m·K
    'SEDIMENTS': 1.8         # W/m·K
}

HEAT_CAPACITY_SOLID = {
    'PLUTON': 850,           # J/kg·K
    'HOST_ROCK': 800,        # J/kg·K
    'SEDIMENTS': 1000        # J/kg·K
}

DENSITY_SOLID = {
    'PLUTON': 2700,          # kg/m³
    'HOST_ROCK': 2500,       # kg/m³
    'SEDIMENTS': 2200        # kg/m³
}

FLUID_HEAT_CAPACITY = 4200   # J/kg·K

# Hydraulic properties
PERMEABILITY = {
    'PLUTON': 1e-18,         # m²
    'HOST_ROCK': 1e-16,      # m²
    'SEDIMENTS': 1e-14       # m²
}

POROSITY = {
    'PLUTON': 0.05,
    'HOST_ROCK': 0.1,
    'SEDIMENTS': 0.25
}

# BOUNDARY CONDITIONS
# ------------------
SURFACE_TEMPERATURE = 20     # °C
GEOTHERMAL_GRADIENT = 25     # °C/km
PLUTON_TEMPERATURE = 600     # °C
BOTTOM_HEAT_FLUX = 0.06      # W/m² (60 mW/m²)

# SIMULATION CONTROL
# -----------------
# Time stepping
MAX_SIMULATION_TIME = 100000   # Maximum simulation time in years
OUTPUT_INTERVAL = 1000         # Output interval in years
MIN_TIME_STEP = 1              # Minimum time step in days
MAX_TIME_STEP = 1000           # Maximum time step in years
INITIAL_TIME_STEP = 1          # Initial time step in years

# Solver parameters
MAX_ITERATIONS = 100           # Maximum iterations per time step
TEMPERATURE_TOLERANCE = 1e-3   # Convergence tolerance for temperature (°C)
PRESSURE_TOLERANCE = 10.0      # Convergence tolerance for pressure (Pa)
RELAXATION_FACTOR = 0.7        # Relaxation factor for iterations

# Debug settings
DEBUG_MAX_STEPS = 10           # Maximum number of steps when debugging

# FRACTURE ZONES
# -------------
# Default fracture zones
DEFAULT_FRACTURE_ZONES = [
    {"x": -PLUTON_WIDTH/4, "width": 200, "permeability_factor": 10},
    {"x": PLUTON_WIDTH/4, "width": 200, "permeability_factor": 10}
]

# VISUALIZATION
# ------------
# Plotting parameters
TEMPERATURE_COLORMAP = 'viridis'
VELOCITY_COLORMAP = 'Blues'
QUIVER_SCALE = 0.02
QUIVER_DENSITY = 20