"""
Physical properties for hydrothermal simulation.
"""
import numpy as np

class Properties:
    def __init__(self):
        # Import configuration values
        from src.config import (
            THERMAL_CONDUCTIVITY, HEAT_CAPACITY_SOLID, DENSITY_SOLID,
            PERMEABILITY, POROSITY, FLUID_HEAT_CAPACITY,
            SURFACE_TEMPERATURE, GEOTHERMAL_GRADIENT, PLUTON_TEMPERATURE, BOTTOM_HEAT_FLUX
        )
        
        # Region identifiers
        self.PLUTON = 1
        self.HOST_ROCK = 2
        self.SEDIMENTS = 3
        
        # Thermal properties
        self.thermal_conductivity = {
            self.PLUTON: THERMAL_CONDUCTIVITY['PLUTON'],
            self.HOST_ROCK: THERMAL_CONDUCTIVITY['HOST_ROCK'],
            self.SEDIMENTS: THERMAL_CONDUCTIVITY['SEDIMENTS']
        }
        
        self.heat_capacity_solid = {
            self.PLUTON: HEAT_CAPACITY_SOLID['PLUTON'],
            self.HOST_ROCK: HEAT_CAPACITY_SOLID['HOST_ROCK'],
            self.SEDIMENTS: HEAT_CAPACITY_SOLID['SEDIMENTS']
        }
        
        self.density_solid = {
            self.PLUTON: DENSITY_SOLID['PLUTON'],
            self.HOST_ROCK: DENSITY_SOLID['HOST_ROCK'],
            self.SEDIMENTS: DENSITY_SOLID['SEDIMENTS']
        }
        
        # Hydraulic properties
        self.permeability = {
            self.PLUTON: PERMEABILITY['PLUTON'],
            self.HOST_ROCK: PERMEABILITY['HOST_ROCK'],
            self.SEDIMENTS: PERMEABILITY['SEDIMENTS']
        }
        
        self.porosity = {
            self.PLUTON: POROSITY['PLUTON'],
            self.HOST_ROCK: POROSITY['HOST_ROCK'],
            self.SEDIMENTS: POROSITY['SEDIMENTS']
        }
        
        # Constant fluid properties
        self.fluid_heat_capacity = FLUID_HEAT_CAPACITY
        
        # Gravity
        self.gravity = 9.81  # m/s²
        
        # Initial/boundary conditions
        self.surface_temperature = SURFACE_TEMPERATURE
        self.geothermal_gradient = GEOTHERMAL_GRADIENT
        self.pluton_temperature = PLUTON_TEMPERATURE
        self.bottom_heat_flux = BOTTOM_HEAT_FLUX
        
    def interpolate_property_to_mesh(self, mesh, property_name):
        """
        Interpolate a property onto the mesh based on geological units
        """
        # This will be implemented using FiPy CellVariable
        # to store properties on the mesh
        pass
    
    def fluid_density(self, temperature):
        """
        Temperature-dependent fluid density
        ρ = 1000 - 0.375·(T-4)² kg/m³
        
        Args:
            temperature: Temperature in °C
            
        Returns:
            Fluid density in kg/m³
        """
        return 1000 - 0.375 * (temperature - 4)**2
    
    def fluid_viscosity(self, temperature):
        """
        Temperature-dependent fluid viscosity
        μ = 2.414×10⁻⁵·10^(247.8/(T+133)) Pa·s
        
        Args:
            temperature: Temperature in °C
            
        Returns:
            Fluid viscosity in Pa·s
        """
        return 2.414e-5 * 10**(247.8 / (temperature + 133))
    
    def initial_temperature(self, y, domain_height):
        """
        Initial temperature distribution based on depth and geothermal gradient
        
        Args:
            y: Vertical coordinate (m), negative downward
            domain_height: Total domain height (m)
            
        Returns:
            Temperature in °C
        """
        # Convert depth to km for geothermal gradient calculation
        depth_km = -y / 1000.0  # y is negative, so negative y is depth
        
        # Apply geothermal gradient
        return self.surface_temperature + depth_km * self.geothermal_gradient
    
    def permeability_with_fractures(self, base_permeability, fracture_factor=1.0):
        """
        Adjust permeability for fracture zones
        
        Args:
            base_permeability: Base permeability from geological unit
            fracture_factor: Factor by which to increase permeability (1.0 = no fracture)
            
        Returns:
            Adjusted permeability
        """
        return base_permeability * fracture_factor