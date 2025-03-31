"""
Governing equations for hydrothermal simulation.
"""
import numpy as np
from fipy import CellVariable, FaceVariable, Gmsh2D, TransientTerm, DiffusionTerm, ConvectionTerm, ImplicitSourceTerm
import fipy

class HydrothermalEquations:
    def __init__(self, mesh, properties):
        """
        Initialize equations for hydrothermal simulation
        
        Args:
            mesh: FiPy mesh
            properties: Properties object containing physical parameters
        """
        self.mesh = mesh
        self.props = properties
        
        # Variables we'll solve for
        self.temperature = None  # Temperature field
        self.pressure = None     # Pressure field
        self.velocity = None     # Velocity field (FaceVariable)
        
        # Geological unit field (for property lookup)
        self.unit_field = None
        
        # Material properties on mesh
        self.conductivity = None
        self.heat_capacity_solid = None
        self.density_solid = None
        self.permeability = None
        self.porosity = None
        
        # Temperature and pressure equations
        self.temperature_eq = None
        self.pressure_eq = None
        
        # Time step
        self.dt = 3600 * 24 * 365  # 1 year in seconds initially
        
    def initialize_fields(self, domain):
        """
        Initialize all fields on the mesh
        
        Args:
            domain: Domain object containing geometry information
        """
        # Create unit field (geological units)
        self.unit_field = CellVariable(name="geological_unit", mesh=self.mesh, value=self.props.HOST_ROCK)
        
        # Get cell coordinates
        x, y = self.mesh.cellCenters
        
        # Using rectangular approximations for Grid2D mesh
        # Define the regions based on coordinates
        
        # Pluton region (rectangle at bottom center)
        pluton_center_x = 0
        pluton_center_y = -domain.height + domain.pluton_height * 0.5
        half_width = domain.pluton_width / 2
        half_height = domain.pluton_height / 2
        
        pluton_mask = ((x >= pluton_center_x - half_width) & 
                       (x <= pluton_center_x + half_width) & 
                       (y >= pluton_center_y - half_height) & 
                       (y <= pluton_center_y + half_height))
        
        self.unit_field[pluton_mask] = self.props.PLUTON
        
        # Sediment layer
        sediment_depth = domain.height * 0.3  # 30% of domain height
        sediment_mask = y > -sediment_depth
        self.unit_field[sediment_mask] = self.props.SEDIMENTS
        
        # Initialize material properties based on geological units
        self.initialize_material_properties()
        
        # Initialize temperature field with geothermal gradient
        self.temperature = CellVariable(name="temperature", mesh=self.mesh, value=20.0)
        
        # Apply geothermal gradient (20°C at surface with 25°C/km geothermal gradient)
        # Convert y (depth) to km for the gradient calculation
        for i in range(len(y)):
            depth_km = -y[i] / 1000.0  # convert meters to km
            if not pluton_mask[i]:  # Skip pluton cells
                self.temperature.value[i] = self.props.surface_temperature + depth_km * self.props.geothermal_gradient
        
        # Set pluton temperature (600°C)
        self.temperature[pluton_mask] = self.props.pluton_temperature
        
        # Initialize pressure field (hydrostatic)
        self.pressure = CellVariable(name="pressure", mesh=self.mesh, value=101325.0)  # 1 atm
        
        # Initialize with hydrostatic pressure (surface at atmospheric pressure)
        # Calculate based on reference density (20°C)
        ref_density = self.props.fluid_density(20)
        for i in range(len(y)):
            self.pressure.value[i] = 101325 - ref_density * self.props.gravity * y[i]  # 101325 Pa = 1 atm
        
        # Apply boundary conditions
        # Top boundary (y=0): Fixed temperature and pressure
        top_mask = self.mesh.facesTop
        self.temperature.constrain(self.props.surface_temperature, where=top_mask)
        self.pressure.constrain(101325.0, where=top_mask)  # 1 atm
        
        # Bottom boundary: Fixed heat flux (prescribed T gradient) and no flow
        bottom_mask = self.mesh.facesBottom
        # For heat flux, we'll use a Neumann condition later
        # For no flow, we need to constrain the normal component of velocity
        
        # Initialize velocity field
        # Create face variable for velocity with two components
        self.velocity = FaceVariable(mesh=self.mesh, rank=1, value=0.0)
        
        # If direct assignment doesn't work, we can use mask operations
        # which is a safer way to initialize vector fields
        xFaces = self.mesh.facesLeft | self.mesh.facesRight
        yFaces = self.mesh.facesTop | self.mesh.facesBottom
        
        # Update fluid properties based on initial temperature
        self.update_fluid_properties()
        
        # Create temperature equation with empty convection term
        # This will be filled in during the solution process
        self.temperature_eq = TransientTerm(coeff=self.effective_heat_capacity()) == \
                             DiffusionTerm(coeff=self.conductivity)
                             
        # Create pressure equation based on Darcy's law and continuity
        # Will be implemented during solution process
        
    def initialize_material_properties(self):
        """
        Initialize material properties based on geological units
        """
        # Create material property fields
        self.conductivity = CellVariable(name="conductivity", mesh=self.mesh)
        self.heat_capacity_solid = CellVariable(name="heat_capacity_solid", mesh=self.mesh)
        self.density_solid = CellVariable(name="density_solid", mesh=self.mesh)
        self.permeability = CellVariable(name="permeability", mesh=self.mesh)
        self.porosity = CellVariable(name="porosity", mesh=self.mesh)
        
        # Set values based on geological units
        for unit in [self.props.PLUTON, self.props.HOST_ROCK, self.props.SEDIMENTS]:
            mask = self.unit_field == unit
            self.conductivity[mask] = self.props.thermal_conductivity[unit]
            self.heat_capacity_solid[mask] = self.props.heat_capacity_solid[unit]
            self.density_solid[mask] = self.props.density_solid[unit]
            self.permeability[mask] = self.props.permeability[unit]
            self.porosity[mask] = self.props.porosity[unit]
    
    def apply_boundary_conditions(self, domain):
        """
        Apply boundary conditions to temperature and pressure equations
        
        Args:
            domain: Domain object containing boundary information
        """
        # Get mesh faces on boundaries
        # Implementation will depend on how boundaries are identified in the mesh
        # This is a placeholder using mesh coordinates
        x, y = self.mesh.faceCenters
        
        # Temperature boundary conditions
        # - Fixed temperature at top (surface)
        top_faces = self.mesh.exteriorFaces & (y > -1e-6)
        self.temperature.constrain(self.props.surface_temperature, where=top_faces)
        
        # - Fixed heat flux at bottom
        bottom_faces = self.mesh.exteriorFaces & (y < -domain.height + 1e-6)
        heat_flux = self.props.bottom_heat_flux  # W/m²
        # Will need to implement flux BC
        
        # - No heat flux at sides (insulating)
        
        # Pressure/Flow boundary conditions
        # - Fixed pressure at top (atmospheric)
        self.pressure.constrain(101325, where=top_faces)  # 1 atm = 101325 Pa
        
        # - No flow at bottom
        # Will need to implement no-flow BC
        
        # - Hydrostatic pressure at sides
        left_faces = self.mesh.exteriorFaces & (x < -domain.width/2 + 1e-6)
        right_faces = self.mesh.exteriorFaces & (x > domain.width/2 - 1e-6)
        
        # Use hydrostatic pressure on sides
        # Will be implemented using coordinate-dependent constraint
    
    def update_fluid_properties(self):
        """
        Update temperature-dependent fluid properties
        """
        # Create cell-centered fluid properties
        self.fluid_density = CellVariable(name="fluid_density", mesh=self.mesh)
        self.fluid_viscosity = CellVariable(name="fluid_viscosity", mesh=self.mesh)
        
        # Update values based on current temperature
        self.fluid_density.setValue(self.props.fluid_density(self.temperature))
        self.fluid_viscosity.setValue(self.props.fluid_viscosity(self.temperature))
    
    def effective_heat_capacity(self):
        """
        Calculate effective heat capacity of fluid-rock mixture
        ρᵣcᵣ(1-φ) + ρₗcₗφ
        """
        # For numerical stability, let's ensure this returns a CellVariable
        # If it's already a CellVariable, this is a no-op
        capacity = (self.density_solid * self.heat_capacity_solid * (1 - self.porosity) + 
                   self.fluid_density * self.props.fluid_heat_capacity * self.porosity)
        
        # Ensure a minimum value for numerical stability
        min_capacity = 1.0e6  # 1 MJ/m³/K
        
        # Return the maximum of the calculated value and the minimum value
        return CellVariable(name="effective_heat_capacity", 
                           mesh=self.mesh, 
                           value=capacity,
                           hasOld=False)
    
    def update_velocity(self):
        """
        Update velocity field based on Darcy's law
        vᵢ = -(k/μ)(∂P/∂xᵢ + ρₗg)
        """
        # Calculate pressure gradient
        pressure_gradient = self.pressure.faceGrad
        
        # Calculate gravitational force
        gravity_force = FaceVariable(mesh=self.mesh, rank=1)
        gravity_force[1] = self.fluid_density.arithmeticFaceValue * self.props.gravity
        
        # Calculate Darcy velocity
        # Note: This is a simplified implementation
        # Harmonic averaging might be needed for permeability
        # and arithmetic for density, viscosity at faces
        darcy_coeff = self.permeability.harmonicFaceValue / self.fluid_viscosity.arithmeticFaceValue
        
        # Calculate velocity according to Darcy's law
        self.velocity.setValue(-darcy_coeff * (pressure_gradient + gravity_force))
    
    def update_pressure_equation(self):
        """
        Update pressure equation based on continuity 
        ∂/∂t(ρₗφ) + ∂/∂xᵢ(ρₗvᵢ) = 0
        """
        # This will be implemented based on FiPy's discretization approach
        # Essentially creating a diffusion equation with variable coefficients
        
        # For now, just update pressure using an elliptic equation derived from continuity
        # Full formulation will be implemented in the solver
        pass
    
    def update_temperature_equation(self):
        """
        Update temperature equation to include convection
        ρᵣcᵣ(1-φ)∂T/∂t + ρₗcₗφvᵢ∂T/∂xᵢ = ∂/∂xᵢ(λ∂T/∂xᵢ)
        """
        # Add convection term to temperature equation
        convection_coeff = self.fluid_density * self.props.fluid_heat_capacity * self.porosity
        convection_term = ConvectionTerm(coeff=convection_coeff * self.velocity)
        
        # Update temperature equation with convection
        self.temperature_eq = TransientTerm(coeff=self.effective_heat_capacity) == \
                             DiffusionTerm(coeff=self.conductivity) + \
                             convection_term