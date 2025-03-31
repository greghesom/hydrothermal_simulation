"""
Numerical solver for hydrothermal simulation.
"""
import numpy as np
import time
import os
import h5py
from fipy import Gmsh2D, Viewer
import fipy

class HydrothermalSolver:
    def __init__(self, domain, properties, equations, output_dir="output"):
        """
        Initialize the solver for hydrothermal simulation
        
        Args:
            domain: Domain object with mesh information
            properties: Properties object with physical parameters
            equations: HydrothermalEquations object with governing equations
            output_dir: Directory to save outputs
        """
        # Import configuration values
        from src.config import (
            MAX_SIMULATION_TIME, OUTPUT_INTERVAL,
            MIN_TIME_STEP, MAX_TIME_STEP, INITIAL_TIME_STEP,
            MAX_ITERATIONS, TEMPERATURE_TOLERANCE, PRESSURE_TOLERANCE, RELAXATION_FACTOR,
            DEBUG_MAX_STEPS
        )
        
        self.domain = domain
        self.props = properties
        self.eqns = equations
        
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Simulation control parameters
        self.t = 0.0  # Current simulation time (seconds)
        self.max_time = 3600 * 24 * 365 * MAX_SIMULATION_TIME  # Convert years to seconds
        self.min_dt = 3600 * 24 * MIN_TIME_STEP  # Convert days to seconds
        self.max_dt = 3600 * 24 * 365 * MAX_TIME_STEP  # Convert years to seconds
        self.dt = 3600 * 24 * 365 * INITIAL_TIME_STEP  # Convert years to seconds
        
        # Output control
        self.output_interval = 3600 * 24 * 365 * OUTPUT_INTERVAL  # Convert years to seconds
        self.next_output_time = 0.0  # Time for next output
        
        # Convergence parameters
        self.temperature_tolerance = TEMPERATURE_TOLERANCE
        self.pressure_tolerance = PRESSURE_TOLERANCE
        self.max_iterations = MAX_ITERATIONS
        self.relaxation_factor = RELAXATION_FACTOR
        
        # Debug parameters
        self.debug_max_steps = DEBUG_MAX_STEPS
        
        # Tracking variables
        self.temperature_history = []  # Store temperature evolution at key points
        self.heat_flux_history = []    # Store heat extraction rates
        self.max_velocity_history = [] # Store maximum velocity
    
    def initialize(self):
        """
        Initialize the simulation
        """
        # Initialize fields on the mesh
        self.eqns.initialize_fields(self.domain)
        
        # Apply boundary conditions
        self.eqns.apply_boundary_conditions(self.domain)
        
        # Update fluid properties based on initial temperature
        self.eqns.update_fluid_properties()
        
        # Initialize time step and simulation time
        self.t = 0.0  # Start at time zero
        self.dt = 3600 * 24 * 365  # 1 year in seconds initially
        self.next_output_time = 0.0  # First output at time zero
        
        # Save initial state
        self.save_output()
    
    def solve_step(self):
        """
        Solve one time step of the coupled equations
        
        Returns:
            converged (bool): Whether the solution converged
            max_error (float): Maximum error in the iteration
        """
        # Initialize error tracking
        max_error = float('inf')
        iteration = 0
        converged = False
        
        # Store previous temperature and pressure for convergence check
        temp_old = self.eqns.temperature.copy()
        pres_old = self.eqns.pressure.copy()
        
        # Iterative solution of coupled equations
        while max_error > max(self.temperature_tolerance, self.pressure_tolerance) and iteration < self.max_iterations:
            # 1. Update fluid properties based on current temperature
            self.eqns.update_fluid_properties()
            
            # 2. Solve pressure equation and update velocity field
            self.solve_pressure()
            self.eqns.update_velocity()
            
            # 3. Solve temperature equation with current velocity field
            self.solve_temperature()
            
            # 4. Calculate errors
            temp_error = max(abs(self.eqns.temperature - temp_old))
            pres_error = max(abs(self.eqns.pressure - pres_old))
            max_error = max(temp_error, pres_error / 100)  # Scale pressure error
            
            # 5. Update old values with relaxation
            temp_old = temp_old * (1 - self.relaxation_factor) + self.eqns.temperature * self.relaxation_factor
            pres_old = pres_old * (1 - self.relaxation_factor) + self.eqns.pressure * self.relaxation_factor
            
            iteration += 1
            
            print(f"Iteration {iteration}: max error = {max_error:.6f}")
        
        converged = (iteration < self.max_iterations)
        return converged, max_error
    
    def solve_pressure(self):
        """
        Solve the pressure equation
        """
        # Create a simple diffusion equation for pressure
        # In a real implementation, this would be derived from the continuity equation
        from fipy import TransientTerm, DiffusionTerm
        
        # Create a simple diffusion coefficient 
        diffusion_coeff = 1e-10  # A small value for stability
        
        # Simple pressure equation with diffusion term
        pressure_eq = DiffusionTerm(coeff=diffusion_coeff)
        
        # Solve the equation
        pressure_eq.solve(var=self.eqns.pressure)
        
        # The solver will automatically respect boundary conditions set on self.eqns.pressure
    
    def solve_temperature(self):
        """
        Solve the temperature equation
        """
        # In a real implementation, this would solve the full temperature equation
        # with conduction and advection terms
        from fipy import TransientTerm, DiffusionTerm
        
        # For now, implement a simple conduction equation
        # using the thermal conductivity from the equations object
        diffusion_coeff = self.eqns.conductivity / self.eqns.effective_heat_capacity()
        
        # Simple conduction equation
        temp_eq = TransientTerm() == DiffusionTerm(coeff=diffusion_coeff)
        
        # Solve the equation
        temp_eq.solve(var=self.eqns.temperature, dt=self.dt)
        
        # The solver will automatically respect boundary conditions set on self.eqns.temperature
    
    def adapt_time_step(self, converged, error):
        """
        Adapt time step based on solution convergence
        
        Args:
            converged: Whether the solution converged
            error: Maximum error in the iteration
        """
        if not converged:
            # Reduce time step if not converged
            self.dt = max(self.min_dt, self.dt * 0.5)
        elif error < self.temperature_tolerance * 0.1:
            # Increase time step if error is small
            self.dt = min(self.max_dt, self.dt * 1.5)
        
        # Ensure time step is at least the minimum value
        # This is critical to avoid getting stuck at dt = 0
        if self.dt < self.min_dt:
            self.dt = self.min_dt
            
        print(f"Adapted time step to {self.dt / (3600 * 24 * 365):.2f} years")
    
    def run_simulation(self):
        """
        Run the simulation to completion
        """
        print(f"Starting simulation, max time: {self.max_time / (3600 * 24 * 365):.0f} years")
        
        # Initialize simulation
        self.initialize()
        
        # Set an initial non-zero time step
        self.dt = 3600 * 24 * 365  # 1 year in seconds initially
        
        # Maximum number of steps (to prevent infinite loop)
        max_steps = 1000
        step_count = 0
        
        # Main time loop
        while self.t < self.max_time and step_count < max_steps:
            # Adjust time step to hit output times exactly
            if self.t + self.dt > self.next_output_time:
                self.dt = self.next_output_time - self.t
            
            # Ensure time step is at least the minimum value
            if self.dt < self.min_dt:
                self.dt = self.min_dt
            
            print(f"\nTime: {self.t / (3600 * 24 * 365):.2f} years, dt: {self.dt / (3600 * 24 * 365):.2f} years")
            
            # Solve one step
            converged, error = self.solve_step()
            
            if converged:
                # Update time
                self.t += self.dt
                
                # Adapt time step for next iteration
                self.adapt_time_step(converged, error)
                
                # Check if output is needed
                if self.t >= self.next_output_time:
                    self.save_output()
                    self.next_output_time += self.output_interval
                
                # Update tracking variables
                self.update_tracking()
            else:
                print("Solution did not converge, reducing time step and retrying")
                self.adapt_time_step(converged, error)
            
            # Increment step counter
            step_count += 1
            
            # For debugging, limit to the configured maximum number of steps
            if step_count >= self.debug_max_steps:
                print(f"Reached maximum step count of {step_count}, stopping for debugging")
                break
        
        print("Simulation complete")
        self.save_final_results()
    
    def update_tracking(self):
        """
        Update tracking variables for analysis
        """
        # Track temperature at key points
        # For example, tracking a point near the top of the pluton
        # This requires finding the closest cell to the point of interest
        
        # Track maximum fluid velocity
        max_vel = np.sqrt(self.eqns.velocity[0]**2 + self.eqns.velocity[1]**2).max()
        self.max_velocity_history.append((self.t, max_vel))
        
        # Track heat extraction rate from pluton
        # This requires calculating heat flux across pluton boundary
        # Will be implemented in the full solver
    
    def save_output(self):
        """
        Save the current state to file
        """
        output_file = os.path.join(self.output_dir, f"hydrothermal_{self.t / (3600 * 24 * 365):.0f}yrs.h5")
        
        with h5py.File(output_file, 'w') as f:
            # Save mesh coordinates
            x, y = self.eqns.mesh.cellCenters
            f.create_dataset('x', data=x)
            f.create_dataset('y', data=y)
            
            # Save temperature field
            f.create_dataset('temperature', data=self.eqns.temperature.value)
            
            # Save velocity field
            # Need to get velocity at cell centers for visualization
            f.create_dataset('velocity_x', data=self.eqns.velocity[0])
            f.create_dataset('velocity_y', data=self.eqns.velocity[1])
            
            # Save pressure field
            f.create_dataset('pressure', data=self.eqns.pressure.value)
            
            # Save material properties
            f.create_dataset('geological_unit', data=self.eqns.unit_field.value)
            f.create_dataset('permeability', data=self.eqns.permeability.value)
            
            # Save simulation metadata
            f.attrs['time'] = self.t
            f.attrs['time_years'] = self.t / (3600 * 24 * 365)
        
        print(f"Saved output to {output_file}")
    
    def save_final_results(self):
        """
        Save final results and analysis
        """
        # Save tracking data
        tracking_file = os.path.join(self.output_dir, "tracking_data.h5")
        with h5py.File(tracking_file, 'w') as f:
            # Convert tracking data to arrays
            times = np.array([t for t, _ in self.max_velocity_history])
            velocities = np.array([v for _, v in self.max_velocity_history])
            
            f.create_dataset('times', data=times)
            f.create_dataset('max_velocities', data=velocities)
            
            # Save temperature history if available
            if self.temperature_history:
                temps = np.array([temp for _, temp in self.temperature_history])
                f.create_dataset('temperature_history', data=temps)
            
            # Save heat flux history if available
            if self.heat_flux_history:
                heat_flux = np.array([flux for _, flux in self.heat_flux_history])
                f.create_dataset('heat_flux_history', data=heat_flux)
        
        print(f"Saved tracking data to {tracking_file}")
        
    def validate_against_analytical(self):
        """
        Validate results against analytical solutions where possible
        """
        # This would implement validation tests
        # For example, compare against conduction-only solution
        # or velocity in simple cases
        pass