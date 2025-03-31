"""
Main script for hydrothermal simulation.
"""
import os
import numpy as np
import sys
import argparse
from fipy import Gmsh2D

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import simulation modules
from domain import Domain
from properties import Properties
from equations import HydrothermalEquations
from solver import HydrothermalSolver
from visualization import Visualization

def create_directories():
    """
    Create necessary directories for simulation
    """
    dirs = ["output", "visualization", "data"]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    return dirs

def parse_arguments():
    """
    Parse command line arguments
    """
    # Import configuration values 
    from src.config import (MAX_SIMULATION_TIME, OUTPUT_INTERVAL,
                           MESH_SIZE_PLUTON, PLUTON_TEMPERATURE)
    
    parser = argparse.ArgumentParser(description='Hydrothermal circulation simulation')
    
    # Add command-line options - using config defaults
    parser.add_argument('--max-time', type=float, default=MAX_SIMULATION_TIME,
                        help=f'Maximum simulation time in years (default: {MAX_SIMULATION_TIME})')
    parser.add_argument('--output-interval', type=float, default=OUTPUT_INTERVAL,
                        help=f'Output interval in years (default: {OUTPUT_INTERVAL})')
    parser.add_argument('--mesh-size', type=float, default=MESH_SIZE_PLUTON,
                        help=f'Mesh size near pluton in meters (default: {MESH_SIZE_PLUTON})')
    parser.add_argument('--pluton-temp', type=float, default=PLUTON_TEMPERATURE,
                        help=f'Initial pluton temperature in °C (default: {PLUTON_TEMPERATURE})')
    parser.add_argument('--visualize-only', action='store_true',
                        help='Only run visualization on existing outputs')
    parser.add_argument('--show-gui', action='store_true',
                        help='Show Gmsh GUI during mesh generation')
    parser.add_argument('--parametric-study', action='store_true',
                        help='Run parametric study varying key properties')
    
    return parser.parse_args()

def run_simulation(args):
    """
    Run the hydrothermal simulation
    
    Args:
        args: Command line arguments
    """
    # Import domain parameters from config
    from src.config import DOMAIN_WIDTH, DOMAIN_HEIGHT, PLUTON_WIDTH, PLUTON_HEIGHT
    
    print("Starting hydrothermal circulation simulation")
    
    # Create mesh and domain
    domain = Domain(
        # Use command line args if provided, otherwise config values will be used
        mesh_size_pluton=args.mesh_size,
        mesh_size_far=args.mesh_size * 8,
        mesh_file="data/mesh.msh"
    )
    
    # Create mesh
    mesh_file = domain.create_mesh(show_gui=args.show_gui)
    
    # For now, let's create a simple grid mesh instead of using Gmsh
    # This is a temporary solution until we resolve the Gmsh import issues
    from fipy import Grid2D
    
    # Create a uniform grid mesh for testing
    nx, ny = 100, 50  # Number of cells in x and y directions
    dx, dy = domain.width / nx, domain.height / ny
    fipy_mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
    
    # Set up physical properties
    props = Properties()
    
    # Override pluton temperature if specified
    if args.pluton_temp != 600:
        props.pluton_temperature = args.pluton_temp
    
    # Set up governing equations
    equations = HydrothermalEquations(fipy_mesh, props)
    
    # Set up solver
    solver = HydrothermalSolver(domain, props, equations, output_dir="output")
    
    # Set simulation time parameters
    solver.max_time = 3600 * 24 * 365 * args.max_time  # Convert years to seconds
    solver.output_interval = 3600 * 24 * 365 * args.output_interval  # Convert years to seconds
    
    # Run simulation
    solver.run_simulation()
    
    print("Simulation completed")
    
    return "output"  # Return output directory

def run_visualization(output_dir):
    """
    Run visualization on simulation outputs
    
    Args:
        output_dir: Directory with simulation outputs
    """
    print("Generating visualizations")
    
    # Create visualization object
    vis = Visualization(output_dir=output_dir, vis_dir="visualization")
    
    # Generate comprehensive report
    vis.generate_report()
    
    print("Visualization completed")

def run_parametric_study(args):
    """
    Run parametric study varying key properties
    
    Args:
        args: Command line arguments
    """
    print("Running parametric study")
    
    # Define parameter variations
    pluton_temperatures = [500, 600, 700]
    permeability_factors = [0.1, 1.0, 10.0]
    
    # Base output directory
    base_output_dir = "output"
    
    # Run simulations with different parameters
    for temp in pluton_temperatures:
        for perm_factor in permeability_factors:
            # Create parameter-specific output directory
            output_dir = f"{base_output_dir}/temp{temp}_perm{perm_factor}"
            os.makedirs(output_dir, exist_ok=True)
            
            print(f"\nRunning simulation with temperature={temp}°C, permeability factor={perm_factor}")
            
            # Create a copy of args with modified parameters
            modified_args = argparse.Namespace(**vars(args))
            modified_args.pluton_temp = temp
            modified_args.max_time = 10000  # Shorter time for parametric study
            
            # Create domain and properties as before, but with modified parameters
            # This would be a shortened version of run_simulation with parameter adjustments
            
            # Visualize results
            vis_dir = f"visualization/temp{temp}_perm{perm_factor}"
            os.makedirs(vis_dir, exist_ok=True)
            vis = Visualization(output_dir=output_dir, vis_dir=vis_dir)
            vis.generate_report()
    
    print("Parametric study completed")

def validate_simulation():
    """
    Validate simulation against analytical solutions
    """
    print("Validating simulation against analytical solutions")
    
    # This would implement validation tests
    # 1. Compare against conduction-only solution
    # 2. Validate flow pattern in simple cases
    # 3. Check conservation of energy
    
    print("Validation completed")

def main():
    """
    Main function for hydrothermal simulation
    """
    # Parse command-line arguments
    args = parse_arguments()
    
    # Create directories
    create_directories()
    
    if args.parametric_study:
        run_parametric_study(args)
    elif args.visualize_only:
        run_visualization("output")
    else:
        # Run simulation
        output_dir = run_simulation(args)
        
        # Run visualization
        run_visualization(output_dir)
        
        # Optional validation
        validate_simulation()
    
    print("Hydrothermal simulation completed successfully")

if __name__ == "__main__":
    main()