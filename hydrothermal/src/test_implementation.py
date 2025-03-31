"""
Test the implementation of the hydrothermal simulation.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import simulation modules
from src.domain import Domain
from src.properties import Properties
from src.equations import HydrothermalEquations

def test_properties():
    """
    Test the properties module.
    """
    print("Testing Properties module:")
    props = Properties()
    
    # Test fluid property functions
    temps = [20, 50, 100, 200, 300]
    
    print("Temperature (°C) | Density (kg/m³) | Viscosity (Pa·s)")
    print("-" * 60)
    
    for temp in temps:
        density = props.fluid_density(temp)
        viscosity = props.fluid_viscosity(temp)
        print(f"{temp:14.1f} | {density:14.2f} | {viscosity:.4e}")
    
    print("\nInitial temperature profile:")
    depths = [0, 1000, 2000, 3000, 4000, 5000]
    domain_height = 5000
    
    print("Depth (m) | Temperature (°C)")
    print("-" * 30)
    
    for depth in depths:
        y = -depth  # Convert depth to y-coordinate (negative downward)
        temp = props.initial_temperature(y, domain_height)
        print(f"{depth:8.0f} | {temp:16.2f}")
    
    print("\nProperties test completed")
    return True

def test_domain():
    """
    Test the domain module.
    """
    print("\nTesting Domain module:")
    
    # Create a small test domain
    test_domain = Domain(
        width=5000,        # 5 km
        height=2500,       # 2.5 km
        pluton_width=1500, # 1.5 km
        pluton_height=750, # 0.75 km
        mesh_size_pluton=100,
        mesh_size_far=400,
        mesh_file="data/test_mesh.msh"
    )
    
    # Create mesh
    os.makedirs("data", exist_ok=True)
    mesh_file = test_domain.create_mesh(show_gui=False)
    
    print(f"Created test mesh at {mesh_file}")
    
    # Calculate and display mesh stats (using meshio)
    import meshio
    mesh = meshio.read(mesh_file)
    
    num_points = len(mesh.points)
    num_cells = sum(len(cells.data) for cells in mesh.cells)
    
    print(f"Mesh statistics:")
    print(f"  Number of points: {num_points}")
    print(f"  Number of cells: {num_cells}")
    
    # Visualize mesh
    plt.figure(figsize=(10, 6))
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            for triangle in cell_block.data:
                x = [mesh.points[i][0] for i in triangle]
                y = [mesh.points[i][1] for i in triangle]
                # Close the triangle
                x.append(x[0])
                y.append(y[0])
                plt.plot(x, y, 'k-', linewidth=0.2)
    
    plt.axis('equal')
    plt.title('Test Mesh')
    plt.savefig('data/test_mesh.png')
    plt.close()
    
    print("  Mesh visualization saved to data/test_mesh.png")
    print("Domain test completed")
    return True

def main():
    """
    Run the implementation tests.
    """
    print("Running implementation tests...")
    
    # Create output directories
    os.makedirs("data", exist_ok=True)
    os.makedirs("output", exist_ok=True)
    
    # Run tests
    tests = [
        test_properties,
        test_domain
    ]
    
    results = [test() for test in tests]
    
    if all(results):
        print("\nAll implementation tests passed!")
    else:
        print("\nSome implementation tests failed.")

if __name__ == "__main__":
    main()