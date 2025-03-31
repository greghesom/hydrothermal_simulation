"""
Initialize the hydrothermal simulation.
"""
import os
import sys
from src.domain import Domain

def main():
    """Create initial mesh for the simulation."""
    print("Initializing hydrothermal simulation...")
    
    # Create data directory
    os.makedirs("data", exist_ok=True)
    
    # Create domain and mesh
    domain = Domain(
        width=10000,         # 10 km
        height=5000,         # 5 km
        pluton_width=3000,   # 3 km
        pluton_height=1500,  # 1.5 km
        mesh_size_pluton=25,
        mesh_size_far=200,
        mesh_file="data/mesh.msh"
    )
    
    # Create mesh
    mesh_file = domain.create_mesh(show_gui=False)
    print(f"Created mesh at {mesh_file}")
    print("Initialization complete")

if __name__ == "__main__":
    main()