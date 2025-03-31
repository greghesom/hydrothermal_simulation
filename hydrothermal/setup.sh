#!/bin/bash
# Setup script for hydrothermal simulation environment

echo "Setting up the hydrothermal simulation environment..."

# Create directory structure
mkdir -p data output visualization

# Check if the virtual environment is active
if [[ -z "${VIRTUAL_ENV}" ]]; then
    echo "Please activate the virtual environment first using:"
    echo "source ../activate_env.sh"
    exit 1
fi

# Check Python version
python_version=$(python --version)
echo "Using $python_version"

# Check required packages
echo "Checking required packages..."
pip list | grep -E 'numpy|scipy|matplotlib|fipy|gmsh|meshio|pyvista|h5py'

# Create a test mesh to verify Gmsh installation
echo "Testing Gmsh installation..."
cat > test_gmsh.py << 'EOF'
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("test")
gmsh.model.geo.addPoint(0, 0, 0, 0.1, 1)
gmsh.model.geo.addPoint(1, 0, 0, 0.1, 2)
gmsh.model.geo.addPoint(1, 1, 0, 0.1, 3)
gmsh.model.geo.addPoint(0, 1, 0, 0.1, 4)
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("data/test.msh")
gmsh.finalize()
print("Gmsh test completed successfully")
EOF

python test_gmsh.py
rm test_gmsh.py

# Create initialization script
echo "Creating initialization script..."
cat > initialize.py << 'EOF'
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
EOF

# Create a simple test to verify the environment works
echo "Creating test script..."
cat > test.py << 'EOF'
"""
Test script for hydrothermal simulation.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def test_numpy():
    """Test NumPy functionality."""
    print("Testing NumPy...")
    a = np.array([1, 2, 3])
    b = np.array([4, 5, 6])
    c = a + b
    print(f"  {a} + {b} = {c}")
    return True

def test_matplotlib():
    """Test Matplotlib functionality."""
    print("Testing Matplotlib...")
    plt.figure(figsize=(6, 4))
    x = np.linspace(0, 2*np.pi, 100)
    y = np.sin(x)
    plt.plot(x, y)
    plt.title('Sine Wave')
    plt.savefig('test_plot.png')
    plt.close()
    print("  Plot saved to test_plot.png")
    return True

def test_imports():
    """Test importing key modules."""
    print("Testing imports...")
    try:
        import fipy
        print("  FiPy imported successfully")
        import gmsh
        print("  Gmsh imported successfully")
        import meshio
        print("  meshio imported successfully")
        import h5py
        print("  h5py imported successfully")
        try:
            import pyvista
            print("  PyVista imported successfully")
        except ImportError:
            print("  PyVista not available (optional)")
        return True
    except ImportError as e:
        print(f"  Import error: {e}")
        return False

def main():
    """Run tests to verify environment."""
    print("Running environment tests...")
    
    tests = [
        test_numpy,
        test_matplotlib,
        test_imports
    ]
    
    results = [test() for test in tests]
    
    if all(results):
        print("\nAll tests passed! Environment is ready.")
    else:
        print("\nSome tests failed. Please check the environment setup.")
    
if __name__ == "__main__":
    main()
EOF

echo "Setup complete!"
echo "To initialize the simulation, run:"
echo "python initialize.py"
echo "To test the environment, run:"
echo "python test.py"
echo "To run the simulation, run:"
echo "python src/main.py"