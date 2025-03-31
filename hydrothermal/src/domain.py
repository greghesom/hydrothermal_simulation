"""
Domain setup and mesh generation for hydrothermal simulation.
"""
import numpy as np
import gmsh
import meshio
import os

class Domain:
    def __init__(self, 
                 width=None,          # Width of domain in meters (x-direction)
                 height=None,         # Height of domain in meters (y-direction)
                 pluton_width=None,   # Width of pluton in meters
                 pluton_height=None,  # Height of pluton in meters
                 mesh_size_pluton=None, # Mesh size near pluton in meters
                 mesh_size_far=None,  # Mesh size far from pluton in meters
                 mesh_file="mesh.msh"):
        
        # Import configuration values
        from src.config import (DOMAIN_WIDTH, DOMAIN_HEIGHT, 
                               PLUTON_WIDTH, PLUTON_HEIGHT,
                               MESH_SIZE_PLUTON, MESH_SIZE_FAR,
                               SEDIMENT_LAYER_DEPTH)
        
        # Use config values as defaults
        self.width = width if width is not None else DOMAIN_WIDTH
        self.height = height if height is not None else DOMAIN_HEIGHT
        self.pluton_width = pluton_width if pluton_width is not None else PLUTON_WIDTH
        self.pluton_height = pluton_height if pluton_height is not None else PLUTON_HEIGHT
        self.mesh_size_pluton = mesh_size_pluton if mesh_size_pluton is not None else MESH_SIZE_PLUTON
        self.mesh_size_far = mesh_size_far if mesh_size_far is not None else MESH_SIZE_FAR
        self.mesh_file = mesh_file
        
        # Store sediment depth for use in mesh generation
        self.sediment_depth = self.height * SEDIMENT_LAYER_DEPTH
        
        # Domain regions (for assigning properties)
        self.PLUTON = 1
        self.HOST_ROCK = 2
        self.SEDIMENTS = 3
        
        # Boundary tags
        self.TOP_BOUNDARY = 10
        self.BOTTOM_BOUNDARY = 11
        self.LEFT_BOUNDARY = 12
        self.RIGHT_BOUNDARY = 13
        
    def create_mesh(self, show_gui=False):
        """
        Create a 2D mesh for the hydrothermal simulation using gmsh.
        - Rectangular domain with pluton
        - Refinement near the pluton
        """
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        
        if show_gui:
            gmsh.fltk.initialize()
        else:
            gmsh.option.setNumber("General.Verbosity", 2)
        
        # Create model
        model_name = "HydrothermalModel"
        gmsh.model.add(model_name)
        
        # Define domain dimensions (in meters)
        x_min, x_max = -self.width/2, self.width/2
        y_min, y_max = -self.height, 0  # Top at y=0, bottom at y=-height
        
        # Create rectangle for the entire domain
        domain_tag = gmsh.model.occ.addRectangle(x_min, y_min, 0, self.width, self.height)
        
        # Create sediment boundary position
        sediment_depth = self.height * 0.3  # 30% of domain height is sediment layer
        
        # Create rectangle for sediment layer
        sediment_tag = gmsh.model.occ.addRectangle(x_min, -sediment_depth, 0, self.width, sediment_depth)
        
        # Create pluton at the bottom center
        pluton_center_x = 0
        pluton_center_y = -self.height + self.pluton_height * 0.5
        
        # Create rectangle for pluton
        pluton_tag = gmsh.model.occ.addRectangle(
            pluton_center_x - self.pluton_width/2, 
            pluton_center_y - self.pluton_height/2, 
            0, 
            self.pluton_width, 
            self.pluton_height
        )
        
        # Fragment the domain to create separate regions
        # First, cut the domain with the sediment (to separate sediment and host rock)
        out_fragment, _ = gmsh.model.occ.fragment([(2, domain_tag)], [(2, sediment_tag)])
        
        # Now cut the result with the pluton
        out_fragment, _ = gmsh.model.occ.fragment(out_fragment, [(2, pluton_tag)])
        
        # Synchronize
        gmsh.model.occ.synchronize()
        
        # Extract all surfaces after fragmentation
        all_surfaces = gmsh.model.getEntities(2)
        
        # Identify different regions by their position
        
        # Find pluton (centered around pluton_center_x, pluton_center_y)
        pluton_surface = None
        
        # Find sediment layer (top portion)
        sediment_surface = None
        
        # Find host rock (middle portion)
        host_rock_surface = None
        
        # Use coordinates to identify regions
        for surface in all_surfaces:
            surface_id = surface[1]
            com = gmsh.model.occ.getCenterOfMass(2, surface_id)
            
            # Check if this is the pluton (by position)
            if (com[0] >= pluton_center_x - self.pluton_width/2 and 
                com[0] <= pluton_center_x + self.pluton_width/2 and
                com[1] >= pluton_center_y - self.pluton_height/2 and
                com[1] <= pluton_center_y + self.pluton_height/2):
                pluton_surface = surface_id
            
            # Check if this is the sediment layer (by y-coordinate)
            elif com[1] > -sediment_depth:
                sediment_surface = surface_id
            
            # Otherwise, it's the host rock
            else:
                host_rock_surface = surface_id
        
        # Add physical groups for domain regions
        gmsh.model.addPhysicalGroup(2, [pluton_surface], self.PLUTON)
        gmsh.model.addPhysicalGroup(2, [sediment_surface], self.SEDIMENTS)
        gmsh.model.addPhysicalGroup(2, [host_rock_surface], self.HOST_ROCK)
        
        # Get boundary curves to set boundary conditions
        boundaries = gmsh.model.getBoundary(all_surfaces, recursive=False)
        
        # Extract curves
        top_curves = []
        bottom_curves = []
        left_curves = []
        right_curves = []
        
        for boundary in boundaries:
            curve_id = boundary[1]
            boundary_points = gmsh.model.getBoundary([(1, curve_id)], recursive=False)
            
            # If the curve has two points
            if len(boundary_points) == 2:
                # Get coordinates of the two points
                p1 = gmsh.model.getValue(0, abs(boundary_points[0][1]), [])
                p2 = gmsh.model.getValue(0, abs(boundary_points[1][1]), [])
                
                # Horizontal lines (top and bottom)
                if abs(p1[1] - p2[1]) < 1e-6:
                    if abs(p1[1] - y_max) < 1e-6:
                        top_curves.append(curve_id)
                    elif abs(p1[1] - y_min) < 1e-6:
                        bottom_curves.append(curve_id)
                
                # Vertical lines (left and right)
                if abs(p1[0] - p2[0]) < 1e-6:
                    if abs(p1[0] - x_min) < 1e-6:
                        left_curves.append(curve_id)
                    elif abs(p1[0] - x_max) < 1e-6:
                        right_curves.append(curve_id)
        
        # Add physical groups for boundaries
        gmsh.model.addPhysicalGroup(1, top_curves, self.TOP_BOUNDARY)
        gmsh.model.addPhysicalGroup(1, bottom_curves, self.BOTTOM_BOUNDARY)
        gmsh.model.addPhysicalGroup(1, left_curves, self.LEFT_BOUNDARY)
        gmsh.model.addPhysicalGroup(1, right_curves, self.RIGHT_BOUNDARY)
        
        # Set mesh size
        # Add a distance field to refine mesh near pluton
        field_id = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_id, "PointsList", [pluton_surface])
        
        threshold_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold_field, "IField", field_id)
        gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", self.mesh_size_pluton)
        gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", self.mesh_size_far)
        gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", 1000)
        
        # Use the threshold field to define the mesh size
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Write mesh file
        gmsh.write(self.mesh_file)
        
        # Display mesh if GUI is enabled
        if show_gui:
            gmsh.fltk.run()
        
        # Finalize gmsh
        gmsh.finalize()
        
        print(f"Mesh generated and saved to {self.mesh_file}")
        return self.mesh_file
    
    def convert_mesh_to_fipy(self, mesh_file=None):
        """
        Convert gmsh mesh to a format suitable for FiPy
        """
        if mesh_file is None:
            mesh_file = self.mesh_file
            
        # Load mesh with meshio
        mesh = meshio.read(mesh_file)
        
        # Extract cells and points for FiPy
        # Implementation depends on the FiPy interface
        # This will be implemented later
        
        return mesh
    
    def define_geological_units(self, mesh):
        """
        Define geological units (pluton, host rock, sediments) on the mesh
        """
        # This will be implemented to create cell variables for FiPy
        # indicating which geological unit each cell belongs to
        pass
    
    def setup_fracture_zones(self, mesh, fracture_zones=None):
        """
        Set up fracture zones with enhanced permeability
        """
        # Define default fracture zones if none provided
        if fracture_zones is None:
            # Example: Two vertical fracture zones from pluton to surface
            fracture_zones = [
                {"x": -self.pluton_width/4, "width": 200, "permeability_factor": 10},
                {"x": self.pluton_width/4, "width": 200, "permeability_factor": 10}
            ]
        
        # Implementation will mark cells that belong to fracture zones
        # for permeability enhancement
        pass