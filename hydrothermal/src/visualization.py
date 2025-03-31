"""
Visualization tools for hydrothermal simulation.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import os
import h5py
import pyvista as pv

class Visualization:
    def __init__(self, output_dir="output", vis_dir="visualization"):
        """
        Initialize visualization tools
        
        Args:
            output_dir: Directory with simulation outputs
            vis_dir: Directory to save visualizations
        """
        # Import configuration values
        from src.config import (TEMPERATURE_COLORMAP, VELOCITY_COLORMAP,
                              QUIVER_SCALE, QUIVER_DENSITY)
        
        self.output_dir = output_dir
        self.vis_dir = vis_dir
        os.makedirs(vis_dir, exist_ok=True)
        
        # Use configured colormap for temperature
        self.temp_cmap = getattr(cm, TEMPERATURE_COLORMAP)
        self.velocity_cmap = getattr(cm, VELOCITY_COLORMAP)
        
        # Use configured quiver parameters for velocity
        self.quiver_scale = QUIVER_SCALE
        self.quiver_density = QUIVER_DENSITY  # Use every nth point
    
    def list_output_files(self):
        """
        List all output files in chronological order
        
        Returns:
            list: List of output file paths
        """
        files = [f for f in os.listdir(self.output_dir) if f.endswith('.h5')]
        
        # Filter for hydrothermal simulation files only (exclude tracking_data.h5 etc.)
        sim_files = [f for f in files if f.startswith('hydrothermal_')]
        
        if sim_files:
            # Sort by time if possible
            try:
                sim_files.sort(key=lambda x: float(x.split('_')[1].split('yrs')[0]))
            except (IndexError, ValueError):
                # If sorting fails, use alphabetical sorting as fallback
                sim_files.sort()
        
        return [os.path.join(self.output_dir, f) for f in sim_files]
    
    def load_output_file(self, file_path):
        """
        Load an output file
        
        Args:
            file_path: Path to HDF5 output file
            
        Returns:
            dict: Data from the file
        """
        data = {}
        with h5py.File(file_path, 'r') as f:
            # Load mesh coordinates
            data['x'] = f['x'][:]
            data['y'] = f['y'][:]
            
            # Load fields
            data['temperature'] = f['temperature'][:]
            data['velocity_x'] = f['velocity_x'][:]
            data['velocity_y'] = f['velocity_y'][:]
            data['pressure'] = f['pressure'][:]
            data['geological_unit'] = f['geological_unit'][:]
            data['permeability'] = f['permeability'][:]
            
            # Load metadata
            data['time'] = f.attrs['time']
            data['time_years'] = f.attrs['time_years']
        
        return data
    
    def plot_temperature(self, data, output_file=None, show=True):
        """
        Plot temperature field
        
        Args:
            data: Data dictionary from load_output_file
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(12, 8))
        
        # Create 2D temperature map
        plt.tripcolor(data['x'], data['y'], data['temperature'], 
                      cmap=self.temp_cmap, shading='gouraud')
        
        # Add contour lines
        contour_levels = np.linspace(data['temperature'].min(), data['temperature'].max(), 10)
        plt.tricontour(data['x'], data['y'], data['temperature'], 
                       levels=contour_levels, colors='k', linewidths=0.5, alpha=0.7)
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('Temperature (°C)')
        
        # Set labels and title
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title(f'Temperature Distribution at {data["time_years"]:.0f} years')
        
        # Equal aspect ratio
        plt.axis('equal')
        plt.grid(True, alpha=0.3)
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved temperature plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_temperature_with_velocity(self, data, output_file=None, show=True):
        """
        Plot temperature field with velocity vectors
        
        Args:
            data: Data dictionary from load_output_file
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(12, 8))
        
        # Create 2D temperature map
        plt.tripcolor(data['x'], data['y'], data['temperature'], 
                      cmap=self.temp_cmap, shading='gouraud')
        
        # For now, skip the velocity vectors since we're having issues with them
        # Just plot the temperature field
        
        # Add temperature colorbar
        cbar_temp = plt.colorbar()
        cbar_temp.set_label('Temperature (°C)')
        
        # Set labels and title
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title(f'Temperature Field at {data["time_years"]:.0f} years')
        
        # Equal aspect ratio
        plt.axis('equal')
        plt.grid(True, alpha=0.3)
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved temperature plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_geological_units(self, data, output_file=None, show=True):
        """
        Plot geological units
        
        Args:
            data: Data dictionary from load_output_file
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(12, 8))
        
        # Create 2D map of geological units
        unit_map = plt.tripcolor(data['x'], data['y'], data['geological_unit'], 
                                 cmap='viridis', shading='flat')
        
        # Add colorbar with custom ticks for geological units
        cbar = plt.colorbar(ticks=[1, 2, 3])
        cbar.ax.set_yticklabels(['Pluton', 'Host Rock', 'Sediments'])
        
        # Set labels and title
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('Geological Units')
        
        # Equal aspect ratio
        plt.axis('equal')
        plt.grid(True, alpha=0.3)
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved geological units plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_permeability(self, data, output_file=None, show=True):
        """
        Plot permeability field (logarithmic scale)
        
        Args:
            data: Data dictionary from load_output_file
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(12, 8))
        
        # Create 2D permeability map (log scale)
        perm_map = plt.tripcolor(data['x'], data['y'], data['permeability'], 
                                 norm=LogNorm(vmin=1e-18, vmax=1e-14),
                                 cmap='plasma', shading='flat')
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('Permeability (m²)')
        
        # Set labels and title
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('Permeability Distribution')
        
        # Equal aspect ratio
        plt.axis('equal')
        plt.grid(True, alpha=0.3)
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved permeability plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_heat_flux(self, data, output_file=None, show=True):
        """
        Plot heat flux vectors
        
        Args:
            data: Data dictionary from load_output_file
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        # This requires computing conductive and advective heat fluxes
        # Will implement in the full version
        pass
    
    def create_temperature_animation(self, output_file=None):
        """
        Create an animation of temperature evolution
        
        Args:
            output_file: Output file for the animation
        """
        # This would create an animation by iterating through timesteps
        # Will implement in the full version
        pass
    
    def create_pyvista_visualization(self, data):
        """
        Create 3D visualization using PyVista
        
        Args:
            data: Data dictionary from load_output_file
        """
        # Create a PyVista mesh from the data
        mesh = pv.PolyData(np.column_stack((data['x'], data['y'], np.zeros_like(data['x']))))
        
        # Add temperature as a point scalar
        mesh.point_data['Temperature'] = data['temperature']
        
        # Add velocity as a point vector
        velocity = np.column_stack((data['velocity_x'], data['velocity_y'], np.zeros_like(data['x'])))
        mesh.point_data['Velocity'] = velocity
        
        # Create plotter
        plotter = pv.Plotter()
        
        # Add temperature surface
        plotter.add_mesh(mesh, scalars='Temperature', cmap='viridis', point_size=10, render_points_as_spheres=True)
        
        # Add velocity glyphs
        glyphs = mesh.glyph(
            orient='Velocity',
            scale='Velocity',
            factor=0.02,
            geom=pv.Arrow()
        )
        plotter.add_mesh(glyphs, color='black')
        
        # Set camera position
        plotter.view_xy()
        
        # Add a title
        plotter.add_title(f"Temperature at {data['time_years']:.0f} years", font_size=15)
        
        # Show the plot
        plotter.show()
    
    def plot_temperature_profiles(self, data_list, location_x=0, output_file=None, show=True):
        """
        Plot temperature profiles along vertical transect at different times
        
        Args:
            data_list: List of data dictionaries for different times
            location_x: x-coordinate for the vertical profile
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(10, 8))
        
        # Plot temperature vs. depth for each time step
        for data in data_list:
            # Find points close to the specified x-coordinate
            tolerance = 100  # meters
            mask = np.abs(data['x'] - location_x) < tolerance
            
            if np.sum(mask) > 0:
                # Extract y-coordinates and temperatures along the profile
                y_profile = data['y'][mask]
                temp_profile = data['temperature'][mask]
                
                # Sort by depth
                sort_idx = np.argsort(y_profile)
                y_profile = y_profile[sort_idx]
                temp_profile = temp_profile[sort_idx]
                
                # Plot the profile
                plt.plot(temp_profile, y_profile, label=f"{data['time_years']:.0f} years")
        
        # Set labels and title
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Depth (m)')
        plt.title(f'Vertical Temperature Profile at x={location_x}m')
        
        # Adjust y-axis direction to show depth properly
        plt.gca().invert_yaxis()
        
        # Add legend and grid
        plt.legend()
        plt.grid(True)
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved temperature profile plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_time_series(self, tracking_data, output_file=None, show=True):
        """
        Plot time series data from tracking
        
        Args:
            tracking_data: Dictionary with tracking data
            output_file: Optional file path to save the plot
            show: Whether to display the plot
        """
        plt.figure(figsize=(12, 8))
        
        # Convert times to years
        times_years = tracking_data['times'] / (3600 * 24 * 365)
        
        # Plot maximum velocity vs. time
        plt.subplot(2, 1, 1)
        plt.plot(times_years, tracking_data['max_velocities'], 'b-', linewidth=2)
        plt.xlabel('Time (years)')
        plt.ylabel('Maximum Velocity (m/s)')
        plt.title('Maximum Fluid Velocity Over Time')
        plt.grid(True)
        
        # Plot heat flux if available
        if 'heat_flux_history' in tracking_data:
            plt.subplot(2, 1, 2)
            plt.plot(times_years, tracking_data['heat_flux_history'], 'r-', linewidth=2)
            plt.xlabel('Time (years)')
            plt.ylabel('Heat Extraction Rate (W)')
            plt.title('Heat Extraction Rate from Pluton')
            plt.grid(True)
        
        plt.tight_layout()
        
        # Save if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved time series plot to {output_file}")
        
        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()
    
    def generate_report(self):
        """
        Generate a comprehensive visualization report
        """
        # List all output files
        files = self.list_output_files()
        
        if not files:
            print("No output files found to visualize.")
            return
        
        # Create report directory
        report_dir = os.path.join(self.vis_dir, "report")
        os.makedirs(report_dir, exist_ok=True)
        
        # Load first time step to show geological units and permeability
        initial_data = self.load_output_file(files[0])
        
        # Plot geological units
        self.plot_geological_units(initial_data, 
                                  output_file=os.path.join(report_dir, "geological_units.png"),
                                  show=False)
        
        # Plot permeability
        self.plot_permeability(initial_data,
                              output_file=os.path.join(report_dir, "permeability.png"),
                              show=False)
        
        # Generate temperature plots for each time step
        for file in files:
            data = self.load_output_file(file)
            time_years = int(data['time_years'])
            
            # Plot temperature
            self.plot_temperature(data,
                                 output_file=os.path.join(report_dir, f"temperature_{time_years}yrs.png"),
                                 show=False)
            
            # Plot temperature with velocity
            self.plot_temperature_with_velocity(data,
                                              output_file=os.path.join(report_dir, f"temp_vel_{time_years}yrs.png"),
                                              show=False)
        
        # Load tracking data
        tracking_file = os.path.join(self.output_dir, "tracking_data.h5")
        if os.path.exists(tracking_file):
            with h5py.File(tracking_file, 'r') as f:
                tracking_data = {
                    'times': f['times'][:],
                    'max_velocities': f['max_velocities'][:]
                }
                
                # Load heat flux if available
                if 'heat_flux_history' in f:
                    tracking_data['heat_flux_history'] = f['heat_flux_history'][:]
            
            # Plot time series
            self.plot_time_series(tracking_data,
                                 output_file=os.path.join(report_dir, "time_series.png"),
                                 show=False)
        
        # Load data for selected time steps for profiles
        key_times = [0, len(files)//4, len(files)//2, 3*len(files)//4, -1]
        data_list = [self.load_output_file(files[i % len(files)]) for i in key_times]
        
        # Plot temperature profiles
        self.plot_temperature_profiles(data_list,
                                      output_file=os.path.join(report_dir, "temperature_profiles.png"),
                                      show=False)
        
        print(f"Generated visualization report in {report_dir}")

def create_paraview_visualization(input_file, output_file=None):
    """
    Create a ParaView visualization from HDF5 output
    
    Args:
        input_file: Input HDF5 file
        output_file: Optional output file for ParaView state
    """
    # This would create a ParaView state file or directly invoke ParaView
    # Will be implemented in the full version
    pass