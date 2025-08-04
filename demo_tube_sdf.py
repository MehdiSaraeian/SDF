import time
import numpy as np
from geomdl import NURBS
from sdf_utils import generate_grid_points, tube_sdf, radius_function, plot_curve_sdf

def main():
    """
    Main function to demonstrate the tube_sdf functionality.
    """
    # --- Demonstrate tube_sdf ---
    print("\nDemonstrating tube_sdf...")
    tube_curve = NURBS.Curve()
    tube_curve.degree = 3
    tube_curve.ctrlpts = [(0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0)]
    tube_curve.knotvector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    tube_curve.delta = 0.001
    tube_curve.evaluate()

    start_radius_val = 0.5
    end_radius_val = 0.05

    # Generate grid points
    bounds = (-1.5, 1.5, -1.5, 1.5, -1.5, 1.5)  # Define bounds for the grid
    resolution = 64  # Define resolution for the grid
    points, (X, Y, Z) = generate_grid_points(bounds, resolution)

    start_time = time.time()
    tube_sdf_values = tube_sdf(tube_curve, points, radius_function, start_radius_val, end_radius_val)
    plot_curve_sdf(tube_sdf_values, [tube_curve], bounds, resolution, X, Y, Z)
    end_time = time.time()
    print(f"Tube SDF computation and plotting time: {end_time - start_time:.2f} seconds")

if __name__ == '__main__':
    main()
