import time
import numpy as np
import plotly.graph_objects as go
from geomdl import NURBS
from scipy.spatial import cKDTree


def generate_grid_points(bounds, resolution):
    """Generate grid points efficiently using numpy operations."""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    # Generate points directly in the right shape
    x = np.linspace(x_min, x_max, resolution)
    y = np.linspace(y_min, y_max, resolution)
    z = np.linspace(z_min, z_max, resolution)
    
    # Create meshgrid
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Reshape points for KDTree query
    # Use numpy's ascontiguousarray for better memory access
    points = np.ascontiguousarray(np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1))
    return points, (X, Y, Z)

def compute_point_distance(curve_obj, points):
    """
    Compute the distance from specific points to the nearest point on the NURBS curve.

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        points: numpy array of shape (N, 3) or single point [x, y, z]

    Returns:
        distances: numpy array of distances from each point to the curve
    """
    # Convert input to numpy array if it's a single point
    if not isinstance(points, np.ndarray):
        points = np.array(points)
    if points.ndim == 1:
        points = points.reshape(1, -1)
    
    # Build KDTree for curve points
    curve_points = np.array(curve_obj.evalpts, dtype=float)
    leaf_size = min(20, len(curve_points) // 2)
    kdtree = cKDTree(curve_points, leafsize=leaf_size)
    
    # Query distances
    distances, _ = kdtree.query(points, workers=-1)
    
    # Return scalar if single point, array otherwise
    return float(distances[0]) if len(points) == 1 else distances

def plot_curve_sdf(curve_obj, bounds, resolution=50):
    """
    Plot the signed distance field of a curve (3D).

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        bounds (tuple): Plot bounds for the SDF grid (x_min, x_max, y_min, y_max, z_min, z_max)
        resolution (int): Grid resolution for SDF calculation
        num_processes (int, optional): Number of processes for parallel computation
    """
    # Generate grid points efficiently
    points, (X, Y, Z) = generate_grid_points(bounds, resolution)
    
    # Build KDTree once for all points with optimized leaf_size
    curve_points = np.array(curve_obj.evalpts, dtype=float)
    leaf_size = min(20, len(curve_points) // 2)  # Optimize leaf_size based on point count
    kdtree = cKDTree(curve_points, leafsize=leaf_size)
    
    print(f"Calculating SDF on a {resolution}x{resolution}x{resolution} grid...")
    
    # For this problem size, sequential execution with threaded KDTree is fastest
    distances, _ = kdtree.query(points, workers=-1)

    # Reshape distances back to grid
    SDF = distances.reshape((resolution, resolution, resolution))
    print("SDF calculation complete.")

    # Plot results
    curve_trace = go.Scatter3d(
        x=curve_points[:, 0],
        y=curve_points[:, 1],
        z=curve_points[:, 2],
        mode='lines',
        line=dict(color='blue', width=5),
        name='NURBS Curve'
    )
    
    # iso_surface_trace = go.Isosurface(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=SDF.flatten(),
    #     isomin=0.01,
    #     # isomax=SDF.max(),
    #     isomax=0.1,
    #     opacity=0.4,
    #     surface_count=5,
    #     caps=dict(x_show=False, y_show=False, z_show=False),
    #     colorscale='viridis',
    #     colorbar=dict(title='Distance'),
    #     name='SDF Isosurfaces'
    # )

    # fig = go.Figure(data=[curve_trace, iso_surface_trace])
    
    volume_trace = go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=SDF.flatten(),
        opacity=0.1,
        surface_count=20,
        colorscale='Viridis',
        colorbar=dict(title='Distance'),
        name='SDF Volume'
    )

    fig = go.Figure(data=[curve_trace, volume_trace])
    
    fig.update_layout(
        title='Interactive 3D Signed Distance Field (SDF) of NURBS Curve',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            xaxis=dict(range=[bounds[0], bounds[1]], autorange=False),
            yaxis=dict(range=[bounds[2], bounds[3]], autorange=False),
            zaxis=dict(range=[bounds[4], bounds[5]], autorange=False),
            aspectratio=dict(x=1, y=1, z=1)
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )
    print("Displaying interactive plot...")
    fig.show()
    return SDF

def demo_point_distance():
    """
    Demonstrate the point distance calculation with a simple example.
    """
    # Create a simple NURBS curve
    curve = NURBS.Curve()
    curve.degree = 3
    curve.ctrlpts = [(-1, -1, -1), (0, -1, 0), (0, 1, 0), (1, 1, 1)]
    curve.weights = [1, 10, 10, 1]
    curve.knotvector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    curve.delta = 0.001
    curve.evaluate()

    # Test points
    test_points = np.array([
        [0, 0, 0],      # Point near the curve
        [2, 2, 2],      # Point far from the curve
        [-1, -1, -1],   # Point at curve endpoint
    ])

    # Calculate and display distances
    distances = compute_point_distance(curve, test_points)
    
    print("\nPoint Distance Demo:")
    for point, dist in zip(test_points, distances):
        print(f"Distance from point {point} to curve: {dist:.4f}")
    return test_points, distances

if __name__ == '__main__':
    # --- NURBS Curve Definition using geomdl ---
    degree = 3
    ctrlpts = [(-1, -1, -1), (0, -1, 0), (0, 1, 0), (1, 1, 1)]
    weights = [1, 10, 10, 1]
    knotvector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    curve = NURBS.Curve()
    curve.degree = degree
    curve.ctrlpts = ctrlpts
    curve.weights = weights
    curve.knotvector = knotvector
    curve.delta = 0.001  # Set the delta for evaluation
    curve.evaluate()   # Precompute evaluation points

    bounds = (-1.5, 1.5, -1.5, 1.5, -1.5, 1.5)  # Define bounds for the plot
    
    start_time = time.time()
    # plot_curve_sdf(curve, bounds, resolution=50)
    demo_point_distance()
    end_time = time.time()
    print(f"Total computation time: {end_time - start_time:.2f} seconds")
