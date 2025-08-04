import time, math
import numpy as np
import plotly.graph_objects as go
from geomdl import NURBS
from scipy.spatial import cKDTree # type: ignore


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

def union_sdf(sdfs):
    """
    Perform union operation on multiple SDFs.
    
    Parameters:
        sdfs: List of numpy arrays containing SDF values
    
    Returns:
        combined_sdf: Numpy array with unified SDF values
    """
    if not sdfs:
        return None
    
    # Start with first SDF
    combined_sdf = sdfs[0].copy()
    
    # Combine with remaining SDFs
    for i in range(1, len(sdfs)):
        combined_sdf = np.minimum(combined_sdf, sdfs[i])
    
    return combined_sdf

def subtract_sdf(sdfs):
    """
    Perform subtraction (difference) operation on multiple SDFs.
    The first SDF is the base; subsequent SDFs are subtracted from it.
    
    Parameters:
        sdfs: List of numpy arrays containing SDF values
    
    Returns:
        result_sdf: Numpy array with subtracted SDF values
    """
    if not sdfs:
        return None
    result_sdf = sdfs[0].copy()
    for sdf in sdfs[1:]:
        result_sdf = np.maximum(result_sdf, -sdf) # Corrected subtraction operation
    return result_sdf

def intersect_sdf(sdfs):
    """
    Perform intersection operation on multiple SDFs.
    This computes the SDF of the intersection by taking the maximum
    among the provided SDF values.

    Parameters:
        sdfs: List of numpy arrays containing SDF values

    Returns:
        result_sdf: Numpy array representing the intersection SDF.
    """
    if not sdfs:
        return None
    result_sdf = sdfs[0].copy()
    for sdf in sdfs[1:]:
        result_sdf = np.maximum(result_sdf, sdf)
    return result_sdf

def curve_sdf(curve_obj, points):
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

def plot_curve_sdf(sdf, curves, bounds, resolution, X, Y, Z):
    """
    Plot the signed distance field of a curve (3D).

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        bounds (tuple): Plot bounds for the SDF grid (x_min, x_max, y_min, y_max, z_min, z_max)
        resolution (int): Grid resolution for SDF calculation
    """
    
    # Reshape distances back to grid
    SDF = sdf.reshape((resolution, resolution, resolution))

    # Plot results
    data =[]
    for curve in curves:
        curve_points = np.array(curve.evalpts, dtype=float)
        curve_trace = go.Scatter3d(
            x=curve_points[:, 0],
            y=curve_points[:, 1],
            z=curve_points[:, 2],
            mode='lines',
            line=dict(color='blue', width=5),
            name='NURBS Curve'
        )
        data.append(curve_trace)
    
    # iso_surface_trace = go.Isosurface(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=SDF.flatten(),
    #     isomin=0.01,
    #     # isomax=SDF.max(),
    #     isomax=0.15,
    #     opacity=0.5,
    #     surface_count=5,
    #     caps=dict(x_show=False, y_show=False, z_show=False),
    #     colorscale='viridis',
    #     colorbar=dict(title='Distance'),
    #     name='SDF Isosurfaces'
    # )
    # data.append(iso_surface_trace)
    # fig = go.Figure(data=data)
    
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
    data.append(volume_trace)
    fig = go.Figure(data=data)
    
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
    return None

def plot_surf_sdf(sdf, bounds, resolution, X, Y, Z):
    """
    Plot the signed distance field of a curve (3D).

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        bounds (tuple): Plot bounds for the SDF grid (x_min, x_max, y_min, y_max, z_min, z_max)
        resolution (int): Grid resolution for SDF calculation
    """
    
    # Reshape distances back to grid
    SDF = sdf.reshape((resolution, resolution, resolution))

    # Plot results
    data =[]

    iso_surface_trace = go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=SDF.flatten(),
        isomin=0.0,
        # isomax=SDF.max(),
        isomax=0.1,
        opacity=0.1,
        surface_count=3,
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorscale='viridis',
        colorbar=dict(title='Distance'),
        name='SDF Isosurfaces'
    )
    data.append(iso_surface_trace)
    fig = go.Figure(data=data)
    
    # volume_trace = go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=SDF.flatten(),
    #     opacity=0.1,
    #     surface_count=20,
    #     colorscale='Viridis',
    #     colorbar=dict(title='Distance'),
    #     name='SDF Volume'
    # )
    # data.append(volume_trace)
    # fig = go.Figure(data=data)
    
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
    return None

def main():
    # --- NURBS Curve Definition using geomdl ---
    curve = NURBS.Curve()
    curve.degree = 3
    curve.ctrlpts = [(-1, -1, -1), (0, -1, 0), (0, 1, 0), (1, 1, 1)]
    curve.weights = [1, 10, 10, 1]
    curve.knotvector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    curve.delta = 0.001 # Set the delta for evaluation
    curve.evaluate()    # Precompute evaluation points
    
    trunk = NURBS.Curve()
    trunk.degree = 3
    trunk.ctrlpts = [(-0.5, 0.25, -0.5), (0, 0, -0.5), (0, 0, 0), (0, 0, 0.5), (0.5, 0.25, 0.5)]
    trunk.weights = [1.0, 2.0, 2.0, 2.0, 1.0]
    trunk.knotvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    trunk.delta = 0.001
    trunk.evaluate()
    
    curves=[curve, trunk]
    
    # ctrlpts = [
    # [[1.0, 0.0, 0.0, 1.0], [0.7071, 0.7071, 0.0, 0.7071], [0.0, 1.0, 0.0, 1.0], [-0.7071, 0.7071, 0.0, 0.7071], [-1.0, 0.0, 0.0, 1.0], [-0.7071, -0.7071, 0.0, 0.7071], [0.0, -1.0, 0.0, 1.0], [0.7071, -0.7071, 0.0, 0.7071], [1.0, 0.0, 0.0, 1.0]],
    # [[1.0, 0.0, 1.0, 1.0], [0.7071, 0.7071, 0.7071, 0.7071], [0.0, 1.0, 1.0, 1.0], [-0.7071, 0.7071, 0.7071, 0.7071], [-1.0, 0.0, 1.0, 1.0], [-0.7071, -0.7071, 0.7071, 0.7071], [0.0, -1.0, 1.0, 1.0], [0.7071, -0.7071, 0.7071, 0.7071], [1.0, 0.0, 1.0, 1.0]]
    # ]
    # # --- Create a NURBS surface ---
    # surf = NURBS.Surface()
    # surf.degree_u = 1
    # surf.degree_v = 2
    # surf.ctrlpts2d = ctrlpts
    # surf.knotvector_u = [0, 0, 1, 1]
    # surf.knotvector_v = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
    # surf.delta = 0.005
    # surf.evaluate()
    
    # Generate grid points
    bounds = (-1.5, 1.5, -1.5, 1.5, -1.5, 1.5)  # Define bounds for the grid
    resolution = 64  # Define resolution for the grid
    points, (X, Y, Z) = generate_grid_points(bounds, resolution)
    
    # start_time = time.time()
    # SDF = curve_sdf(surf, points)
    # plot_surf_sdf(SDF, bounds, resolution)
    # end_time = time.time()
    # print(f"Total computation time: {end_time - start_time:.2f} seconds")
    
    start_time = time.time()
    SDF1 = curve_sdf(curve, points)
    SDF2 = curve_sdf(trunk, points)
    U_SDF = union_sdf([SDF1, SDF2])  # Combine SDFs using union operation
    S_SDF = subtract_sdf([SDF1, SDF2])  # Combine SDFs using subtraction operation
    I_SDF = intersect_sdf([SDF1, SDF2])  # Combine SDFs using intersection operation
    plot_curve_sdf(U_SDF, curves, bounds, resolution, X, Y, Z)
    print("S_SDF:", U_SDF)
    end_time = time.time()
    print(f"Total computation time: {end_time - start_time:.2f} seconds")
    
    # --- Example usage of curve_sdf to compute the distance of a point to the curve ---
    point = [-1,-1,-1]
    dist = curve_sdf(curve, point)
    print(f"Distance from point {point} to curve: {dist:.4f}")

if __name__ == '__main__':
    main()
