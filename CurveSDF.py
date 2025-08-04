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

def tube_sdf(curve_obj, query_points, radius_func, start_radius, end_radius):
    """
    Compute the signed distance field for a tube with varying radius and capped ends.

    Parameters:
        curve_obj: geomdl.NURBS.Curve() - The centerline curve.
        query_points: numpy array of shape (N, 3) - The grid points.
        radius_func: function - A function that takes 't' and returns the radius.
        start_radius (float): Radius at the start of the tube (t=0).
        end_radius (float): Radius at the end of the tube (t=1).

    Returns:
        numpy.ndarray: Signed distance field values for the tube.
    """
    # 1. Get distances and t-parameters to the closest point on the curve
    distances, t_values = curve_sdf(curve_obj, query_points)

    # 2. Calculate varying radius for each query point
    radii_at_closest_points = radius_func(t_values, start_radius, end_radius)

    # 3. Calculate the signed distance for the main tube body
    # This is the distance to the curve minus the radius at that point
    sdf_body = distances - radii_at_closest_points

    # 4. Calculate SDFs for the end caps
    # Get start and end points of the curve
    curve_start_pt = np.array(curve_obj.evaluate_single(0.0))
    curve_end_pt = np.array(curve_obj.evaluate_single(1.0))

    # Get tangent vectors at start and end
    # geomdl.NURBS.Curve.tangent() method returns a list of [point, tangent_vector]
    start_eval_res = curve_obj.tangent(0.0)
    end_eval_res = curve_obj.tangent(1.0)
    
    # Ensure start_eval_res and end_eval_res are not None and have expected structure
    if start_eval_res and len(start_eval_res) > 1:
        start_tangent = np.array(start_eval_res[1])
    else:
        # Fallback if tangent evaluation fails or returns unexpected format
        # This might happen for very simple curves or specific configurations
        # For a robust solution, one might approximate or handle this more gracefully
        print("Warning: Could not get start tangent. Approximating.")
        start_tangent = np.array(curve_obj.evaluate_single(0.01)) - curve_start_pt
    
    if end_eval_res and len(end_eval_res) > 1:
        end_tangent = np.array(end_eval_res[1])
    else:
        print("Warning: Could not get end tangent. Approximating.")
        end_tangent = curve_end_pt - np.array(curve_obj.evaluate_single(0.99))
    
    # Ensure tangents are normalized for plane SDF calculation
    start_tangent = start_tangent / np.linalg.norm(start_tangent)
    end_tangent = end_tangent / np.linalg.norm(end_tangent)

    # SDF for the plane at the start (points "before" the start tangent are positive)
    sdf_cap_start = signed_plane_sdf(query_points, curve_start_pt, -start_tangent)

    # SDF for the plane at the end (points "after" the end tangent are positive)
    sdf_cap_end = signed_plane_sdf(query_points, curve_end_pt, end_tangent)

    # 5. Combine SDFs: The tube is the intersection of the body and the two half-spaces defined by the caps.
    # We use np.maximum to implement intersection (min of positive values, max of negative values)
    # The body SDF is combined with the cap SDFs.
    # A point is inside the tube if it's inside the body AND inside both caps.
    # So, we want the max of (sdf_body, sdf_cap_start, sdf_cap_end)
    # This assumes the plane SDFs are defined such that positive is "outside" the desired volume
    # and negative is "inside".
    
    # The tube is defined by:
    # 1. Being "inside" the cylinder (sdf_body <= 0)
    # 2. Being "after" the start plane (sdf_cap_start <= 0)
    # 3. Being "before" the end plane (sdf_cap_end <= 0)
    
    # If sdf_cap_start is positive, it means the point is "outside" the start cap.
    # If sdf_cap_end is positive, it means the point is "outside" the end cap.
    
    # To combine, we want the maximum of the "outside" distances.
    # If any of these is positive, the point is outside the final shape.
    # If all are negative, the point is inside.
    
    # Combined SDF for the watertight tube
    sdf_combined = np.maximum(sdf_body, np.maximum(sdf_cap_start, sdf_cap_end))

    return sdf_combined

def radius_function(t, start_radius, end_radius):
    """
    Defines a linearly varying radius along the curve's t-parameter.
    
    Parameters:
        t (float or np.ndarray): Parameter along the curve (0.0 to 1.0).
        start_radius (float): Radius at t=0.0.
        end_radius (float): Radius at t=1.0.
        
    Returns:
        float or np.ndarray: The interpolated radius.
    """
    return start_radius + (end_radius - start_radius) * t

def signed_plane_sdf(points, plane_origin, plane_normal):
    """
    Computes the signed distance from points to a plane.
    
    Parameters:
        points (np.ndarray): Array of points (N, 3).
        plane_origin (np.ndarray): A point on the plane (3,).
        plane_normal (np.ndarray): The normal vector of the plane (3,).
        
    Returns:
        np.ndarray: Signed distances from each point to the plane.
    """
    plane_normal = plane_normal / np.linalg.norm(plane_normal) # Normalize normal
    return np.dot(points - plane_origin, plane_normal)

def curve_sdf(curve_obj, query_points):
    """
    Compute the distance from specific points to the nearest point on the NURBS curve
    and return the 't' parameter of that closest point.

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        query_points: numpy array of shape (N, 3) or single point [x, y, z]

    Returns:
        distances: numpy array of distances from each query point to the curve
        t_values: numpy array of 't' parameters corresponding to the closest points on the curve
    """
    # Convert input to numpy array if it's a single point
    if not isinstance(query_points, np.ndarray):
        query_points = np.array(query_points)
    if query_points.ndim == 1:
        query_points = query_points.reshape(1, -1)
    
    # Generate a finer set of curve points and their corresponding 't' values
    # The 'delta' property controls the step size for evaluation
    # We need to manually generate t_values as evalpts doesn't store them directly
    num_eval_points = int(1.0 / curve_obj.delta) + 1
    t_params = np.linspace(0.0, 1.0, num_eval_points)
    
    curve_points_with_t = []
    for t in t_params:
        pt = curve_obj.evaluate_single(t)
        curve_points_with_t.append(pt)
    
    curve_points_np = np.array(curve_points_with_t, dtype=float)
    
    # Build KDTree for curve points
    leaf_size = min(20, len(curve_points_np) // 2) if len(curve_points_np) > 0 else 1
    kdtree = cKDTree(curve_points_np, leafsize=leaf_size)
    
    # Query distances and indices of closest points
    distances, indices = kdtree.query(query_points, workers=-1)
    
    # Get the 't' values corresponding to the closest points
    t_values = t_params[indices]
    
    # Return scalar if single point, array otherwise
    if len(query_points) == 1:
        return float(distances[0]), float(t_values[0])
    else:
        return distances, t_values

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
    # # --- NURBS Curve Definition using geomdl ---
    # curve = NURBS.Curve()
    # curve.degree = 3
    # curve.ctrlpts = [(-1, -1, -1), (0, -1, 0), (0, 1, 0), (1, 1, 1)]
    # curve.weights = [1, 10, 10, 1]
    # curve.knotvector = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    # curve.delta = 0.001 # Set the delta for evaluation
    # curve.evaluate()    # Precompute evaluation points
    
    # trunk = NURBS.Curve()
    # trunk.degree = 3
    # trunk.ctrlpts = [(-0.5, 0.25, -0.5), (0, 0, -0.5), (0, 0, 0), (0, 0, 0.5), (0.5, 0.25, 0.5)]
    # trunk.weights = [1.0, 2.0, 2.0, 2.0, 1.0]
    # trunk.knotvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    # trunk.delta = 0.001
    # trunk.evaluate()
    
    # curves=[curve, trunk]
    
    # # ctrlpts = [
    # # [[1.0, 0.0, 0.0, 1.0], [0.7071, 0.7071, 0.0, 0.7071], [0.0, 1.0, 0.0, 1.0], [-0.7071, 0.7071, 0.0, 0.7071], [-1.0, 0.0, 0.0, 1.0], [-0.7071, -0.7071, 0.0, 0.7071], [0.0, -1.0, 0.0, 1.0], [0.7071, -0.7071, 0.0, 0.7071], [1.0, 0.0, 0.0, 1.0]],
    # # [[1.0, 0.0, 1.0, 1.0], [0.7071, 0.7071, 0.7071, 0.7071], [0.0, 1.0, 1.0, 1.0], [-0.7071, 0.7071, 0.7071, 0.7071], [-1.0, 0.0, 1.0, 1.0], [-0.7071, -0.7071, 0.7071, 0.7071], [0.0, -1.0, 1.0, 1.0], [0.7071, -0.7071, 0.7071, 0.7071], [1.0, 0.0, 1.0, 1.0]]
    # # ]
    # # # --- Create a NURBS surface ---
    # # surf = NURBS.Surface()
    # # surf.degree_u = 1
    # # surf.degree_v = 2
    # # surf.ctrlpts2d = ctrlpts
    # # surf.knotvector_u = [0, 0, 1, 1]
    # # surf.knotvector_v = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
    # # surf.delta = 0.005
    # # surf.evaluate()
    
    # Generate grid points
    bounds = (-1.5, 1.5, -1.5, 1.5, -1.5, 1.5)  # Define bounds for the grid
    resolution = 64  # Define resolution for the grid
    points, (X, Y, Z) = generate_grid_points(bounds, resolution)
    
    # # start_time = time.time()
    # # SDF = curve_sdf(surf, points)
    # # plot_surf_sdf(SDF, bounds, resolution)
    # # end_time = time.time()
    # # print(f"Total computation time: {end_time - start_time:.2f} seconds")
    
    # start_time = time.time()
    # SDF1_dist, _ = curve_sdf(curve, points)
    # SDF2_dist, _ = curve_sdf(trunk, points)
    # U_SDF = union_sdf([SDF1_dist, SDF2_dist])  # Combine SDFs using union operation
    # S_SDF = subtract_sdf([SDF1_dist, SDF2_dist])  # Combine SDFs using subtraction operation
    # I_SDF = intersect_sdf([SDF1_dist, SDF2_dist])  # Combine SDFs using intersection operation
    # plot_curve_sdf(U_SDF, curves, bounds, resolution, X, Y, Z)
    # print("U_SDF (union of original curves):", U_SDF)
    # end_time = time.time()
    # print(f"Total computation time: {end_time - start_time:.2f} seconds")
    
    # # --- Example usage of curve_sdf to compute the distance of a point to the curve ---
    # point = [-1,-1,-1]
    # # The curve_sdf now returns distance and t_value, so unpack them
    # dist, t_val = curve_sdf(curve, point)
    # print(f"Distance from point {point} to curve: {dist:.4f}, t-value: {t_val:.4f}")

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
    
    start_time = time.time()
    tube_sdf_values = tube_sdf(tube_curve, points, radius_function, start_radius_val, end_radius_val)
    plot_curve_sdf(tube_sdf_values, [tube_curve], bounds, resolution, X, Y, Z)
    end_time = time.time()
    print(f"Tube SDF computation and plotting time: {end_time - start_time:.2f} seconds")


if __name__ == '__main__':
    main()
