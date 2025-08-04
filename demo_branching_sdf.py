import time
import numpy as np
from geomdl import NURBS
from sdf_utils import generate_grid_points, tube_sdf, radius_function, plot_curve_sdf, union_sdf



def create_branch_curve(trunk_curve, branch_point_param, branch_direction, branch_length):
    """
    Create a branch curve that starts at a point on the trunk curve.

    Parameters:
        trunk_curve: geomdl.NURBS.Curve - The parent curve
        branch_point_param: float - Parameter value on trunk where branch connects
        branch_direction: numpy array - Direction vector for the branch
        branch_length: float - Length of the branch curve

    Returns:
        geomdl.NURBS.Curve - The branch curve
    """
    # Evaluate trunk at branch point
    point = np.array(trunk_curve.evaluate_single(branch_point_param))

    # Normalize direction vector
    direction = np.array(branch_direction)
    direction = direction / np.linalg.norm(direction)

    # Create perpendicular vectors for more interesting curve shapes
    # Find two vectors perpendicular to the direction vector
    if abs(direction[0]) < abs(direction[1]) and abs(direction[0]) < abs(direction[2]):
        perp1 = np.array([0, -direction[2], direction[1]])
    elif abs(direction[1]) < abs(direction[0]) and abs(direction[1]) < abs(direction[2]):
        perp1 = np.array([-direction[2], 0, direction[0]])
    else:
        perp1 = np.array([-direction[1], direction[0], 0])

    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # Create branch curve
    branch = NURBS.Curve()
    branch.degree = 3

    # Create 5 control points for the branch (cubic curve)
    # First control point is at the branch connection
    # Use perpendicular vectors to create curved shape
    factor = branch_length * 0.25  # Scale factor for control point displacement

    ctrlpts = [
        point.tolist(),
        (point + 0.2 * branch_length * direction + factor * perp1).tolist(),
        (point + 0.5 * branch_length * direction + factor * perp2).tolist(),
        (point + 0.8 * branch_length * direction - factor * perp1).tolist(),
        (point + branch_length * direction).tolist()
    ]

    branch.ctrlpts = ctrlpts
    branch.weights = [1.0, 1.0, 1.0, 1.0, 1.0]  # Equal weights

    # Clamped knot vector for cubic curve with 5 control points
    branch.knotvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    branch.delta = 0.001  # Set evaluation delta
    branch.evaluate()

    return branch


if __name__ == '__main__':

    # --- Create trunk NURBS curve ---
    trunk = NURBS.Curve()
    trunk.degree = 3
    trunk.ctrlpts = [(-0.5, 0.25, -0.5), (0, 0, -0.5), (0, 0, 0), (0, 0, 0.5), (0.5, 0.25, 0.5)]
    trunk.weights = [1.0, 2.0, 2.0, 2.0, 1.0]
    trunk.knotvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    trunk.delta = 0.001  # Set the delta for evaluation
    trunk.evaluate()

    # Generate grid points

    bounds = (-1.5, 1.5, -1.5, 1.5, -1.5, 1.5)

    resolution = 128
    points, (X, Y, Z) = generate_grid_points(bounds, resolution)

    # --- Create two branch curves ---
    # Branch point is at parameter 0.5 (middle of the trunk)
    branch_point_param_1 = 0.5
    branch_point_param_2 = 0.5

    # First branch going upward
    branch1 = create_branch_curve(trunk, 
                                 branch_point_param_1, 
                                 branch_direction=[0, 1, 0.8], 
                                 branch_length=0.8)
    
    # Second branch going downward
    branch2 = create_branch_curve(trunk, 
                                 branch_point_param_2, 
                                 branch_direction=[0, -0.5, 1], 
                                 branch_length=0.9)
    
    curves = [trunk, branch1, branch2]
        
    # --- Calculate Tube SDFs ---
    # Define radii for the tubes
    trunk_start_radius = 0.2
    trunk_end_radius = 0.1
    branch_start_radius = 0.1
    branch_end_radius = 0.05

    # Calculate tube SDF for each curve
    trunk_sdf = tube_sdf(trunk, points, radius_function, trunk_start_radius, trunk_end_radius)
    branch1_sdf = tube_sdf(branch1, points, radius_function, branch_start_radius, branch_end_radius)
    branch2_sdf = tube_sdf(branch2, points, radius_function, branch_start_radius, branch_end_radius)
    
    sdfs = [trunk_sdf, branch1_sdf, branch2_sdf]

    # Perform union operation on all SDFs
    combined_sdf = union_sdf(sdfs)
    
    # Plot the combined SDF
    start_time = time.time()
    plot_curve_sdf(combined_sdf, curves, bounds, resolution, X, Y, Z)
    end_time = time.time()
    print(f"Total computation time: {end_time - start_time:.2f} seconds")