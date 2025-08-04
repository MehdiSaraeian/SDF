import time
import numpy as np
from geomdl import NURBS
from sdf_utils import generate_grid_points, curve_sdf, union_sdf


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

def plot(sdf, curves, bounds, resolution, X, Y, Z):
    """
    Plot the signed distance field of a curve (3D).

    Parameters:
        curve_obj: geomdl.NURBS.Curve()
        bounds (tuple): Plot bounds for the SDF grid (x_min, x_max, y_min, y_max, z_min, z_max)
        resolution (int): Grid resolution for SDF calculation
    """
    import plotly.graph_objects as go

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

    iso_surface_trace = go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=SDF.flatten(),
        isomin=0.01,
        # isomax=SDF.max(),
        isomax=0.05,
        opacity=0.1,
        surface_count=5,
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

    # fig = go.Figure(data=[curve_trace, volume_trace])

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
    bounds = (-1, 1, -1, 1, -1, 1)
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
                                 branch_length=0.5)

    # Second branch going downward
    branch2 = create_branch_curve(trunk,
                                 branch_point_param_2,
                                 branch_direction=[0, -0.5, 1],
                                 branch_length=0.7)

    curves = [trunk, branch1, branch2]

    # Calculate SDFs
    trunk_sdf, _ = curve_sdf(trunk, points)
    branch1_sdf, _ = curve_sdf(branch1, points)
    branch2_sdf, _ = curve_sdf(branch2, points)
    sdfs = [trunk_sdf, branch1_sdf, branch2_sdf]

    # Perform union operation on all SDFs
    combined_sdf = union_sdf(sdfs)

    # Plot the combined SDF
    start_time = time.time()
    sdf = plot(combined_sdf, curves, bounds, resolution, X, Y, Z)
    end_time = time.time()
    print(f"Total computation time: {end_time - start_time:.2f} seconds")
