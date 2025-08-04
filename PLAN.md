# Plan for Modularization of SDF Codebase

This document outlines the plan for modularizing the SDF codebase by introducing a Python package structure, using an `__init__.py` file, and reorganizing existing and new demonstration scripts.

## Goals

*   Create a clean, organized, and extensible package for SDF-related functionalities.
*   Separate core SDF functions from demonstration scripts.
*   Ensure all existing and new functionalities work correctly after modularization.

## Detailed Steps

### Step 1: Create a New Package Directory

*   Create a new top-level directory named `sdf_utils`. This directory will serve as the root for our new Python package.

### Step 2: Move and Rename `CurveSDF.py`

*   Move the existing `CurveSDF.py` file into the newly created `sdf_utils` directory.
*   Rename `CurveSDF.py` to `core.py`. This file will contain the core SDF calculation and utility functions.

### Step 3: Create `__init__.py`

*   Create an empty file named `__init__.py` inside the `sdf_utils` directory. This file is essential for Python to recognize `sdf_utils` as a package.
*   Add import statements to `__init__.py` to expose necessary functions from `core.py` directly under the `sdf_utils` namespace. For example:
    ```python
    # sdf_utils/__init__.py
    from .core import (
        generate_grid_points,
        union_sdf,
        subtract_sdf,
        intersect_sdf,
        tube_sdf,
        radius_function,
        signed_plane_sdf,
        curve_sdf,
        plot_curve_sdf,
        plot_surf_sdf
    )
    ```

### Step 4: Refactor `core.py` (formerly `CurveSDF.py`)

*   **Remove `main` Function**: The `if __name__ == '__main__': main()` block and the `main()` function itself will be removed from `core.py`. This file should now only contain reusable functions.
*   **Review Imports**: Ensure all imports within `core.py` are absolute for external libraries (e.g., `import numpy as np`, `from geomdl import NURBS`, `from scipy.spatial import cKDTree`) and that there are no remaining relative imports that would break after the move.

### Step 5: Adapt Existing Demo (`demo_branching_sdf.py`)

*   Update the import statements in `demo_branching_sdf.py` to reflect the new package structure. For example, instead of `from CurveSDF import generate_grid_points, ...`, it will become `from sdf_utils import generate_grid_points, ...`.

### Step 6: Create New Demo for Tube SDF (`demo_tube_sdf.py`)

*   Create a new Python file named `demo_tube_sdf.py` at the root level of the project.
*   Move the demonstration code for `tube_sdf` (currently located in the `main` function of `CurveSDF.py`, specifically the section under `--- Demonstrate tube_sdf ---`) into this new file.
*   Add necessary imports to `demo_tube_sdf.py` from the `sdf_utils` package (e.g., `from sdf_utils import NURBS, generate_grid_points, tube_sdf, radius_function, plot_curve_sdf`).
*   Include a `if __name__ == '__main__':` block to run the demonstration when the script is executed.

### Step 7: Testing

*   After completing the refactoring, execute both `demo_branching_sdf.py` and `demo_tube_sdf.py` to verify that all functionalities work correctly and that imports are resolved as expected.

### Step 8: Documentation (Optional but Recommended)

*   Consider creating a `README.md` file within the `sdf_utils` package directory or updating the main project `README.md` to explain the new modular structure, the purpose of each module, and how to use the exposed functions.

This modularization will enhance the project's maintainability, readability, and scalability.
