# 2D Viscous Finite Volume Solver for Structured Quadrilateral Meshes

> **Note**
> Some of the code and documentation in this repository was generated with ChatGPT.

## calc_dU

This code is an implementation of a numerical solver called "Roe's approximate Riemann solver" for two-dimensional problems. The solver is used to calculate the change in the primitive variables (delta-U) at the interfaces between cells in a mesh, given the gradient of the primitive variables in each element (dUdX) and the vector from the cell center to each cell edge (dXE).

The input matrices dUdX and dXE are expected to have the following dimensions:

    dUdX is a 3D matrix of size n_fields x n_dims x n_cells_tot, where n_fields is the number of fields (primitive variables), n_dims is the number of dimensions (2 in this case), and n_cells_tot is the total number of cells in the mesh (including ghost cells).
    dXE is a 2D matrix of size n_dims x n_faces x n_cells, where n_faces is the number of faces (edges) per cell and n_cells is the number of non-ghost cells in the mesh.

The output matrix dU will have dimensions n_fields x n_faces x n_cells, and will contain the change in the primitive variables at the interfaces between cells.

The solver works by looping over all cells, faces, fields, and dimensions, and calculating the change in each field at each face by summing the products of the gradient of the field in the element and the vector from the cell center to the cell edge.

The final output matrix dU will contain the calculated change in the primitive variables at each face of each cell.

### ghost cells

In the context of finite volume methods, ghost cells are additional cells that are added to the computational mesh to facilitate the implementation of boundary conditions. These cells are typically located just outside the physical domain of the problem and are used to store the values of the primitive variables at the boundary.

In this case, it appears that ghost cells are not included in the input matrix dXE, but are included in the input matrix dUdX. This means that the solver is only used to calculate the change in the primitive variables at the interfaces between non-ghost cells, but the input data for the solver (the gradient of the primitive variables) includes values from the ghost cells as well.

Ghost cells can be useful in situations where the boundary conditions are complex and cannot be easily implemented using other methods. They allow the solver to use the same numerical method to calculate the change in the primitive variables at both the internal interfaces between cells and the external interfaces between cells and the boundary.
