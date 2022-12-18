# 2D Viscous Finite Volume Solver for Structured Quadrilateral Meshes

> **Note**
> Some of the code and documentation in this repository was generated with ChatGPT.

## calc_dU

This code is an implementation of a numerical solver called "Roe's approximate Riemann solver" for two-dimensional problems. The solver is used to calculate the change in the primitive variables (delta-U) at the interfaces between cells in a mesh, given the gradient of the primitive variables in each element (dUdX) and the vector from the cell center to each cell edge (dXE).

The input matrices `dUdX` and `dXE` are expected to have the following dimensions:

- `dUdX` is a 3D matrix of size `n_fields x n_dims x n_cells_tot`, where `n_fields` is the number of fields (primitive variables), `n_dims` is the number of dimensions (2 in this case), and `n_cells_tot` is the total number of cells in the mesh (including ghost cells).
- `dXE` is a 2D matrix of size `n_dims x n_faces x n_cells`, where `n_faces` is the number of faces (edges) per cell and `n_cells` is the number of non-ghost cells in the mesh.

The output matrix `dU` will have dimensions `n_fields x n_faces x n_cells`, and will contain the change in the primitive variables at the interfaces between cells.

The solver works by looping over all cells, faces, fields, and dimensions, and calculating the change in each field at each face by summing the products of the gradient of the field in the element and the vector from the cell center to the cell edge.

The final output matrix `dU` will contain the calculated change in the primitive variables at each face of each cell.

### ghost cells

In the context of finite volume methods, ghost cells are additional cells that are added to the computational mesh to facilitate the implementation of boundary conditions. These cells are typically located just outside the physical domain of the problem and are used to store the values of the primitive variables at the boundary.

In this case, it appears that ghost cells are not included in the input matrix dXE, but are included in the input matrix dUdX. This means that the solver is only used to calculate the change in the primitive variables at the interfaces between non-ghost cells, but the input data for the solver (the gradient of the primitive variables) includes values from the ghost cells as well.

Ghost cells can be useful in situations where the boundary conditions are complex and cannot be easily implemented using other methods. They allow the solver to use the same numerical method to calculate the change in the primitive variables at both the internal interfaces between cells and the external interfaces between cells and the boundary.


## calc_grad_dxinv

This Matlab function calculates the gradient of a set of primitive variables `U` in each cell of a grid. The gradient is computed using the QR method, which involves encircling each cell with a set of neighboring cells and using an inverse QR matrix `dXinv` to compute the gradient.

The input matrices are:

- `U`: a `n_cells_total x n_fields` matrix containing the solution vector in each cell, including ghost cells.
- `c2ac`: a `n_cells x max_nc` matrix containing a list of all cells encircling each interior cell.
- `c2nac`: a `n_cells x 1 matrix` containing the number of cells encircling each interior cell.
- `dXinv`: a `n_dims x max_nc x n_cells` matrix containing the inverse QR matrix at each cell (R\Q').

The output is a `n_fields x n_dims x n_cells_total` matrix `dUdX` containing the gradient of primitive variables in each cell.

The algorithm works by looping over the interior cells of the grid and:

1. Computing the difference `dU` between the solution vector `U` in the current cell and each of its encircling cells.
1. Using the inverse QR matrix `dXinv` to compute the gradient of `U` in the current cell as `dUdX = dXinv * dU`.
1. Assigning the computed gradient to the corresponding element of the output matrix `dUdX`.

### QR method

The QR method is a numerical technique for computing the gradient of a function using finite differences. It involves encircling a cell with a set of neighboring cells and using the differences between the function values at these cells to approximate the gradient at the center of the cell.

The method involves constructing a matrix Q whose columns are the difference vectors between the function values at the encircling cells and the center cell, and a matrix R whose columns are a set of orthonormal vectors that span the same space as the columns of Q. The gradient can then be computed as the product of the inverse of R and the transpose of Q (R\Q').

The QR method is an alternative to the traditional finite difference method, which involves directly computing the gradient as the difference between the function values at the encircling cells and the center cell, divided by the distance between the cells. The QR method has the advantage that it is more accurate and less sensitive to round-off errors, especially when the encircling cells are not equally spaced. However, it requires the computation of the QR decomposition, which can be more computationally expensive than the finite difference method.

## convective_flux

This function calculates the convective flux at a given interface between two fluid regions in a 2D computational fluid dynamics (CFD) simulation. The input variables are:

- u_l: a pointer to an array of 4 doubles representing the conservative variables of the left state at the interface
- u_r: a pointer to an array of 4 doubles representing the conservative variables of the right state at the interface
- norm: a pointer to an array of 2 doubles representing the outward-pointing normal vector of the interface
- area: a double representing the area of the interface
- fn: a pointer to an array of 4 doubles where the result will be stored

The conservative variables of the fluid are density, momentum, and energy.
These are stored in the arrays u_l and u_r in the following order: density, x-momentum, y-momentum, energy.

The function first calculates the velocities of the left and right states using the density and momentum.
It then calculates the pressure of the left and right states using the ideal gas law and the energy and velocity.
It then calculates the enthalpy of the left and right states using the density, pressure, and velocity.

Next, the function calculates the face-normal momentum of the left and right states by taking the dot product of the momentum and the normal vector.

Finally, the function calculates the Euler flux at the interface by averaging the left and right states using the face-normal momentum, velocities, and pressure.
The result is stored in the array pointed to by fn.
