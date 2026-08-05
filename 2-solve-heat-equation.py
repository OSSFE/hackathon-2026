import h5py
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from pathlib import Path

from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    create_vector,
    set_bc,
)
from dolfinx.io import vtkhdf

tally_file = "tally.vtkhdf"
comm = MPI.COMM_WORLD

# OpenMC writes the points as float32, reading as float64 matches the PETSc scalar type
domain = vtkhdf.read_mesh(comm=comm, filename=tally_file, dtype=np.float64)

tdim = domain.topology.dim
num_local_cells = domain.topology.index_map(tdim).size_local

# dolfinx renumbers and repartitions the cells while reading, so the tally values need
# gathering through original_cell_index before they line up with the local cells.
# VTKHDF is just HDF5, so the tally array itself is a plain dataset read.
original_cells = domain.topology.original_cell_index[:num_local_cells]

with h5py.File(tally_file, "r") as fh:
    order = np.argsort(original_cells)  # h5py wants increasing indices
    values = fh["VTKHDF/CellData/mean"][original_cells[order]]
    heating = np.empty_like(values)
    heating[order] = values

# the neutron heating is a single value per cell, so it belongs in DG0
Q = fem.functionspace(domain, ("DG", 0))
f = fem.Function(Q, name="heating")
f.x.array[:num_local_cells] = heating
f.x.scatter_forward()

t = 0.0  # Start time
T = 10.0  # Final time
num_steps = 100
dt = T / num_steps  # time step size

V = fem.functionspace(domain, ("CG", 1))

u_n = fem.Function(V)
u_n.name = "u_n"


# Create boundary condition
fdim = tdim - 1
boundary_facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
)
bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)


# write_point_data wants the values at the local mesh vertices, but a CG1 dof array is
# neither vertex ordered nor free of ghosts, so build the permutation between them once
geometry_map = domain.geometry.index_map()
num_local_vertices = geometry_map.size_local
permutation = np.zeros(num_local_vertices + geometry_map.num_ghosts, dtype=np.int32)
permutation[domain.geometry.dofmaps[0].reshape(-1)] = V.dofmap.list.reshape(-1)

xdmf = io.XDMFFile(domain.comm, "diffusion.xdmf", "w")
xdmf.write_mesh(domain)

filename = Path("diffusion.vtkhdf")
vtkhdf.write_mesh(filename=filename, mesh=domain)

uh = fem.Function(V)
uh.name = "uh"

u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
a = u * v * ufl.dx + dt * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = (u_n + dt * f) * v * ufl.dx

bilinear_form = fem.form(a)
linear_form = fem.form(L)


A = assemble_matrix(bilinear_form, bcs=[bc])
A.assemble()
b = create_vector(fem.extract_function_spaces(linear_form))


solver = PETSc.KSP().create(domain.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

for i in range(num_steps):
    t += dt
    if domain.comm.rank == 0:
        print(f"Time step {i + 1}/{num_steps}, Time: {t:.2f}")

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, linear_form)

    # Apply Dirichlet boundary condition to the vector
    apply_lifting(b, [bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, [bc])

    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()

    # Update solution at previous time step (u_n)
    u_n.x.array[:] = uh.x.array

    # Write solution to file
    xdmf.write_function(u=uh, t=t)
    vtkhdf.write_point_data(
        filename=filename,
        mesh=domain,
        data=uh.x.array[permutation][:num_local_vertices],
        time=t,
    )

xdmf.close()

A.destroy()
b.destroy()
solver.destroy()
