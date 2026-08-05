from mpi4py import MPI

import numpy as np

import io4dolfinx

MPI.COMM_WORLD.barrier()

domain = io4dolfinx.read_mesh(
    "tally.vtkhdf",
    MPI.COMM_WORLD,
    backend="pyvista",
    backend_args={"dtype": np.float64},
)
f = io4dolfinx.read_cell_data(
    "tally.vtkhdf",
    "mean",
    domain,
    backend="pyvista",
)


names = io4dolfinx.read_function_names("tally.vtkhdf", domain.comm, backend="vtkhdf")
print(names)

import matplotlib as mpl
import dolfinx
import pyvista
import ufl
import numpy as np

from petsc4py import PETSc
from mpi4py import MPI

from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc,
)
from pathlib import Path

t = 0.0  # Start time
T = 10.0  # Final time
num_steps = 100
dt = T / num_steps  # time step size

V = dolfinx.fem.functionspace(domain, ("CG", 1))

u_n = fem.Function(V)
u_n.name = "u_n"
# u_n.interpolate(initial_condition)


# Create boundary condition
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
)
bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)


xdmf = io.XDMFFile(domain.comm, "diffusion.xdmf", "w")
xdmf.write_mesh(domain)

filename = Path("diffusion.vtkhdf")
backend_args = {"name": "MyGrid"}

io4dolfinx.write_mesh(
    filename, domain, backend="vtkhdf", mode=io4dolfinx.FileMode.write, backend_args=backend_args
)

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
    xdmf.write_function(uh, t)

    io4dolfinx.write_point_data(
        filename,
        uh,
        time=t,
        backend="vtkhdf",
        mode=io4dolfinx.FileMode.append,
        backend_args=backend_args,
    )

xdmf.close()

A.destroy()
b.destroy()
solver.destroy()
