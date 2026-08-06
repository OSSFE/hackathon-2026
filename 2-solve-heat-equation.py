import h5py
import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import vtkhdf

heating_file = "heating.vtkhdf"

# OpenMC writes the points as float32, reading as float64 matches the PETSc scalar type
domain = vtkhdf.read_mesh(comm=MPI.COMM_WORLD, filename=heating_file, dtype=np.float64)

tdim = domain.topology.dim
num_local_cells = domain.topology.index_map(tdim).size_local

# dolfinx renumbers and repartitions the cells while reading, so the tally values need
# gathering through original_cell_index before they line up with the local cells.
# VTKHDF is just HDF5, so the tally array itself is a plain dataset read.
original_cells = domain.topology.original_cell_index[:num_local_cells]

with h5py.File(heating_file, "r") as fh:
    order = np.argsort(original_cells)  # h5py wants increasing indices
    values = fh["VTKHDF/CellData/mean"][original_cells[order]]
    heating = np.empty_like(values)
    heating[order] = values

# the neutron heating is a single value per cell, so it belongs in DG0
Q = fem.functionspace(domain, ("DG", 0))
f = fem.Function(Q, name="heating")
f.x.array[:num_local_cells] = heating
f.x.scatter_forward()

V = fem.functionspace(domain, ("CG", 1))

# hold the outside of the model at zero temperature
boundary_facets = mesh.locate_entities_boundary(
    domain, tdim - 1, lambda x: np.full(x.shape[1], True, dtype=bool)
)
bc = fem.dirichletbc(0.0, fem.locate_dofs_topological(V, tdim - 1, boundary_facets), V)

u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options_prefix="heat_",
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
)
temperature = problem.solve()

# write_point_data wants the values at the local mesh vertices, but a CG1 dof array is
# neither vertex ordered nor free of ghosts, so build the permutation between them
geometry_map = domain.geometry.index_map()
num_local_vertices = geometry_map.size_local
permutation = np.zeros(num_local_vertices + geometry_map.num_ghosts, dtype=np.int32)
permutation[domain.geometry.dofmaps[0].reshape(-1)] = V.dofmap.list.reshape(-1)

output_file = "temperature.vtkhdf"
vtkhdf.write_mesh(filename=output_file, mesh=domain)
vtkhdf.write_point_data(
    filename=output_file,
    mesh=domain,
    data=temperature.x.array[permutation][:num_local_vertices],
    time=0.0,
)
