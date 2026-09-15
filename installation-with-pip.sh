#!/bin/bash
set -e

# Everything goes into one virtual environment, all from wheels. openmc, dolfinx, basix
# and petsc all come from the same custom index, the rest from PyPI. Nothing is compiled,
# so the apt list is only what cadquery and pyvista need to open a window.

sudo apt update
sudo apt install -y \
  python3-venv \
  libgl1 libglu1-mesa libxrender1 libxext6 \
  libxft2 libxinerama1 libxcursor1 libxfixes3 libfontconfig1

python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip

# openmc wheel, which bundles dagmc and moab so no separate install is needed
pip install --extra-index-url https://shimwell.github.io/wheels openmc

# the cad and neutronics side, these all have wheels. openmc-data-downloader is here so
# that getting the nuclear data for the first script is one command, see the README.
pip install cadquery cad_to_dagmc dagmc_h5m_file_inspector numpy h5py pyvista openmc-data-downloader

# mpich supplies the MPI runtime. The petsc wheels are built against it, so it has to be
# the MPI in use. openmc is a serial build and links no MPI, so nothing clashes.
pip install mpich

# dolfinx, basix, ffcx, ufl, petsc and petsc4py, all as wheels. --pre is required because
# these are development versions, and without it pip ignores them and falls back to
# building petsc from the source distribution on PyPI, which is the slow path this
# replaces. petsc4py comes in through the fenics-dolfinx extra.
pip install --pre --extra-index-url https://shimwell.github.io/wheels \
  "fenics-dolfinx[petsc4py]" petsc

python -c "import openmc, cadquery, cad_to_dagmc, dolfinx; print('dolfinx', dolfinx.__version__)"
