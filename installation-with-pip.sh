#!/bin/bash
set -e

# Everything here goes into one virtual environment. openmc comes from a wheel, the cad
# and dagmc packages come from PyPI, and fenics-dolfinx is built from source because it
# has no wheel on PyPI. The apt packages are the c++ libraries that build needs, and they
# are all the openmpi flavour so that petsc, hdf5 and mpi4py agree with each other.

sudo apt update
sudo apt install -y \
  build-essential cmake ninja-build pkg-config git \
  python3-dev python3-venv \
  libopenmpi-dev libhdf5-openmpi-dev libopenblas-dev \
  libboost-dev libpugixml-dev libspdlog-dev libptscotch-dev \
  libgl1 libglu1-mesa libxrender1 libxext6 \
  libxft2 libxinerama1 libxcursor1 libxfixes3 libfontconfig1

python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip

# openmc wheel, which bundles dagmc and moab so no separate install is needed
pip install --extra-index-url https://shimwell.github.io/wheels openmc

# the cad and neutronics side, these all have wheels. openmc-data-downloader is here so
# that getting the nuclear data for the first script is one command, see the README
pip install cadquery cad_to_dagmc dagmc_h5m_file_inspector numpy h5py pyvista \
  openmc-data-downloader cad-to-dagmc-mesher

# petsc is compiled while pip installs the sdist, and this is the slow part. dolfinx 0.11
# wants 3.25 or newer, which is why apt petsc is not used. 3.25.5 is the first release on
# PyPI to keep the pkg-config files in the wheel (petsc merge request 9488), and that is
# how dolfinx finds petsc, so the git clone of the release branch is no longer needed.
# petsc4py is pinned to the same version to keep the two in step.
pip install mpi4py
pip install "petsc==3.25.5"

# petsc lands in site-packages, and both petsc4py and dolfinx find it from here. dolfinx
# prepends $PETSC_DIR/lib/pkgconfig to the pkg-config search path itself, which is exactly
# where the wheel now keeps PETSc.pc, so no PKG_CONFIG_PATH juggling is needed
export PETSC_DIR="$(python -c 'import petsc; print(petsc.get_petsc_dir())')"

# petsc4py takes PETSC_DIR from the environment (setup.cfg has petsc_dir = $PETSC_DIR) and
# records it, so the export above has to happen first, otherwise dolfinx cannot find
# libpetsc.so at run time
pip install "petsc4py==3.25.5"

# the petsc wheel ships libpetsc.so as a linker script rather than a library, which the
# loader dolfinx uses cannot read, so make it a symlink to the real thing instead. dolfinx
# pull request 4468 fixes this on the dolfinx side, but it is only on main, so this stays
# for as long as we build 0.11
ln -sf "$(cd "$PETSC_DIR/lib" && ls libpetsc.so.* | head -1)" "$PETSC_DIR/lib/libpetsc.so"

# the pure python parts of fenics, plus the tools needed to build dolfinx. nanobind is
# pinned because it only shares its type registry between extensions built with the same
# internal ABI version, and the fenics-basix wheel below is prebuilt. 2.12.0 was the
# newest nanobind when that wheel was published. Building dolfinx against 3.x instead
# leaves basix elements unrecognisable to dolfinx, which surfaces at run time as
# "incompatible function arguments" from fem.functionspace rather than at build time
pip install fenics-ufl fenics-ffcx scikit-build-core "nanobind==2.12.0" cffi

export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"

# the basix wheel carries libbasix.so, the headers and BasixConfig.cmake inside the python
# package, and dolfinx asks the interpreter where basix lives (DOLFINX_BASIX_PYTHON, on by
# default), so there is no c++ library to build here. The version has to track dolfinx.
pip install "fenics-basix==0.11.0"

# dolfinx 0.11 is needed for the native VTKHDF reader (dolfinx.io.vtkhdf)
git clone --branch v0.11.0 --depth 1 https://github.com/FEniCS/dolfinx.git
cmake -G Ninja -B dolfinx/build -S dolfinx/cpp -DCMAKE_INSTALL_PREFIX="$VIRTUAL_ENV"
cmake --build dolfinx/build
cmake --install dolfinx/build
pip install --no-build-isolation ./dolfinx/python

python -c "import openmc, cadquery, cad_to_dagmc, dolfinx; print('dolfinx', dolfinx.__version__)"
