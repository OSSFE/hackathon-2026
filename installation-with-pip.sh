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

# the cad and neutronics side, these all have wheels
pip install cadquery "cad_to_dagmc>=0.14.1" dagmc_h5m_file_inspector numpy h5py pyvista

# petsc is built from source into the venv, this is the slow part. dolfinx 0.11 wants
# 3.25 or newer, which is why apt petsc is not used. It comes from the release branch
# rather than the PyPI wheel because petsc merge request 9488 keeps the pkg-config files
# in the wheel, and that is how dolfinx finds petsc. Building petsc4py from the same
# clone is also what keeps the two versions in step.
pip install mpi4py
git clone --branch release --depth 1 https://gitlab.com/petsc/petsc.git
pip install ./petsc
# petsc4py has to be built without isolation, otherwise it records the temporary pip
# build directory as PETSC_DIR and dolfinx cannot find libpetsc.so at run time
pip install setuptools wheel cython
pip install --no-build-isolation ./petsc/src/binding/petsc4py

# petsc ships PETSc.pc inside its own install now, it just has to be on the search path
petsc_dir="$(python -c 'import petsc; print(petsc.get_petsc_dir())')"
export PKG_CONFIG_PATH="$petsc_dir/lib/pkgconfig:$PKG_CONFIG_PATH"

# the petsc wheel ships libpetsc.so as a linker script rather than a library, which the
# loader dolfinx uses cannot read, so make it a symlink to the real thing instead
ln -sf "$(cd "$petsc_dir/lib" && ls libpetsc.so.* | head -1)" "$petsc_dir/lib/libpetsc.so"

# the pure python parts of fenics, plus the tools needed to build dolfinx
pip install fenics-ufl fenics-ffcx scikit-build-core nanobind cffi

export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"

# basix, the c++ library then the python bindings
git clone --branch v0.11.0 --depth 1 https://github.com/FEniCS/basix.git
cmake -G Ninja -B basix/build -S basix/cpp -DCMAKE_INSTALL_PREFIX="$VIRTUAL_ENV"
cmake --build basix/build
cmake --install basix/build
pip install --no-build-isolation ./basix/python

# dolfinx 0.11 is needed for the native VTKHDF reader (dolfinx.io.vtkhdf)
git clone --branch v0.11.0 --depth 1 https://github.com/FEniCS/dolfinx.git
cmake -G Ninja -B dolfinx/build -S dolfinx/cpp -DCMAKE_INSTALL_PREFIX="$VIRTUAL_ENV"
cmake --build dolfinx/build
cmake --install dolfinx/build
pip install --no-build-isolation ./dolfinx/python

python -c "import openmc, cadquery, cad_to_dagmc, dolfinx; print('dolfinx', dolfinx.__version__)"
