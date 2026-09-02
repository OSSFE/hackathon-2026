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
# that getting the nuclear data for the first script is one command, see the README.
# TODO: drop the cad-to-dagmc-mesher bound once fusion-energy/cad-to-dagmc-mesher#157 is
# fixed. 0.3.0 panics on this geometry with "orient_3d_sos: all projections degenerate"
pip install cadquery cad_to_dagmc dagmc_h5m_file_inspector numpy h5py pyvista \
  openmc-data-downloader "cad-to-dagmc-mesher<0.3.0"

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

# the pure python parts of fenics, plus the tools needed to build dolfinx. basix and
# dolfinx are both built from source below so they agree with each other whatever nanobind
# is installed, but it is pinned anyway because the 0.11 releases predate nanobind 3 and
# were never built against it
pip install fenics-ufl fenics-ffcx scikit-build-core "nanobind==2.12.0" cffi

export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"

# basix, the c++ library then the python bindings. The fenics-basix wheel cannot be used
# here even though it carries libbasix.so and BasixConfig.cmake: it is a Py_LIMITED_API
# build, and nanobind puts that in its ABI tag (NB_STABLE_ABI in nb_abi.h), so its types
# are invisible to the dolfinx bindings built below. That shows up at run time as
# "incompatible function arguments" from fem.functionspace. Building both from source
# keeps them on one ABI tag. The version has to track dolfinx.
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
