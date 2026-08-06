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
  libgl1 libglu1-mesa libxrender1 libxext6

python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip

# openmc wheel, which bundles dagmc and moab so no separate install is needed
pip install --extra-index-url https://shimwell.github.io/wheels openmc

# the cad and neutronics side, these all have wheels
pip install cadquery cad_to_dagmc dagmc_h5m_file_inspector numpy h5py

# petsc is built from source into the venv, this is the slow part. The version has to
# match petsc4py, and dolfinx 0.11 wants 3.25 or newer, which is why apt petsc is not used.
pip install mpi4py petsc
# petsc4py has to be built without isolation, otherwise it records the temporary pip
# build directory as PETSC_DIR and dolfinx cannot find libpetsc.so at run time
pip install setuptools wheel cython
pip install --no-build-isolation petsc4py

# dolfinx looks for petsc with pkg-config but the petsc wheel ships no .pc file
mkdir -p "$VIRTUAL_ENV/lib/pkgconfig"
cat > "$VIRTUAL_ENV/lib/pkgconfig/PETSc.pc" <<EOF
prefix=$(python -c "import petsc; print(petsc.get_petsc_dir())")
Name: PETSc
Description: Portable Extensible Toolkit for Scientific Computation
Version: $(python -c "import petsc4py; print(petsc4py.__version__)")
Cflags: -I\${prefix}/include
Libs: -L\${prefix}/lib -lpetsc
EOF
export PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"

# the petsc wheel ships libpetsc.so as a linker script rather than a library, which the
# loader dolfinx uses cannot read, so make it a symlink to the real thing instead
petsc_lib="$(python -c 'import petsc; print(petsc.get_petsc_dir())')/lib"
ln -sf "$(cd "$petsc_lib" && ls libpetsc.so.* | head -1)" "$petsc_lib/libpetsc.so"

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

python -c "import openmc, cadquery, dolfinx; print('dolfinx', dolfinx.__version__)"
