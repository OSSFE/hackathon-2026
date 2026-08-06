

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
source ~/miniforge3/bin/activate

conda create --name demo -y
conda activate demo
# fenics-dolfinx 0.11 is needed for the native VTKHDF reader (dolfinx.io.vtkhdf) and
# openmc takes its dagmc_nompi build because nothing here runs in parallel.
# NOTE: this does not currently resolve. cadquery and fenics-dolfinx cannot be installed
# into the same conda environment, so for now the two scripts need separate environments.
conda install -c conda-forge \
  "fenics-dolfinx>=0.11" \
  "cadquery>=2.8.0" \
  "cad_to_dagmc>=0.13.2" \
  "openmc=*=dagmc_nompi*" \
  pyvista dagmc_h5m_file_inspector -y

# conda-forge only has cad_to_dagmc 0.13.2 and does not carry cad-to-dagmc-mesher at all,
# so the newer version and the mesher come from pip. --no-deps stops pip replacing the
# conda builds of cadquery and the rest that they depend on.
pip install --no-deps \
  "cad_to_dagmc>=0.14.0" \
  "cad-to-dagmc-mesher>=0.2.0" \
  "cadquery_direct_mesh_plugin>=0.2.0"
