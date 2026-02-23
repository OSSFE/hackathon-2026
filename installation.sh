

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
source ~/miniforge3/bin/activate

conda create --name demo -y
conda activate demo
conda install conda-forge::fenics-dolfinx cadquery=2.7.0 cad_to_dagmc==0.11.2 openmc=0.15.3  mpich pyvista dagmc_h5m_file_inspector -y
