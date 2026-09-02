# OSSFE-Hackathon-2026

[![Install with pip](https://github.com/OSSFE/hackathon-2026/actions/workflows/install-with-pip.yml/badge.svg)](https://github.com/OSSFE/hackathon-2026/actions/workflows/install-with-pip.yml)

Material for the hackathon at the end of OSSFE 2026

![OpenMC to FEniCS workflow](flowchart.svg)

This is not as easy as I would like. `installation-with-pip.sh` does the lot with pip in
one environment. `installation-with-conda.sh` needs pip on top for the parts conda-forge
does not carry, and cadquery and fenics-dolfinx cannot share a conda environment, so that
route needs two, see #8.

Run the scripts in order. The first one needs nuclear data, and the pip install brings
`openmc_data_downloader` with it, so fetching the two nuclides this example uses is one
command. If you already have a library, point `OPENMC_CROSS_SECTIONS` at its
`cross_sections.xml` and skip the download.

```bash
openmc_data_downloader -l FENDL-3.1d -i H1 Fe56 -d nuclear_data
export OPENMC_CROSS_SECTIONS=$PWD/nuclear_data/cross_sections.xml

python 1-make-cad-get-neutron-heating.py
python 2-solve-heat-equation.py
python 3-make-images.py
```

The CAD, a box with a sphere cut out of it and the sphere sitting in the hole, made with
CadQuery. The two are separate volumes with their own materials.

![cad](cad.png)

The same geometry as OpenMC sees it on the xz plane, with the materials labelled and the
sampled source positions marked, made with `openmc.Model.plot`

![materials xz](materials-xz.png)

The tetrahedral mesh that the tally is scored on, which fills both the box and the sphere,
with half of it cut away. The three images below this one are made with PyVista.

![mesh](mesh.png)

The tally, which scores fast neutron flux between 10 and 20 MeV and stands in for the
heating here, sliced through the middle. It is on a log scale because it covers about
eight orders of magnitude.

![tally](tally.png)

The steady state temperature from FEniCS with the tally as the heat source, sliced through
the middle. The box and the sphere share nodes on the face where they meet, so heat
crosses from one material into the other.

![temperature](temperature.png)
