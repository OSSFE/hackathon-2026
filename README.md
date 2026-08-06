# OSSFE-Hackathon-2026

Material for the hackathon at the end of OSSFE 2026

![OpenMC to FEniCS workflow](flowchart.svg)

This is not as easy as I would like but you can install partly with conda and partly with pip

The CAD, a box with a sphere cut out of it, made with CadQuery

![cad](cad.png)

The same geometry as OpenMC sees it, with the materials labelled and the source position
marked, made with `openmc.Model.plot`

![materials xz](materials-xz.png)

The three images below are slices through the middle of the model, made with PyVista.

The tetrahedral mesh that the tally is scored on, which fills both the box and the sphere

![mesh](mesh.png)

The OpenMC tally, on a log scale because the flux covers several orders of magnitude

![tally](tally.png)

The steady state temperature from FEniCS, with the tally as the heat source

![temperature](temperature.png)


