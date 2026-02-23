import openmc
from pathlib import Path
from cad_to_dagmc import CadToDagmc
import cadquery as cq

# from cadquery.vis import show
import dagmc_h5m_file_inspector as di

assembly = cq.Assembly()
box = cq.Workplane("XY").box(30, 30, 30)
sphere = cq.Workplane("XY").moveTo(20, 0).sphere(10)
assembly.add(box.cut(sphere), name="box")
assembly.add(sphere, name="sphere")
# show(assembly)

# assembly = cq.Assembly()
# box = cq.Workplane("XY").box(30, 30, 30)
# sphere = cq.Workplane("XY").sphere(10)
# assembly.add(sphere, name="sphere")
# assembly.add(box.cut(sphere), name="box")

# my_model = CadToDagmc()
# my_model.add_cadquery_object(
#     cadquery_object=assembly,
#     material_tags="assembly_names"
# )

# my_model.export_dagmc_h5m_file(
#     filename="dagmc.h5m",
#     set_size={"sphere": 0.5,"box": 15,},
#     unstructured_volumes=["sphere"],
#     umesh_filename="umesh.vtk",
#     meshing_backend="gmsh",
# )

# di.convert_h5m_to_vtkhdf('dagmc.h5m', 'dagmc.vtkhdf')
# # Setting the cross section path to the correct location in the docker image.
# # If you are running this outside the docker image you will have to change this path to your local cross section path.
# openmc.config['cross_sections'] = Path.home() / 'nuclear_data' / 'endf-b8.0-hdf5/' / 'cross_sections.xml'

# mat1 = openmc.Material(name='box')
# mat1.add_nuclide("H1", 1, percent_type="ao")
# mat1.set_density("g/cm3", 0.001)

# mat2 = openmc.Material(name='sphere')
# mat2.add_nuclide("Fe56", 1, percent_type="ao")
# mat2.set_density("g/cm3", 7)
# my_materials = openmc.Materials([mat1, mat2])


# dag_univ = openmc.DAGMCUniverse(filename="dagmc.h5m").bounded_universe(padding_distance=10)
# my_geometry = openmc.Geometry(root=dag_univ)


# umesh = openmc.UnstructuredMesh("umesh.vtk", library="moab")
# mesh_filter = openmc.MeshFilter(umesh)
# energy_filter = openmc.EnergyFilter([10.0e6, 20.0e6])  # filtering neutrons between 10MeV and 20MeV, this is just an example and your definition of fast neutron may vary

# tally = openmc.Tally(name="unstructured_mesh_tally")
# tally.filters = [mesh_filter, energy_filter]
# tally.scores = ["flux"]
# my_tallies = openmc.Tallies([tally])


# my_settings = openmc.Settings()
# my_settings.batches = 10
# my_settings.particles = 5000
# my_settings.run_mode = "fixed source"


# my_source = openmc.IndependentSource()
# source_location = my_geometry.bounding_box.center
# source_location[2] =+ 9  # moving the source above the geometry so more tehraherals get hit
# my_source.space = openmc.stats.Point(source_location)
# my_source.angle = openmc.stats.Isotropic()
# my_source.energy = openmc.stats.Discrete([14e6], [1])
# my_settings.source = my_source

# model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)
# sp_filename = model.run()


# sp = openmc.StatePoint(sp_filename)

# tally_result = sp.get_tally(name="unstructured_mesh_tally")

# flux_mean = tally_result.get_values(scores=["flux"], value="mean").flatten()


# umesh_from_sp = tally_result.find_filter(openmc.MeshFilter).mesh


# umesh_from_sp.write_data_to_vtk(filename="tally.vtkhdf", datasets={"mean": flux_mean})
