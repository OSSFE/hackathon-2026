# Makes the images in the README. This needs ParaView's own python, so run it with
#   pvpython 3-make-paraview-images.py
# and on a machine with no display use
#   xvfb-run -a pvpython 3-make-paraview-images.py

from paraview.simple import *

size = [1400, 1000]

# --------------------------------------------------------------------- the mesh
tally = OpenDataFile("tally.vtkhdf")
UpdatePipeline(proxy=tally)

view = CreateRenderView()
view.ViewSize = size
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0

# slice through the middle, the outside of the box is only a few large triangles
mesh_slice = Slice(Input=tally)
mesh_slice.SliceType = "Plane"
mesh_slice.SliceType.Origin = [0.0, 0.0, 0.0]
mesh_slice.SliceType.Normal = [0.0, 1.0, 0.0]

display = Show(mesh_slice, view)
display.SetRepresentationType("Surface With Edges")
display.ColorArrayName = ["POINTS", ""]
display.DiffuseColor = [0.62, 0.75, 0.90]
display.EdgeColor = [0.16, 0.20, 0.26]

ResetCamera(view)
GetActiveCamera().Elevation(-90)  # look straight down onto the slice
ResetCamera(view)
SaveScreenshot("mesh.png", view, ImageResolution=size)

# -------------------------------------------------------------------- the tally
view = CreateRenderView()
view.ViewSize = size
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0

tally_slice = Slice(Input=tally)
tally_slice.SliceType = "Plane"
tally_slice.SliceType.Origin = [0.0, 0.0, 0.0]
tally_slice.SliceType.Normal = [0.0, 1.0, 0.0]
UpdatePipeline(proxy=tally_slice)

display = Show(tally_slice, view)
display.SetRepresentationType("Surface")
ColorBy(display, ("CELLS", "mean"))

# the flux covers several orders of magnitude, so colour it on a log scale
highest = tally_slice.CellData["mean"].GetRange()[1]
colours = GetColorTransferFunction("mean")
colours.ApplyPreset("Viridis", True)
colours.RescaleTransferFunction(highest * 1e-4, highest)
colours.MapControlPointsToLogSpace()
colours.UseLogScale = 1
display.SetScalarBarVisibility(view, True)

ResetCamera(view)
GetActiveCamera().Elevation(-90)
ResetCamera(view)
SaveScreenshot("tally.png", view, ImageResolution=size)

# -------------------------------------------------------------- the temperature
temperature = OpenDataFile("temperature.vtkhdf")
UpdatePipeline(proxy=temperature)

view = CreateRenderView()
view.ViewSize = size
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0

temperature_slice = Slice(Input=temperature)
temperature_slice.SliceType = "Plane"
temperature_slice.SliceType.Origin = [0.0, 0.0, 0.0]
temperature_slice.SliceType.Normal = [0.0, 1.0, 0.0]

display = Show(temperature_slice, view)
display.SetRepresentationType("Surface")
ColorBy(display, ("POINTS", "u"))
display.RescaleTransferFunctionToDataRange(True, False)
GetColorTransferFunction("u").ApplyPreset("Inferno", True)
display.SetScalarBarVisibility(view, True)

ResetCamera(view)
GetActiveCamera().Elevation(-90)
ResetCamera(view)
SaveScreenshot("temperature.png", view, ImageResolution=size)
