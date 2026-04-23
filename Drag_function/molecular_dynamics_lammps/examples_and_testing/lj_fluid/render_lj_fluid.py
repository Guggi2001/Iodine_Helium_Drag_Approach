from ovito.io import import_file
from ovito.modifiers import ColorByTypeModifier
from ovito.vis import Viewport, TachyonRenderer

# 1. Load the LAMMPS trajectory
# Ensure this matches the filename you defined in your dump command
pipeline = import_file("dump.lj_restart")

# 2. Assign different colors to Type 1 and Type 2 atoms to visualize demixing
pipeline.modifiers.append(ColorByTypeModifier())
pipeline.add_to_scene()

# 3. Configure the viewport and camera
vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.zoom_all(size=(800, 600))

# 4. Render the animation
output_video = "lj_fluid_restart.mp4"
print(f"Starting render of {pipeline.source.num_frames} frames...")

vp.render_anim(
    filename=output_video,
    size=(800, 600),
    fps=15,
    background=(1.0, 1.0, 1.0),
    renderer=TachyonRenderer(shadows=True, direct_light_intensity=1.2)
)

print(f"Render complete: {output_video}")