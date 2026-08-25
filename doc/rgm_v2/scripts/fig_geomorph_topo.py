#
# 3D topography renderings (pyvista) of the erosional surfaces of the
# four geomorphological unconformity types, from the depth maps produced
# by rgm/tmp/gen_paper_extras.f90 at the 3D mapping settings.
# The depth map (downward-positive) is rendered as an elevation surface
# so the channel/canyon/drainage morphology is directly visible.
#
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.font.family = 'arial'

ny, nx = 201, 251
cases = [('topo_meander_channel_201x251.bin', '(a) Meandering channels'),
         ('topo_meander_canyon_201x251.bin', '(b) Meandering canyon'),
         ('topo_drainage_channel_201x251.bin', '(c) Drainage channels'),
         ('topo_drainage_canyon_201x251.bin', '(d) Drainage canyon')]

p = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2200, 1500),
               border=False)

for i, (fname, title) in enumerate(cases):
    d = np.fromfile(f'../tmp/{fname}',
                    dtype=np.float32).reshape((ny, nx), order='F')
    # render as an elevation model: interfluves high, channels carved down
    zscale = 0.22*max(ny, nx)/max(d.max(), 1.0)
    elev = (d.max() - d)*zscale
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    grid = pv.StructuredGrid(x.astype(np.float32), y.astype(np.float32),
                             elev.astype(np.float32))
    grid.point_data['depth'] = d.flatten(order='F')

    p.subplot(i // 2, i % 2)
    p.set_background('white')
    p.add_mesh(grid, scalars='depth', cmap='gist_earth_r',
               smooth_shading=True, show_scalar_bar=False,
               specular=0.3, diffuse=0.9, ambient=0.25)
    p.add_text(title, font_size=14, color='black', position='upper_left')
    p.camera_position = 'iso'
    p.camera.elevation = 8
    p.camera.azimuth = 8
    p.camera.zoom(1.45)

p.screenshot('../fig_geomorph_topo.png')
print('saved fig_geomorph_topo.png')
