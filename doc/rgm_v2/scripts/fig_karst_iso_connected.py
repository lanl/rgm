#
# 3D isosurface rendering (pyvista) of the high-connectivity karst
# example (karst_connect = 0.95, karst_npassage = 20), produced by
# rgm/tmp/gen_paper_extras.f90, with a labeled axis box (depth
# increasing downward).
#
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.font.family = 'arial'

n1, n2, n3 = 151, 201, 251
vol = np.fromfile('../tmp/karst_connected_mask_151x201x251.bin',
                  dtype=np.float32).reshape((n1, n2, n3), order='F')

v = np.transpose(vol, (2, 1, 0)).copy()

grid = pv.ImageData(dimensions=v.shape)
grid.point_data['karst'] = v.flatten(order='F')
grid.point_data['depth'] = np.broadcast_to(
    np.arange(v.shape[2])[None, None, :], v.shape
).flatten(order='F').astype(np.float32)

surf = grid.contour([0.5], scalars='karst')
surf = surf.smooth(n_iter=30)

p = pv.Plotter(off_screen=True, window_size=(2000, 1350))
p.set_background('white')
p.add_mesh(surf, scalars='depth', cmap='jet', smooth_shading=True,
           scalar_bar_args=dict(title='Z (Grid Number)', color='black',
                                vertical=True, position_x=0.90,
                                position_y=0.25, height=0.5, width=0.05,
                                title_font_size=26, label_font_size=22))
p.show_bounds(bounds=(0, n3 - 1, 0, n2 - 1, 0, n1 - 1),
              xtitle='X (Grid Number)', ytitle='Y (Grid Number)',
              ztitle='Z (Grid Number)', color='black', font_size=18,
              location='outer', ticks='outside', grid='back', n_xlabels=6, n_ylabels=5, n_zlabels=4)

focal = (0.5*(n3 - 1), 0.5*(n2 - 1), 0.62*(n1 - 1))
pos = (focal[0] + 385, focal[1] + 305, focal[2] - 255)
p.camera_position = [pos, focal, (0, 0, -1)]
# fit the view to the scene while keeping the viewing direction
p.reset_camera()
p.camera.zoom(1.12)
p.screenshot('../fig_karst_iso_connected.png')
print('saved fig_karst_iso_connected.png')
