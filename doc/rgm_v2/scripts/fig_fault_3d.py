#
# 3D rendering (pyvista) of the fault surfaces of the models shown in the
# strike-varying and varying-displacement slice examples (the vstrike_1
# and vdisp_1 models of example/example_rgm3_curved.f90):
# (a) strike-varying, through-going faults colored by the local strike;
# (b) the same configuration with spatially varying displacement,
#     colored by the local displacement.
#
import numpy as np
import pyvista as pv
from scipy.ndimage import distance_transform_edt

pv.OFF_SCREEN = True
pv.global_theme.font.family = 'arial'

n1, n2, n3 = 171, 251, 251

def load(name):
    return np.fromfile(f'../tmp/example_3d_{name}_1.bin',
                       dtype=np.float32).reshape((n1, n2, n3), order='F')

def vol3(a):
    # axes: x = x3, y = x2, z = depth (down; camera view-up flipped)
    return np.transpose(a, (2, 1, 0)).copy()

def fault_surface(fid, attr):
    # propagate fault-voxel attribute values outward (nearest fault voxel)
    # so that isosurface interpolation picks up correct values
    mask = fid > 0.5
    _, idx = distance_transform_edt(~mask, return_indices=True)
    fill = attr[tuple(idx)]
    grid = pv.ImageData(dimensions=(n3, n2, n1))
    grid.point_data['fault'] = vol3(mask.astype(np.float32)).flatten(order='F')
    grid.point_data['attr'] = vol3(fill).flatten(order='F')
    surf = grid.contour([0.5], scalars='fault')
    return surf.smooth(n_iter=20)

surf_a = fault_surface(load('fault_vstrike'), load('fstrike_vstrike'))
surf_b = fault_surface(load('fault_vdisp'), load('fdisp_vdisp'))

# focal shifted along the screen-right direction so the scene sits left
# of the scalar bar
focal = (0.5*(n3 - 1) + 18, 0.5*(n2 - 1) - 14, 0.58*(n1 - 1))
pos = (focal[0] + 360, focal[1] + 290, focal[2] - 520)
# panel (b): viewpoint rotated about the vertical axis and raised so the
# bull's-eye displacement patterns face the camera
pos_b = (focal[0] + 60, focal[1] + 534, focal[2] - 330)

p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1250),
               border=False)

for i, (surf, title, btitle, clim) in enumerate(
        [(surf_a, '(a)', 'Strike (degree)', (0, 180)),
         (surf_b, '(b)', 'Displacement\n(Grid Number)', (0, 25))]):
    p.subplot(0, i)
    p.set_background('white')
    p.add_mesh(surf, scalars='attr', cmap='jet', smooth_shading=True,
               ambient=0.35, diffuse=0.85, specular=0.1, clim=clim,
               scalar_bar_args=dict(title=btitle + ' '*i, color='black',
                                    vertical=True, position_x=0.89,
                                    position_y=0.25, height=0.5,
                                    width=0.045, title_font_size=30,
                                    label_font_size=26, fmt='%.0f',
                                    n_labels=5 if i == 0 else 6))
    p.add_text(title, font_size=20, color='black', position='upper_left')
    p.show_bounds(bounds=(0, n3 - 1, 0, n2 - 1, 0, n1 - 1),
                  xtitle='X (Grid Number)', ytitle='Y (Grid Number)',
                  ztitle='Z (Grid Number)', color='black', font_size=22,
                  location='outer', ticks='outside', grid=False,
                  n_xlabels=6, n_ylabels=5, n_zlabels=4)
    p.camera_position = [pos if i == 0 else pos_b, focal, (0, 0, -1)]
    p.camera.zoom(0.80)

p.screenshot('../fig_fault_3d.png')
print('saved fig_fault_3d.png')
