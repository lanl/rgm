#
# 3D rendering of the spatially varying fault displacement:
# (a) the strike-varying, listric fault surface v = g(u) + dev(z),
#     colored by the elliptical slip patch alpha(u, z);
# (b) volume rendering of the displacement magnitude
#     delta = alpha(u, z) * beta(v_perp) in the hanging-wall block,
#     showing the Gaussian decay away from the fault surface.
#
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True

# ---- fault geometry (same construction as the 2D schematics) ----
nu, nzz = 221, 121
u = np.linspace(-1.0, 1.0, nu)
z = np.linspace(0.0, 1.0, nzz)

dphi_max = np.deg2rad(20.0)
dphi = dphi_max*(0.8*np.sin(2.1*np.pi*0.5*(u + 0.35))
                 + 0.4*np.sin(2.0*np.pi*1.15*(u - 0.2)))
dphi = dphi/np.max(np.abs(dphi))*dphi_max
g = np.cumsum(np.tan(dphi))*(u[1] - u[0])
g -= np.interp(0.0, u, g)

theta = np.deg2rad(np.linspace(72.0, 50.0, nzz))
dev = np.cumsum(1.0/np.tan(theta))*(z[1] - z[0])

# elliptical slip patch
uc, zc, au, az = 0.05, 0.45, 0.75, 0.42
Ug, Zg = np.meshgrid(u, z, indexing='ij')          # (nu, nzz)
alpha = np.maximum(0.0, 1.0 - ((Ug - uc)/au)**2 - ((Zg - zc)/az)**2)
Vg = g[:, None] + dev[None, :]

# depth positive down; camera view-up flipped so shallow is up on screen
surf = pv.StructuredGrid(Ug[:, :, None], Vg[:, :, None], Zg[:, :, None])
surf.point_data['disp'] = alpha.flatten(order='F')

sbar = dict(title='Disp / max disp', color='black', vertical=True,
            position_x=0.84, position_y=0.25, height=0.5, width=0.045,
            title_font_size=30, label_font_size=26)

# focal shifted along the screen-right direction so the scene sits left
# of the scalar bar
focal = (0.30, 0.5*(Vg.min() + Vg.max()) + 0.06, 0.5)
pos = (focal[0] + 2.0, focal[1] + 2.5, focal[2] - 1.4)

p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(2400, 1050),
               border=False)

# ------------------------------------------- (a) fault surface colored by alpha
p.subplot(0, 0)
p.set_background('white')
p.add_mesh(surf, scalars='disp', cmap='jet', clim=[0.0, 1.0],
           smooth_shading=True, ambient=0.45, diffuse=0.85, specular=0.1,
           scalar_bar_args=dict(sbar))
p.add_text('(a)',
           font_size=18, color='black', position='upper_left')
p.show_bounds(xtitle='u (along base strike)', ytitle='v (strike normal)',
              ztitle='depth z', color='black', font_size=20,
              location='outer', ticks='outside', grid=False,
              n_xlabels=5, n_ylabels=4, n_zlabels=3)
p.camera_position = [pos, focal, (0, 0, -1)]
p.camera.zoom(0.95)

# --------------------------------- (b) volume rendering of the decaying field
nvv = 141
vmin, vmax = float(Vg.min()) - 0.15, float(Vg.max()) + 0.55
v3 = np.linspace(vmin, vmax, nvv)
w = 0.28

# delta(u, v, z) = alpha(u, z) * exp(-vperp^2 / 2w^2) in the hanging wall
Vperp = v3[None, :, None] - Vg[:, None, :]         # (nu, nvv, nzz)
beta = np.exp(-0.5*(Vperp/w)**2)
beta[Vperp < 0] = 0.0                              # footwall stays fixed
delta = alpha[:, None, :]*beta

vol = pv.ImageData(dimensions=(nu, nvv, nzz),
                   spacing=(u[1] - u[0], v3[1] - v3[0], z[1] - z[0]),
                   origin=(u[0], v3[0], 0.0))
vol.point_data['disp'] = delta.flatten(order='F')

p.subplot(0, 1)
p.set_background('white')
sbar_b = dict(sbar)
sbar_b['title'] = 'Disp / max disp '     # distinct name: separate lookup table

# solid cut-block rendering: remove the volume above the horizontal cut
# plane through the slip-patch center, so the cut face shows the map-view
# Gaussian decay of the displacement away from the curved fault trace
zcut = zc
block = vol.clip(normal=(0, 0, 1), origin=(0.0, 0.0, zcut), invert=False)
p.add_mesh(block, scalars='disp', cmap='jet', clim=[0.0, 1.0],
           ambient=0.4, diffuse=0.85, scalar_bar_args=sbar_b)
# fault surface above the cut, for context
surf_up = surf.clip(normal=(0, 0, 1), origin=(0.0, 0.0, zcut), invert=True)
p.add_mesh(surf_up, color='lightgray', opacity=0.55, smooth_shading=True,
           ambient=0.45)
p.add_text('(b)',
           font_size=18, color='black', position='upper_left')
p.show_bounds(xtitle='u (along base strike)', ytitle='v (strike normal)',
              ztitle='depth z', color='black', font_size=20,
              location='outer', ticks='outside', grid=False,
              n_xlabels=5, n_ylabels=4, n_zlabels=3)
p.camera_position = [pos, focal, (0, 0, -1)]
p.camera.zoom(0.88)

p.screenshot('../fig_disp_3d.png')
print('saved fig_disp_3d.png')
