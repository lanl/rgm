#
# Schematic of the strike-varying fault construction:
# (a) map view: fault-local frame (u, v), base strike, curved trace
#     v = g(u) with local strike phi + dphi(u);
# (b) 3D fault surface v = g(u) + dev(z): the curved trace extruded
#     down-depth with the listric deviation (generalized cylinder).
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D  # noqa

# Arial for both text and math, with TrueType embedding for the PDF
mpl.rcParams.update({'font.size': 15,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans'],
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial:italic',
                     'mathtext.bf': 'Arial:bold',
                     'mathtext.cal': 'Arial',
                     'mathtext.default': 'it',
                     'pdf.fonttype': 42,
                     'ps.fonttype': 42})

# smooth low-wavenumber strike deviation and trace
nu = 401
u = np.linspace(-1, 1, nu)
dphi_max = np.deg2rad(20.0)
dphi = dphi_max*(0.8*np.sin(2.1*np.pi*0.5*(u + 0.35))
                 + 0.4*np.sin(2.0*np.pi*1.15*(u - 0.2)))
dphi = dphi/np.max(np.abs(dphi))*dphi_max
g = np.cumsum(np.tan(dphi))*(u[1] - u[0])
g -= np.interp(0.0, u, g)

phi = np.deg2rad(35.0)          # base strike in map view
R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

fig = plt.figure(figsize=(13, 5.2))

# ---------------------------------------------------------------- (a) map view
ax = fig.add_subplot(1, 2, 1)
trace = R @ np.vstack([u, g])
base = R @ np.vstack([u, np.zeros_like(u)])

ax.plot(base[0], base[1], '--', color='gray', lw=1.2)
ax.plot(trace[0], trace[1], color='crimson', lw=2.5)

# u and v axes of the fault-local frame
for vec, name, dx, dy, va in [((1.12, 0.0), '$u$ (base strike, $\\phi$)', 0.04, -0.09, 'top'),
                              ((0.0, 0.55), '$v$', 0.02, 0.03, 'baseline')]:
    tip = R @ np.array(vec)
    ax.annotate('', xy=tip, xytext=(0, 0),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.4))
    ax.text(tip[0] + dx, tip[1] + dy, name, fontsize=16, va=va)

# local tangent at u0
i0 = np.argmin(np.abs(u - 0.42))
t0 = np.array([1.0, np.tan(dphi[i0])])
t0 /= np.linalg.norm(t0)
p0 = R @ np.array([u[i0], g[i0]])
tt = R @ (0.42*t0)
ax.annotate('', xy=p0 + tt, xytext=p0,
            arrowprops=dict(arrowstyle='-|>', color='navy', lw=1.6))
bb = R @ (0.42*np.array([1.0, 0.0]))
ax.plot([p0[0], p0[0] + bb[0]], [p0[1], p0[1] + bb[1]], ':', color='navy', lw=1.1)
ax.text(p0[0] + tt[0] - 0.06, p0[1] + tt[1] + 0.16,
        '$\\phi_{\\mathrm{loc}} = \\phi + \\Delta\\phi(u)$',
        color='navy', fontsize=16, ha='right')

ax.plot(0, 0, 'ko', ms=5)
ax.text(0.03, -0.10, 'fault center', fontsize=14)
ax.text(*(R @ np.array([0.72, -0.16])), '$v = g(u)$', color='crimson', fontsize=16)
ax.text(*(R @ np.array([-0.95, 0.05])), 'base strike', color='gray', fontsize=14,
        rotation=np.rad2deg(phi), rotation_mode='anchor')
ax.set_xlabel('$x_3$ (normalized)')
ax.set_ylabel('$x_2$ (normalized)')
ax.set_title('(a)')
ax.set_aspect('equal')
ax.set_xlim(-1.15, 1.35)
ax.set_ylim(-0.85, 1.25)

# ---------------------------------------------------------------- (b) surface
ax = fig.add_subplot(1, 2, 2, projection='3d')
nz = 60
z = np.linspace(0, 1, nz)
theta = np.deg2rad(np.linspace(72.0, 50.0, nz))      # listric dip
dev = np.cumsum(1.0/np.tan(theta))*(z[1] - z[0])
U, Z = np.meshgrid(u, z)
V = np.tile(g, (nz, 1)) + np.tile(dev[:, None], (1, nu))

surf = ax.plot_surface(U, V, -Z, rstride=2, cstride=8, cmap='viridis',
                       linewidth=0.2, edgecolor='gray', alpha=0.95)
ax.plot(u, g, 0.0*u, color='crimson', lw=2.5, zorder=5)
ax.text(1.02, g[-1] + 0.05, 0.12, '$v = g(u)$', color='crimson', fontsize=15)
ax.text(-1.05, 0.42, -0.78,
        '$v = g(u) + \\mathrm{dev}(z)$', fontsize=15)
ax.set_xlabel('$u$ (along base strike)', labelpad=18)
ax.set_ylabel('$v$ (strike normal)', labelpad=14)
# mplot3d does not reliably draw the z-axis label for this view, so place
# it explicitly beside the z ticks
ax.text2D(1.20, 0.52, 'depth $z$ (normalized)', transform=ax.transAxes,
          rotation=90, va='center', ha='center')
ax.set_zticks([0.0, -0.25, -0.5, -0.75, -1.0])
ax.set_zticklabels(['0', '0.25', '0.5', '0.75', '1'])
ax.set_title('(b)')
ax.view_init(elev=22, azim=-58)
ax.set_box_aspect((1.6, 1.0, 0.7))

plt.tight_layout()
plt.savefig('../fig_fault_schematic.png', dpi=300, bbox_inches='tight')
print('saved fig_fault_schematic.png')
