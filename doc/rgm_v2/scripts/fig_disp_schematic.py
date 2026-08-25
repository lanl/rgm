#
# Schematic of the spatially varying displacement construction:
# (a) elliptical slip patch alpha(u, z) on the fault surface, with the
#     tip line, patch center, and semi-axes;
# (b) cross-section: the displaced block moves by delta_max*alpha at the
#     fault and the motion optionally decays away from the fault with the
#     Gaussian factor beta(v_perp), producing drag/rollover.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe

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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.0))

# ------------------------------------------------------------- (a) slip patch
nu, nz = 400, 300
u = np.linspace(0, 1, nu)
z = np.linspace(0, 1, nz)
U, Z = np.meshgrid(u, z)
uc, zc, au, az = 0.52, 0.45, 0.38, 0.34
alpha = np.maximum(0.0, 1.0 - ((U - uc)/au)**2 - ((Z - zc)/az)**2)

im = ax1.contourf(U, Z, alpha, levels=21, cmap='inferno')
ax1.contour(U, Z, alpha, levels=[1e-6], colors='w', linewidths=2.0)
ax1.plot(uc, zc, 'w+', ms=12, mew=2.2)
ax1.annotate('', xy=(uc + au, zc), xytext=(uc, zc),
             arrowprops=dict(arrowstyle='-|>', color='w', lw=1.6))
ax1.annotate('', xy=(uc, zc + az), xytext=(uc, zc),
             arrowprops=dict(arrowstyle='-|>', color='w', lw=1.6))
ax1.text(uc + 0.55*au, zc - 0.05, '$a_u$', color='w', fontsize=18,
         path_effects=[pe.withStroke(linewidth=2.5, foreground='black')])
ax1.text(uc + 0.02, zc + 0.55*az, '$a_z$', color='w', fontsize=18,
         path_effects=[pe.withStroke(linewidth=2.5, foreground='black')])
ax1.text(uc - 0.16, zc - 0.06, '$(u_c, z_c)$', color='w', fontsize=16,
         path_effects=[pe.withStroke(linewidth=2.5, foreground='black')])
ax1.annotate('tip line ($\\alpha = 0$)',
             xy=(uc - 0.72*au, zc - 0.66*az), xytext=(0.03, 0.08),
             color='w', fontsize=16,
             path_effects=[pe.withStroke(linewidth=2.5, foreground='black')],
             arrowprops=dict(arrowstyle='-|>', color='w', lw=1.2))
ax1.set_xlabel('$u$ (along strike)')
ax1.set_ylabel('depth $z$')
ax1.invert_yaxis()
ax1.set_title('(a)')
cb = fig.colorbar(im, ax=ax1, shrink=0.9, pad=0.02)
cb.set_label('$\\alpha = \\delta/\\delta_{\\max}$')

# --------------------------------------------------- (b) fault-normal decay
w = 0.30
theta = np.deg2rad(65.0)
zz = np.linspace(0, 1, 200)
xf = 0.5 + (zz - 0.5)/np.tan(theta)         # dipping fault trace

ax2.plot(xf, zz, color='crimson', lw=3, zorder=4)

# layer markers: footwall flat, hanging wall offset + decayed
delta0 = 0.16
for z0 in np.linspace(0.15, 0.85, 5):
    xfz = 0.5 + (z0 - 0.5)/np.tan(theta)
    xl = np.linspace(0.0, xfz, 60)
    ax2.plot(xl, np.full_like(xl, z0), color='steelblue', lw=1.8)
    xr = np.linspace(xfz, 1.25, 120)
    vperp = (xr - xfz)*np.sin(theta)
    beta = np.exp(-0.5*(vperp/w)**2)
    ax2.plot(xr, z0 + delta0*beta, color='steelblue', lw=1.8)

# slip arrows on the hanging wall, decaying with distance
z0 = 0.32
xfz = 0.5 + (z0 - 0.5)/np.tan(theta)
for xa in [0.06, 0.28, 0.55]:
    vperp = xa*np.sin(theta)
    beta = np.exp(-0.5*(vperp/w)**2)
    ax2.annotate('', xy=(xfz + xa, z0 + delta0*beta),
                 xytext=(xfz + xa, z0),
                 arrowprops=dict(arrowstyle='-|>', color='k', lw=1.5))
ax2.text(0.88, 0.055, '$\\delta(\\mathbf{y}) = \\delta_{\\max}\\,'
         '\\alpha(u, z)\\,\\beta(v_\\perp)$', fontsize=17)
ax2.text(0.92, 0.445, '$\\beta(v_\\perp) = e^{-v_\\perp^2/2w^2}$',
         fontsize=17, color='k')

ax2.text(0.03, 0.28, 'footwall (fixed)', fontsize=16, color='steelblue')
ax2.text(0.83, 0.965, 'hanging wall', fontsize=16, color='steelblue')
ax2.text(0.5 + (0.86 - 0.5)/np.tan(theta) - 0.105, 0.855, 'fault',
         color='crimson', fontsize=17, ha='center', va='center',
         rotation=-np.rad2deg(np.pi/2 - theta) - 14,
         bbox=dict(facecolor='white', alpha=0.85, edgecolor='none',
                   boxstyle='round,pad=0.15'))

ax2.set_xlabel('distance normal to strike')
ax2.set_ylabel('depth $z$')
ax2.set_xlim(0, 1.45)
ax2.set_ylim(0, 1)
ax2.invert_yaxis()
ax2.set_title('(b)')

plt.tight_layout()
plt.savefig('../fig_disp_schematic.png', dpi=200, bbox_inches='tight')
print('saved fig_disp_schematic.png')
