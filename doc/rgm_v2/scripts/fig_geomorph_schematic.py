#
# Schematics of the geomorphological generators (2 x 2):
# (a) meander migration: sine-generated centerlines of increasing
#     maturity, curvature-driven migration, and neck cutoff;
# (b) channel cross-section profiles: the four W_shape options
#     (sharp V, rounded U, flat floor, Gaussian bell);
# (c) meandering canyon: centerline snapshots stacked at depth levels
#     along a V-shaped incision trajectory, aperture widening from the
#     thalweg width W to the canyon-mouth width W_canyon;
# (d) drainage: D8 steepest-descent flow on synthetic terrain, channel
#     network selected by flow accumulation, Hack-scaling widths.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.ndimage import gaussian_filter

# Arial for both text and math, with TrueType embedding for the PDF
mpl.rcParams.update({'font.size': 16,
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

fig, axes = plt.subplots(2, 2, figsize=(16, 11))

# ------------------------------------------------- (a) meander migration
ax = axes[0, 0]
ns = 2000
s = np.linspace(0, 1, ns)
nwav = 2.5

def sine_generated(omega_deg):
    th = np.deg2rad(omega_deg)*np.sin(2*np.pi*nwav*s)
    x = np.cumsum(np.cos(th))/ns
    y = np.cumsum(np.sin(th))/ns
    # uniform scaling (shape-preserving) so every stage spans [0, 1] in x
    f = 1.0/x[-1]
    return f*x, f*(y - y.mean())

stages = [(45, 'lightsteelblue', 1.8, 'seed'),
          (78, 'cornflowerblue', 2.2, 'migrating'),
          (110, 'crimson', 2.8, 'mature')]
curves = []
for w, c, lw, lab in stages:
    x, y = sine_generated(w)
    curves.append((x, y))
    ax.plot(x, y, color=c, lw=lw, label=lab)

# migration arrows from mid stage toward mature stage at bend apexes
x1, y1 = curves[1]
x2, y2 = curves[2]
for frac in [0.1, 0.3, 0.5, 0.7, 0.9]:
    i = int(frac*ns)
    ax.annotate('', xy=(x2[i], y2[i]), xytext=(x1[i], y1[i]),
                arrowprops=dict(arrowstyle='-|>', color='gray', lw=1.5))

# neck cutoff annotation on the mature stage: closest inter-limb points
xm, ym = curves[2]
pts = np.column_stack([xm, ym])
best = (1e9, 0, 0)
for i in range(0, ns - 400, 10):
    j0 = i + 350
    d = np.hypot(pts[j0:, 0] - pts[i, 0], pts[j0:, 1] - pts[i, 1])
    k = np.argmin(d)
    if d[k] < best[0]:
        best = (d[k], i, j0 + k)
_, i, j = best
ax.annotate('', xy=(xm[j], ym[j]), xytext=(xm[i], ym[i]),
            arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.8))
ax.text(0.5*(xm[i] + xm[j]), max(ym[i], ym[j]) + 0.06,
        'neck cutoff when $< W$', fontsize=15, ha='center',
        bbox=dict(facecolor='white', alpha=0.85, edgecolor='none',
                  boxstyle='round,pad=0.25'))

# horizontal legend placed in a clear band below the centerlines
ax.legend(loc='lower center', ncol=3, fontsize=14, frameon=False,
          handlelength=1.8, columnspacing=2.0, borderaxespad=0.4)
ax.set_xlabel('Along-valley distance (normalized)')
ax.set_ylabel('Lateral position')
ax.set_title('(a)')
ax.set_aspect('equal')
ax.set_ylim(-0.44, 0.36)

# ----------------------------------- (b) channel cross-section profiles
ax = axes[0, 1]
t = np.linspace(0, 1, 300)                          # 0 = floor, 1 = banks
Wb = 0.07                                           # floor width fraction

def prof_power(e):
    return t**e

def prof_smoothstep():
    return 3*t**2 - 2*t**3

def prof_bell(sig=0.5):
    e0 = np.exp(-1.0/sig**2)
    return sig*np.sqrt(-np.log((1 - t)*(1 - e0) + e0))

profiles = [(prof_power(1.0), 'sharp V bed\n(power law,\nexponent 1)'),
            (prof_power(0.5), 'rounded U bed\n(power law,\nexponent 0.5)'),
            (prof_smoothstep(), 'flat bed,\nsteep banks\n(smoothstep)'),
            (prof_bell()/prof_bell()[-1], 'bell\n(Gaussian)')]

for k, (f, lab) in enumerate(profiles):
    xc = 1.7*k
    half = 0.5*(Wb + (1 - Wb)*f)
    xs = np.concatenate([xc - half[::-1], xc + half])
    ys = np.concatenate([(1 - t)[::-1], 1 - t])
    ax.plot(xs, ys, color='saddlebrown', lw=2.4)
    ax.fill_between(xs, ys, 0.0, color='#f1e3c9', zorder=0)
    ax.plot([xc - 0.85, xc - half[-1]], [0, 0], color='saddlebrown', lw=2.4)
    ax.plot([xc + half[-1], xc + 0.85], [0, 0], color='saddlebrown', lw=2.4)
    ax.text(xc, 1.12, lab, fontsize=14, ha='center', va='top')

ax.annotate('', xy=(0.5, -0.08), xytext=(-0.5, -0.08),
            arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.6))
ax.text(0, -0.14, '$W$', fontsize=16, ha='center', va='bottom')
ax.set_ylim(1.55, -0.3)
ax.set_xlim(-1.0, 1.7*3 + 1.0)
ax.set_xticks([])
ax.set_xlabel('Across-channel distance')
ax.set_ylabel('Depth below banks')
ax.set_title('(b)')
ax.set_yticks([0, 0.5, 1.0])

# ------------------------------------- (c) canyon: stacked snapshots on V
ax = axes[1, 0]
rng = np.random.default_rng(3)
K = 9
tt = np.linspace(0, 1, K)
zfrac = 1 - 2*np.abs(tt - 0.5)                     # V-trajectory
walk = np.cumsum(rng.normal(0, 0.3, K))
walk -= walk.mean()
walk *= 0.55/max(1e-9, np.max(np.abs(walk)))
W, Wc = 0.08, 0.36

ny, nz = 700, 400
yy = np.linspace(-1, 1, ny)
zz = np.linspace(0, 1, nz)                          # depth, 0 = surface
eroded = np.zeros((nz, ny), dtype=bool)
outlines = []
for tk, zf, c0 in zip(tt, zfrac, walk):
    zi = zf                                         # incision depth of snapshot
    if zi < 0.06:
        continue
    # aperture widens from W at the incision depth to Wc at the surface
    half = 0.5*(W + (Wc - W)*(1 - zz/zi).clip(0, 1))
    m = (np.abs(yy[None, :] - c0) <= half[:, None]) & (zz[:, None] <= zi)
    eroded |= m
    outlines.append((c0, zi, half))

ax.imshow(np.where(eroded, 1.0, 0.0), extent=(-1, 1, 1, 0), aspect='auto',
          cmap=mpl.colors.ListedColormap(['#f1e3c9', '#a8c4de']),
          interpolation='nearest')
for c0, zi, half in outlines:
    zs = zz[zz <= zi]
    hs = half[:len(zs)]
    ax.plot(c0 - hs, zs, color='steelblue', lw=0.9, alpha=0.85)
    ax.plot(c0 + hs, zs, color='steelblue', lw=0.9, alpha=0.85)

imax = np.argmax(zfrac)
c_deep = walk[imax]
ax.annotate('', xy=(c_deep + 0.5*W, 0.96), xytext=(c_deep - 0.5*W, 0.96),
            arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.5))
ax.annotate('thalweg width $W$', xy=(c_deep + 0.5*W + 0.01, 0.96),
            xytext=(0.40, 0.82), fontsize=15,
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='none',
                      boxstyle='round,pad=0.25'),
            arrowprops=dict(arrowstyle='-|>', color='black', lw=1.3))
ax.annotate('', xy=(c_deep + 0.5*Wc, 0.04), xytext=(c_deep - 0.5*Wc, 0.04),
            arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.5))
ax.text(c_deep + 0.5*Wc + 0.06, 0.07, '$W_{\\mathrm{canyon}}$', fontsize=16,
        bbox=dict(facecolor='white', alpha=0.85, edgecolor='none',
                  boxstyle='round,pad=0.25'))
ax.plot(walk, zfrac, 'o-', color='crimson', ms=5, lw=1.6)
ax.text(-0.96, 0.32, 'snapshots on the\nV-trajectory\n$z_{\\mathrm{frac}}(t) = 1 - 2|t - 1/2|$',
        fontsize=15, color='crimson',
        bbox=dict(facecolor='white', alpha=0.85, edgecolor='none',
                  boxstyle='round,pad=0.3'))
ax.set_xlabel('Lateral position (normalized)')
ax.set_ylabel('Depth (normalized)')
ax.set_title('(c)')
ax.set_xlim(-1, 1)
ax.set_ylim(1, 0)

# --------------------------------------------------- (d) drainage pipeline
ax = axes[1, 1]
n = 18
rng = np.random.default_rng(12)
dem = gaussian_filter(rng.normal(size=(n, n)), 1.8)
dem = (dem - dem.min())/(dem.max() - dem.min())
dem += 0.8*np.linspace(1, 0, n)[:, None]            # regional tilt (drain downward)

DX8 = [1, 1, 0, -1, -1, -1, 0, 1]
DY8 = [0, 1, 1, 1, 0, -1, -1, -1]
flow = -np.ones((n, n, 2), dtype=int)
for i in range(n):
    for j in range(n):
        best, bi, bj = 0.0, -1, -1
        for dx, dy in zip(DX8, DY8):
            i2, j2 = i + dy, j + dx
            if 0 <= i2 < n and 0 <= j2 < n:
                drop = (dem[i, j] - dem[i2, j2])/np.hypot(dx, dy)
                if drop > best:
                    best, bi, bj = drop, i2, j2
        flow[i, j] = (bi, bj)

acc = np.ones((n, n))
order = np.argsort(dem, axis=None)[::-1]
for k in order:
    i, j = divmod(k, n)
    bi, bj = flow[i, j]
    if bi >= 0:
        acc[bi, bj] += acc[i, j]

ax.imshow(dem, cmap='gist_earth', origin='upper', alpha=0.75,
          extent=(-0.5, n - 0.5, n - 0.5, -0.5))
for i in range(0, n, 2):
    for j in range(0, n, 2):
        bi, bj = flow[i, j]
        if bi >= 0:
            ax.annotate('', xy=(bj, bi), xytext=(j, i),
                        arrowprops=dict(arrowstyle='-|>', color='k',
                                        lw=0.8, alpha=0.45))
thresh = np.quantile(acc, 0.86)
amax = acc.max()
for i in range(n):
    for j in range(n):
        if acc[i, j] >= thresh:
            bi, bj = flow[i, j]
            if bi >= 0:
                ax.plot([j, bj], [i, bi], color='crimson',
                        lw=1.0 + 5.5*np.sqrt(acc[i, j]/amax),
                        solid_capstyle='round')
ax.text(0.2, 1.0, 'width, depth $\\propto \\sqrt{A/A_{\\max}}$',
        fontsize=16, color='crimson',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none',
                  boxstyle='round,pad=0.3'))
ax.set_xlabel('$x$ (cells)')
ax.set_ylabel('$y$ (cells)')
ax.set_title('(d)')
ax.set_xlim(-0.5, n - 0.5)
ax.set_ylim(n - 0.5, -0.5)

plt.tight_layout()
plt.savefig('../fig_geomorph_schematic.png', dpi=300, bbox_inches='tight')
print('saved fig_geomorph_schematic.png')
