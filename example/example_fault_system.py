#
# Fault systems that combine thrust faults, flower-like branched faults, and
# faults, in three 2D models. The panels (velocity, the fault indices on the
# relative geologic time, and the fault dip on the seismic image) are
# assembled into the figure fault_system_2d.png.
#

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm
import rgm
from pymplot import showmatrix

os.makedirs('fault_system', exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 16})

# ==================================================================
# Models

models = [
    ('Thrusts, flower-like\nbranched faults, faults,\nand an unconformity',
     rgm.rgm2(n1=160, n2=320, nl=40, nf=3, seed=3026, lwv=0.4, lwh=0.1, delta_v=1000.0,
              refl_shape='random', refl_smooth=40, refl_height=[0.0, 12.0], refl_height_top=[0.0, 6.0],
              dip=[55.0, 125.0], disp=[6.0, 15.0],
              nthrust=3,
              nflower=2, flower_type='positive', flower_nbranch=[2, 5], flower_master_disp=[1.0, 1.5],
              unconf=1, unconf_z=[0.12, 0.18], unconf_height=[5.0, 10.0],
              yn_rgt=True, noise_level=0.0)),
    ('Thrusts, flower-like\nbranched faults, and\nfaults that die out',
     rgm.rgm2(n1=160, n2=320, nl=40, nf=3, seed=3028, lwv=0.4, lwh=0.1, delta_v=1000.0,
              refl_shape='random', refl_smooth=40, refl_height=[0.0, 12.0], refl_height_top=[0.0, 6.0],
              dip=[55.0, 125.0], disp=[6.0, 15.0], yn_vary_disp=True,
              nthrust=2, thrust_sequence='random',
              nflower=2, flower_type='hybrid', flower_nbranch=[2, 5], flower_master_disp=[1.0, 2.25],
              yn_rgt=True, noise_level=0.0)),
    ('Thrusts and flower-like\nbranched faults',
     rgm.rgm2(n1=160, n2=320, nl=40, nf=0, seed=3015, lwv=0.4, lwh=0.1, delta_v=1000.0,
              refl_shape='random', refl_smooth=40, refl_height=[0.0, 12.0], refl_height_top=[0.0, 6.0],
              nthrust=3, thrust_sequence='random',
              nflower=1, flower_type='positive', flower_nbranch=[1, 2], flower_master_disp=[1.0, 1.75],
              yn_rgt=True, noise_level=0.0)),
]

rgt = []
fault = []
for i, (name, p) in enumerate(models):
    p.generate()
    showmatrix(p.vp, colormap='jet', cmin=2000, cmax=4000, imageonly=True,
               outfile='fault_system/2d_%02d_vp.png' % i)
    # The fault dip on the faults only (not a number elsewhere, which is not drawn)
    showmatrix(np.where(p.fault > 0, p.fault_dip, np.nan), background=p.image, backcolormap='binary',
               backclip=2*p.image.std(), colormap='jet', cmin=0, cmax=180, imageonly=True,
               outfile='fault_system/2d_%02d_dip.png' % i)
    rgt.append(p.rgt.copy())
    fault.append(p.fault.copy())
    p.free()

# ==================================================================
# Figure: one row per model; the velocity, the fault indices in pure colors on
# the relative geologic time in light gray (pymplot takes only named
# colormaps, so this panel is drawn with matplotlib), and the fault dip on the
# seismic image. The quantities and units are at the colorbars. The layout is
# in inches: panel width pw and height ph, row-label width lw, top margin tm,
# height cb of the band with the colorbars, and gap between the panels; at 192
# dots per inch, each grid point of the fault-index panels covers 3 by 3
# pixels.

colors = ['#ff0000', '#0000ff', '#00a000', '#ff00ff', '#00c0c0', '#ff8000', '#8000ff', '#000000', '#a05000',
          '#008080', '#c0c000', '#ff0080', '#0080ff', '#80ff00', '#800000', '#000080', '#008000', '#808000']
nid = max(int(f.max()) for f in fault)
fmap = ListedColormap([colors[k % len(colors)] for k in range(nid)])
fnorm = BoundaryNorm(np.arange(0.5, nid + 1.0), nid)

nrow, ncol = len(models), 3
pw = 5.0
ph = pw*fault[0].shape[0]/fault[0].shape[1]
lw, tm, cb, gap = 1.25, 0.125, 1.25, 0.0625
fw = lw + ncol*pw + (ncol - 1)*gap + 0.125
fh = tm + nrow*ph + (nrow - 1)*gap + cb
fig = plt.figure(figsize=(fw, fh))
for i in range(nrow):
    y = fh - tm - (i + 1)*ph - i*gap
    for j in range(ncol):
        ax = fig.add_axes([(lw + j*(pw + gap))/fw, y/fh, pw/fw, ph/fh])
        if j == 0:
            ax.imshow(plt.imread('fault_system/2d_%02d_vp.png' % i), interpolation='antialiased', aspect='auto')
        elif j == 1:
            ax.imshow(rgt[i], cmap='gray', alpha=0.3, interpolation='nearest', aspect='auto')
            ax.imshow(np.ma.masked_where(fault[i] == 0, fault[i]), cmap=fmap, norm=fnorm, interpolation='nearest',
                      aspect='auto')
        else:
            ax.imshow(plt.imread('fault_system/2d_%02d_dip.png' % i), interpolation='antialiased', aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
    fig.text((lw - 0.5)/fw, (y + 0.5*ph)/fh, models[i][0], rotation=90, ha='center', va='center')

bars = [(Normalize(2000, 4000), 'jet', [2000, 3000, 4000], 'Vp (m/s)'),
        (fnorm, fmap, [1] + list(range(5, nid + 1, 5)), 'Fault index'),
        (Normalize(0, 180), 'jet', [0, 60, 120, 180], 'Fault dip (degrees)')]
for j, (norm, cmap, ticks, label) in enumerate(bars):
    cax = fig.add_axes([(lw + j*(pw + gap) + 0.1*pw)/fw, 0.8/fh, 0.8*pw/fw, 0.15/fh])
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', ticks=ticks)
    cbar.minorticks_off()
    cbar.set_label(label)
fig.savefig('fault_system_2d.png', dpi=192)
plt.close(fig)
