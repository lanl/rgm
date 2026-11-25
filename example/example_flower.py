#
# Models with flower-like branched fault systems: two systems, two dense
# systems, and three systems on folded layers. The panels (velocity, and the
# image with the fault dip on top) are plotted with pymplot into the folder
# flower, and assembled into the figure flower_2d.png.
#

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import rgm
from pymplot import showmatrix

os.makedirs('flower', exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 16})

# ==================================================================
# Models

models = [
    ('Two flower-like branched\nfault systems',
     rgm.rgm2(n1=240, n2=300, nl=40, nf=0, seed=9104, lwv=0.3, lwh=0.1,
              refl_shape='random', refl_smooth=40, refl_height=[0.0, 25.0], refl_height_top=[0.0, 5.0],
              nflower=2, flower_type='negative', flower_z=[0.45, 0.55], flower_nbranch=[4, 4],
              flower_branch_dip=[68.0, 78.0], flower_disp=[5.0, 9.0], flower_master_disp=[2.0, 3.0],
              noise_level=0.0)),
    ('Two dense flower-like\nbranched fault systems',
     rgm.rgm2(n1=240, n2=300, nl=40, nf=0, seed=9105, lwv=0.3, lwh=0.1,
              refl_shape='random', refl_smooth=40, refl_height=[0.0, 20.0], refl_height_top=[0.0, 5.0],
              nflower=2, flower_type='negative', flower_z=[0.6, 0.7], flower_nbranch=[6, 6],
              flower_spacing=[0.06, 0.09], flower_decay_width=[0.06, 0.09],
              flower_disp=[5.0, 8.0], flower_master_disp=[2.0, 3.0],
              noise_level=0.0)),
    ('Three flower-like branched\nfault systems on folded layers',
     rgm.rgm2(n1=240, n2=300, nl=40, nf=0, seed=9106, lwv=0.3, lwh=0.1,
              refl_shape='random', refl_smooth=60, refl_height=[0.0, 40.0], refl_height_top=[0.0, 5.0],
              nflower=3, flower_type='random', flower_z=[0.5, 0.9], flower_nbranch=[2, 5],
              flower_disp=[5.0, 10.0], flower_master_disp=[1.0, 2.0],
              noise_level=0.0)),
]

for i, (name, p) in enumerate(models):
    p.generate()
    clip = 2*p.image.std()
    showmatrix(p.vp, colormap='jet', cmin=2000, cmax=4000, imageonly=True,
               outfile='flower/2d_%02d_vp.png' % i)
    # The fault dip on the faults only (not a number elsewhere, which is not drawn)
    showmatrix(np.where(p.fault > 0, p.fault_dip, np.nan), background=p.image, backcolormap='binary', backclip=clip,
               colormap='jet', cmin=50, cmax=130, imageonly=True,
               outfile='flower/2d_%02d_dip.png' % i)
    p.free()

# ==================================================================
# Figure: one row per model and one column per quantity, with the quantity
# and its unit at the colorbar. The layout is in inches: panel width pw and
# height ph, row-label width lw, top margin tm, height cb of the band with
# the colorbars, and gap between the panels.

tags = ['vp', 'dip']
units = ['Vp (m/s)', 'Fault dip (degrees)']
ranges = [[2000, 3000, 4000], [50, 70, 90, 110, 130]]

nrow, ncol = len(models), len(tags)
panels = [[plt.imread('flower/2d_%02d_%s.png' % (i, tag)) for tag in tags] for i in range(nrow)]

pw = 4.5
ph = pw*panels[0][0].shape[0]/panels[0][0].shape[1]
lw, tm, cb, gap = 1.0, 0.1, 1.25, 0.08
fw = lw + ncol*pw + (ncol - 1)*gap + 0.1
fh = tm + nrow*ph + (nrow - 1)*gap + cb
fig = plt.figure(figsize=(fw, fh))
for i in range(nrow):
    y = fh - tm - (i + 1)*ph - i*gap
    for j in range(ncol):
        ax = fig.add_axes([(lw + j*(pw + gap))/fw, y/fh, pw/fw, ph/fh])
        ax.imshow(panels[i][j], interpolation='antialiased')
        ax.set_xticks([])
        ax.set_yticks([])
    fig.text((lw - 0.35)/fw, (y + 0.5*ph)/fh, models[i][0], rotation=90, ha='center', va='center')
for j in range(ncol):
    x = lw + j*(pw + gap)
    cax = fig.add_axes([(x + 0.1*pw)/fw, 0.8/fh, 0.8*pw/fw, 0.15/fh])
    cbar = fig.colorbar(ScalarMappable(norm=Normalize(ranges[j][0], ranges[j][-1]), cmap='jet'), cax=cax,
                        orientation='horizontal', ticks=ranges[j])
    cbar.set_label(units[j])
fig.savefig('flower_2d.png', dpi=150)
plt.close(fig)
