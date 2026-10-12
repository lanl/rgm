#
# Elastic models of thrust faults in 2D and 3D: a pop-up, below an unconformity
# in 3D, and a duplex. Each layer has its own random Vp/Vs. The panels (Vp, Vs,
# Vp/Vs, and the PP, PS and SS images) are plotted with pymplot into the folder
# elastic, and assembled into the figures elastic_2d.png and elastic_3d.png.
#

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import rgm
from pymplot import showmatrix, showslice

os.makedirs('elastic', exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 16})

# ==================================================================
# 2D models

models2 = [
    ('Pop-up',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=702, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=2, thrust_back=1.0, thrust_tip_z=[0.25, 0.45],
              yn_elastic=True, yn_dip_independent=True, noise_level=0.0)),
    ('Duplex',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=713, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=4, thrust_vergence='right', thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.5, 0.55],
              yn_elastic=True, yn_dip_independent=True, noise_level=0.0)),
]

for i, (name, p) in enumerate(models2):
    p.generate()
    showmatrix(p.vp, colormap='jet', cmin=2000, cmax=4000, imageonly=True,
               outfile='elastic/2d_%02d_vp.png' % i)
    showmatrix(p.vs, colormap='jet', cmin=1200, cmax=2400, imageonly=True,
               outfile='elastic/2d_%02d_vs.png' % i)
    showmatrix(p.vp/p.vs, colormap='jet', cmin=1.5, cmax=1.8, imageonly=True,
               outfile='elastic/2d_%02d_vpvs.png' % i)
    showmatrix(p.image_pp, colormap='binary', clip=2*p.image_pp.std(), imageonly=True,
               outfile='elastic/2d_%02d_pp.png' % i)
    showmatrix(p.image_ps, colormap='binary', clip=2*p.image_ps.std(), imageonly=True,
               outfile='elastic/2d_%02d_ps.png' % i)
    showmatrix(p.image_ss, colormap='binary', clip=2*p.image_ss.std(), imageonly=True,
               outfile='elastic/2d_%02d_ss.png' % i)
    p.free()

# ==================================================================
# 3D models

models3 = [
    ('Pop-up',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=361, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=2, thrust_back=1.0, thrust_tip_z=[0.25, 0.4], thrust_delta_strike=[5.0, 10.0],
              unconf=1, unconf_z=[0.2, 0.3],
              yn_elastic=True, yn_dip_independent=True, noise_level=0.0)),
    ('Duplex',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=733, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=3, thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.5, 0.55], thrust_disp=[20.0, 25.0],
              yn_elastic=True, yn_dip_independent=True, noise_level=0.0)),
]

# The 3D panels show a depth slice at three quarters of the depth in front of
# the two vertical slices at the far sides of the volume. showslice with render='3d'
# ends with exit() after writing the file, which is caught here to continue.
for i, (name, q) in enumerate(models3):
    q.generate()
    for data, colormap, cmin, cmax, tag in [(q.vp, 'jet', 2000, 4000, 'vp'),
                                            (q.vs, 'jet', 1000, 3000, 'vs'),
                                            (q.vp/q.vs, 'jet', 1.5, 1.8, 'vpvs'),
                                            (q.image_pp, 'binary', -2*q.image_pp.std(), 2*q.image_pp.std(), 'pp'),
                                            (q.image_ps, 'binary', -2*q.image_ps.std(), 2*q.image_ps.std(), 'ps'),
                                            (q.image_ss, 'binary', -2*q.image_ss.std(), 2*q.image_ss.std(), 'ss')]:
        try:
            showslice(data, colormap=colormap, cmin=cmin, cmax=cmax, render='3d',
                      slice1=3*q.vp.shape[0]//4, slice2=q.vp.shape[1] - 2, slice3=q.vp.shape[2] - 2,
                      outfile='elastic/3d_%02d_%s.png' % (i, tag))
        except SystemExit:
            pass
    q.free()

# ==================================================================
# Figures: one row per model and one column per quantity, with the 3D
# renderings cropped to the box of their content. The colorbars span the ticks
# given for each quantity. The layout is in inches: panel width pw and height
# ph, row-label width lw, title height th, height cb of the band with the
# colorbars, and gap between the panels.

tags = ['vp', 'vs', 'vpvs', 'pp', 'ps', 'ss']
titles = ['Vp (m/s)', 'Vs (m/s)', 'Vp/Vs', 'PP image', 'PS image', 'SS image']

for dim, models, ranges in [
        ('2d', models2, [[2000, 3000, 4000], [1200, 1800, 2400], [1.5, 1.6, 1.7, 1.8], None, None, None]),
        ('3d', models3, [[2000, 3000, 4000], [1000, 2000, 3000], [1.5, 1.6, 1.7, 1.8], None, None, None])]:

    nrow, ncol = len(models), len(tags)
    panels = [[plt.imread('elastic/%s_%02d_%s.png' % (dim, i, tag)) for tag in tags] for i in range(nrow)]
    if dim == '3d':
        r0, r1, c0, c1 = 10000, 0, 10000, 0
        for row in panels:
            for img in row:
                r, c = (img[:, :, :3].min(axis=2) < 0.98).nonzero()
                r0, r1, c0, c1 = min(r0, r.min()), max(r1, r.max()), min(c0, c.min()), max(c1, c.max())
        panels = [[img[r0:r1 + 1, c0:c1 + 1] for img in row] for row in panels]

    pw = 3.0
    ph = pw*panels[0][0].shape[0]/panels[0][0].shape[1]
    lw, th, cb, gap = 1.0, 0.5, 0.9, 0.08
    fw = lw + ncol*pw + (ncol - 1)*gap + 0.1
    fh = th + nrow*ph + (nrow - 1)*gap + cb
    fig = plt.figure(figsize=(fw, fh))
    for i in range(nrow):
        y = fh - th - (i + 1)*ph - i*gap
        for j in range(ncol):
            ax = fig.add_axes([(lw + j*(pw + gap))/fw, y/fh, pw/fw, ph/fh])
            ax.imshow(panels[i][j], interpolation='antialiased')
            ax.set_xticks([])
            ax.set_yticks([])
            if dim == '3d':
                ax.axis('off')
        fig.text((lw - 0.35)/fw, (y + 0.5*ph)/fh, models[i][0], rotation=90, ha='center', va='center')
    for j in range(ncol):
        x = lw + j*(pw + gap)
        fig.text((x + 0.5*pw)/fw, (fh - 0.5*th)/fh, titles[j], ha='center', va='center')
        if ranges[j] is not None:
            cax = fig.add_axes([(x + 0.1*pw)/fw, 0.45/fh, 0.8*pw/fw, 0.15/fh])
            fig.colorbar(ScalarMappable(norm=Normalize(ranges[j][0], ranges[j][-1]), cmap='jet'), cax=cax,
                         orientation='horizontal', ticks=ranges[j])
    fig.savefig('elastic_%s.png' % dim, dpi=150)
    plt.close(fig)
