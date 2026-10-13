#
# One model for each type of fault that rgm2 and rgm3 generate. The panels
# (velocity, image, fault dip, fault displacement, in 3D fault strike, and RGT)
# are plotted with pymplot into the folder fault_types, and assembled into the
# figures fault_types_2d.png and fault_types_3d.png.
#

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import rgm
from pymplot import showmatrix, showslice

os.makedirs('fault_types', exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 16})

# ==================================================================
# 2D models

models2 = [
    ('Planar\nfaults',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=4, seed=211, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], disp=[10.0, 20.0], delta_dip=[0.0, 0.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Listric\nfaults',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=4, seed=212, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], disp=[10.0, 20.0], delta_dip=[20.0, 30.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Quasi-parallel\nfaults',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=6, seed=213, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[65.0, 115.0], disp=[8.0, 12.0], delta_dip=[0.0, 10.0], yn_regular_fault=True,
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Faults that\ndie out',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=4, seed=214, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], disp=[15.0, 25.0], delta_dip=[0.0, 20.0],
              yn_vary_disp=True, disp_radius_dip=[0.4, 0.7],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Displacement\ndecay',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=3, seed=215, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], disp=[15.0, 25.0], delta_dip=[0.0, 20.0],
              yn_disp_decay=True, disp_decay_width=[0.3, 0.5],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Fault-propagation\nfold',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=216, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=1, thrust_vergence='right', thrust_z=[1.1, 1.2], thrust_tip_z=[0.3, 0.4],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Fault-bend\nfold',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=217, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=1, thrust_vergence='right', thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.35, 0.4],
              thrust_disp=[50.0, 60.0], yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Imbricate\nfan',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=251, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=3, thrust_vergence='right', thrust_z=[1.1, 1.2], thrust_tip_z=[0.1, 0.5],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Duplex',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=212, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=4, thrust_vergence='right', thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.5, 0.55],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Pop-up',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=0, seed=213, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=2, thrust_back=1.0, thrust_tip_z=[0.25, 0.45],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Thrusts cut\nby faults',
     rgm.rgm2(n1=151, n2=301, nl=30, nf=2, seed=218, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[55.0, 70.0], disp=[8.0, 12.0],
              nthrust=2, thrust_vergence='right', thrust_z=[1.1, 1.2], thrust_tip_z=[0.1, 0.4],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
]

for i, (name, p) in enumerate(models2):
    p.generate()
    clip = 2*p.image.std()
    showmatrix(p.vp, colormap='jet', cmin=2000, cmax=4000, imageonly=True,
               outfile='fault_types/2d_%02d_vp.png' % i)
    showmatrix(p.image, colormap='binary', clip=clip, imageonly=True,
               outfile='fault_types/2d_%02d_image.png' % i)
    showmatrix(p.fault_dip, background=p.image, backcolormap='binary', backclip=clip,
               colormap='jet', cmin=0, cmax=180, alphas='0:0,0.5:1', imageonly=True,
               outfile='fault_types/2d_%02d_dip.png' % i)
    showmatrix(abs(p.fault_disp), background=p.image, backcolormap='binary', backclip=clip,
               colormap='jet', cmin=0, cmax=40, alphas='0:0,0.5:1', imageonly=True,
               outfile='fault_types/2d_%02d_disp.png' % i)
    showmatrix(p.rgt, colormap='magma', cmin=0, cmax=1, imageonly=True,
               outfile='fault_types/2d_%02d_rgt.png' % i)
    p.free()

# ==================================================================
# 3D models

models3 = [
    ('Planar\nfaults',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=4, seed=311, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], strike=[0.0, 180.0], rake=[70.0, 110.0], disp=[10.0, 20.0], delta_dip=[0.0, 0.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Listric\nfaults',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=4, seed=312, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], strike=[0.0, 180.0], rake=[70.0, 110.0], disp=[10.0, 20.0], delta_dip=[20.0, 30.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Curved\nstrike',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=4, seed=313, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], strike=[0.0, 180.0], rake=[70.0, 110.0], disp=[10.0, 20.0], delta_dip=[0.0, 15.0],
              delta_strike=[20.0, 30.0], yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Quasi-parallel\nfaults',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=6, seed=314, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[65.0, 115.0], strike=[40.0, 40.0], rake=[90.0, 90.0], disp=[8.0, 12.0], delta_dip=[0.0, 10.0],
              yn_regular_fault=True, yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Faults that\ndie out',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=4, seed=315, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], strike=[0.0, 180.0], rake=[70.0, 110.0], disp=[15.0, 25.0], delta_dip=[0.0, 15.0],
              yn_vary_disp=True, disp_radius_strike=[0.4, 0.6], disp_radius_dip=[0.5, 0.8],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Displacement\ndecay',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=3, seed=316, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 120.0], strike=[0.0, 180.0], rake=[70.0, 110.0], disp=[15.0, 25.0], delta_dip=[0.0, 15.0],
              yn_disp_decay=True, disp_decay_width=[0.3, 0.5], yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Strike-slip\nfaults',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=2, seed=317, lwv=0.4, lwh=0.2,
              refl_shape='perlin', refl_smooth=3, refl_height=[0.0, 30.0], refl_height_top=[0.0, 4.0],
              dip=[80.0, 100.0], strike=[20.0, 40.0], rake=[0.0, 0.0], disp=[20.0, 30.0], delta_dip=[0.0, 10.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Blind\nthrust',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=214, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=1, thrust_z=[1.1, 1.2], thrust_tip_z=[0.3, 0.4], thrust_disp=[25.0, 30.0],
              thrust_radius=[0.35, 0.4], thrust_delta_strike=[10.0, 15.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Fault-bend\nfold',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=318, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=1, thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.35, 0.4], thrust_disp=[30.0, 35.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Imbricate\nfan',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=319, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=3, thrust_z=[1.1, 1.2], thrust_tip_z=[0.1, 0.5], thrust_delta_strike=[8.0, 12.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Duplex',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=320, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=3, thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.5, 0.55], thrust_disp=[20.0, 25.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Pop-up',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=0, seed=321, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              nthrust=2, thrust_back=1.0, thrust_tip_z=[0.25, 0.4], thrust_delta_strike=[5.0, 10.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
    ('Thrusts cut\nby faults',
     rgm.rgm3(n1=128, n2=160, n3=160, nl=30, nf=2, seed=215, lwv=0.4, lwh=0.1,
              refl_shape='random', refl_smooth=30, refl_height=[0.0, 8.0], refl_height_top=[0.0, 4.0],
              dip=[60.0, 75.0], disp=[6.0, 10.0],
              nthrust=3, thrust_roof=1.0, thrust_z=[0.85, 0.9], thrust_roof_z=[0.5, 0.55], thrust_disp=[20.0, 25.0],
              yn_dip_independent=True, yn_rgt=True, noise_level=0.0)),
]

# The 3D panels show a depth slice at three quarters of the depth in front of
# the two vertical slices at the far sides of the volume. showslice with render='3d'
# ends with exit() after writing the file, which is caught here to continue.
for i, (name, q) in enumerate(models3):
    q.generate()
    clip = 2*q.image.std()
    for data, colormap, cmin, cmax, tag in [(q.vp, 'jet', 2000, 4000, 'vp'),
                                            (q.image, 'binary', -clip, clip, 'image'),
                                            (q.fault_dip, 'jet', 0, 180, 'dip'),
                                            (abs(q.fault_disp), 'jet', 0, 40, 'disp'),
                                            (q.fault_strike, 'jet', 0, 180, 'strike'),
                                            (q.rgt, 'magma', 0, 1, 'rgt')]:
        try:
            if tag in ('vp', 'image', 'rgt'):
                showslice(data, colormap=colormap, cmin=cmin, cmax=cmax, render='3d',
                          slice1=3*q.vp.shape[0]//4, slice2=q.vp.shape[1] - 2, slice3=q.vp.shape[2] - 2,
                          outfile='fault_types/3d_%02d_%s.png' % (i, tag))
            else:
                showslice(data, background=q.image, backcolormap='binary', backclip=clip,
                          colormap=colormap, cmin=cmin, cmax=cmax, alphas='0:0,0.5:1', render='3d',
                          slice1=3*q.vp.shape[0]//4, slice2=q.vp.shape[1] - 2, slice3=q.vp.shape[2] - 2,
                          outfile='fault_types/3d_%02d_%s.png' % (i, tag))
        except SystemExit:
            pass
    q.free()

# ==================================================================
# Figures: one row per fault type and one column per quantity, with the 3D
# renderings cropped to the box of their content. The colorbars span the ticks
# given for each quantity; RGT uses magma, the other quantities jet. The
# layout is in inches: panel width pw and height ph, row-label width lw, title
# height th, height cb of the band with the colorbars, and gap between the panels.

for dim, models, tags, titles, ranges in [
        ('2d', models2, ['vp', 'image', 'dip', 'disp', 'rgt'],
         ['Vp (m/s)', 'Image', 'Fault dip (degrees)', 'Fault displacement', 'RGT'],
         [[2000, 3000, 4000], None, [0, 60, 120, 180], [0, 10, 20, 30, 40], [0, 0.5, 1]]),
        ('3d', models3, ['vp', 'image', 'dip', 'disp', 'strike', 'rgt'],
         ['Vp (m/s)', 'Image', 'Fault dip (degrees)', 'Fault displacement', 'Fault strike (degrees)', 'RGT'],
         [[2000, 3000, 4000], None, [0, 60, 120, 180], [0, 10, 20, 30, 40], [0, 60, 120, 180], [0, 0.5, 1]])]:

    nrow, ncol = len(models), len(tags)
    panels = [[plt.imread('fault_types/%s_%02d_%s.png' % (dim, i, tag)) for tag in tags] for i in range(nrow)]
    if dim == '3d':
        r0, r1, c0, c1 = 10000, 0, 10000, 0
        for row in panels:
            for img in row:
                r, c = (img[:, :, :3].min(axis=2) < 0.98).nonzero()
                r0, r1, c0, c1 = min(r0, r.min()), max(r1, r.max()), min(c0, c.min()), max(c1, c.max())
        panels = [[img[r0:r1 + 1, c0:c1 + 1] for img in row] for row in panels]

    pw = 3.6 if dim == '2d' else 3.0
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
            cmap = 'magma' if tags[j] == 'rgt' else 'jet'
            fig.colorbar(ScalarMappable(norm=Normalize(ranges[j][0], ranges[j][-1]), cmap=cmap), cax=cax,
                         orientation='horizontal', ticks=ranges[j])
    fig.savefig('fault_types_%s.png' % dim, dpi=150)
    plt.close(fig)
