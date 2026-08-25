#
# Paired velocity/seismic-image figures for the karst examples:
#   fig_karst_3d.png  3D model, three-slice views, (a) Vp and (b) image
#   fig_karst_2d.png  2D section,                  (a) Vp and (b) image
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.image as mpimg

mpl.rcParams.update({'font.size': 12,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

exdir = '../tmp'
outdir = '../tmp'
os.makedirs(outdir, exist_ok=True)

clip = 0.06   # symmetric amplitude clip for the seismic images

# ---------------------------------------------------------------- 3D pair
n1, n2, n3 = 151, 201, 251
common = (f'-n1 {n1} -n2 {n2} -n3 {n3} -slice1 104 -slice2 106 -slice3 22 '
          f'-label1 "Z (Grid Number)" -label2 "Y (Grid Number)" '
          f'-label3 "X (Grid Number)"')
os.system(f'x_showslice -in {exdir}/example_3d_vp_karst_1_151x201x251.bin '
          f'{common} -color jet -legend 1 -unit "Vp (m/s)" '
          f'-out {outdir}/karst3d_vp.png >/dev/null 2>&1')
os.system(f'x_showslice -in {exdir}/example_3d_image_karst_1_151x201x251.bin '
          f'{common} -color binary -legend 1 -unit "Amplitude" '
          f'-cmin {-clip:.6f} -cmax {clip:.6f} '
          f'-out {outdir}/karst3d_img.png >/dev/null 2>&1')

fig, axes = plt.subplots(2, 1, figsize=(9.0, 12.6))
for ax, kind, lab in [(axes[0], 'vp', '(a)'), (axes[1], 'img', '(b)')]:
    ax.imshow(mpimg.imread(f'{outdir}/karst3d_{kind}.png'))
    ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=18, va='bottom')
    ax.axis('off')
plt.subplots_adjust(hspace=0.06, left=0.01, right=0.99, top=0.97, bottom=0.01)
plt.savefig('../fig_karst_3d.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print('saved fig_karst_3d.png')

# ---------------------------------------------------------------- 2D pair
# rendered with x_showmatrix, which auto-scales to the natural data
# aspect ratio when the plot sizes are omitted
m1, m2 = 151, 301
common2d = (f'-n1 {m1} -n2 {m2} -label1 "Z (Grid Number)" '
            f'-label2 "X (Grid Number)"')
os.system(f'x_showmatrix -in {exdir}/example_2d_vp_karst_2_151x301.bin '
          f'{common2d} -color jet -legend 1 -unit "Vp (m/s)" '
          f'-out {outdir}/karst2d_vp.png >/dev/null 2>&1')
os.system(f'x_showmatrix -in {exdir}/example_2d_image_karst_2_151x301.bin '
          f'{common2d} -color binary -legend 1 -unit "Amplitude" '
          f'-cmin {-clip:.6f} -cmax {clip:.6f} '
          f'-out {outdir}/karst2d_img.png >/dev/null 2>&1')

fig, axes = plt.subplots(2, 1, figsize=(9.0, 8.4))
for ax, kind, lab in [(axes[0], 'vp', '(a)'), (axes[1], 'img', '(b)')]:
    ax.imshow(mpimg.imread(f'{outdir}/karst2d_{kind}.png'))
    ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=18, va='bottom')
    ax.axis('off')
plt.subplots_adjust(hspace=0.05, left=0.01, right=0.99, top=0.97, bottom=0.01)
plt.savefig('../fig_karst_2d.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print('saved fig_karst_2d.png')
