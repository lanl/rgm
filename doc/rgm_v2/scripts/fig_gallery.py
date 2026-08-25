#
# 4 x 4 gallery figures of 16 random RGM v2.0 models, from the batch
# produced by rgm/tmp/gen_gallery.f90:
#   - velocity models (common color scale)
#   - seismic images (fixed amplitude clip)
#   - local fault strike overlaid on the seismic images
#   - local fault displacement overlaid on the seismic images
# Each model is rendered as a three-slice view with pymplot x_showslice,
# then montaged with matplotlib.
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.image as mpimg

mpl.rcParams.update({'font.size': 16,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

n1, n2, n3 = 128, 160, 200
tmpdir = '../tmp'
outdir = '../tmp'
os.makedirs(outdir, exist_ok=True)

# depth slice per model: through the karst zone for karst-bearing models,
# through the unconformity interval otherwise
karsty = {9, 10, 11, 12, 16}

# fixed symmetric amplitude clip so weaker reflections are visible
clip = 0.05
# larger (softer) clip for the overlay backgrounds so the fault-attribute
# fields stand out more vividly
clip_bg = 0.075
# common displacement color scale (labels reach ~14; most models 6-13)
dmax = 12

for i in range(1, 17):
    s1 = 78 if i in karsty else 40
    vp = f'{tmpdir}/gallery_vp_{i}_128x160x200.bin'
    im = f'{tmpdir}/gallery_img_{i}_128x160x200.bin'
    fs = f'{tmpdir}/gallery_fstrike_{i}_128x160x200.bin'
    fd = f'{tmpdir}/gallery_fdisp_{i}_128x160x200.bin'
    common = (f'-n1 {n1} -n2 {n2} -n3 {n3} -slice1 {s1} -slice2 80 '
              f'-slice3 100 -tick1d 40 -tick2d 40 -tick3d 50 '
              f'-label1 "Z (Grid Number)" -label2 "Y (Grid Number)" '
              f'-label3 "X (Grid Number)"')
    os.system(f'x_showslice -in {vp} {common} -color jet '
              f'-cmin 2000 -cmax 4000 -out {outdir}/vp_{i}.png >/dev/null 2>&1')
    os.system(f'x_showslice -in {im} {common} -color binary '
              f'-cmin {-clip:.6f} -cmax {clip:.6f} '
              f'-out {outdir}/img_{i}.png >/dev/null 2>&1')
    # fault attribute overlays on the seismic image background
    overlay = (f'-background {im} -backcolor binary '
               f'-backcmin {-clip_bg:.6f} -backcmax {clip_bg:.6f}')
    os.system(f'x_showslice -in {fs} {common} {overlay} -color jet '
              f'-cmin 0 -cmax 180 -alphas 0:0,0.5:0,1:1 '
              f'-out {outdir}/fstrike_{i}.png >/dev/null 2>&1')
    os.system(f'x_showslice -in {fd} {common} {overlay} -color jet '
              f'-cmin 0 -cmax {dmax} -alphas 0:0,0.4:0,0.5:1 '
              f'-out {outdir}/fdisp_{i}.png >/dev/null 2>&1')
    print(f'rendered {i}')

# horizontal colorbars for the overlay montages
os.system(f'x_showcolorbar -color jet -cmin 0 -cmax 180 -ld 30 '
          f'-unit "Strike (degree)" -lloc bottom '
          f'-out {outdir}/cb_fstrike.png >/dev/null 2>&1')
os.system(f'x_showcolorbar -color jet -cmin 0 -cmax {dmax} -ld 2 '
          f'-unit "Displacement (Grid Number)" -lloc bottom '
          f'-out {outdir}/cb_fdisp.png >/dev/null 2>&1')

letters = 'abcdefghijklmnop'
montages = {'vp': 'fig_gallery_vp.png',
            'img': 'fig_gallery_image.png',
            'fstrike': 'fig_gallery_strike.png',
            'fdisp': 'fig_gallery_disp.png'}
for kind, name in montages.items():
    fig, axes = plt.subplots(4, 4, figsize=(17, 13.5))
    for k, ax in enumerate(axes.flat):
        png = f'{outdir}/{kind}_{k + 1}.png'
        ax.imshow(mpimg.imread(png))
        ax.text(0.01, 1.02, f'({letters[k]})', transform=ax.transAxes,
                fontsize=17, va='bottom')
        ax.axis('off')
    if kind in ('fstrike', 'fdisp'):
        plt.subplots_adjust(wspace=0.04, hspace=0.10, left=0.01, right=0.99,
                            top=0.97, bottom=0.09)
        cax = fig.add_axes([0.30, 0.005, 0.40, 0.07])
        cax.imshow(mpimg.imread(f'{outdir}/cb_{kind}.png'))
        cax.axis('off')
    else:
        plt.subplots_adjust(wspace=0.04, hspace=0.10, left=0.01, right=0.99,
                            top=0.97, bottom=0.01)
    plt.savefig(f'../{name}', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'saved {name}')
