#
# Paired velocity/seismic-image figures for the four geomorphological
# unconformity types. For each type the P-wave velocity model and the
# corresponding synthetic seismic image are rendered as three-slice views
# with pymplot x_showslice at identical slice positions, then montaged
# as panels (a) and (b).
#
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.image as mpimg

mpl.rcParams.update({'font.size': 16,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

n1, n2, n3 = 151, 201, 251
exdir = '../tmp'
outdir = '../tmp'
os.makedirs(outdir, exist_ok=True)

# depth slice per type: through the interval where the incisions are
# expressed (the drainage incisions are thinner and shallower)
cases = [('meander_channel', 44), ('meander_canyon', 53),
         ('drainage_channel', 49), ('drainage_canyon', 43)]

clip = 0.06   # symmetric amplitude clip for the seismic images

for shape, s1 in cases:
    vp = f'{exdir}/example_3d_vp_unconf_{shape}_151x201x251.bin'
    im = f'{exdir}/example_3d_image_unconf_{shape}_151x201x251.bin'
    common = (f'-n1 {n1} -n2 {n2} -n3 {n3} -slice1 {s1} -slice2 100 '
              f'-slice3 125 -label1 "Z (Grid Number)" '
              f'-label2 "Y (Grid Number)" -label3 "X (Grid Number)"')
    os.system(f'x_showslice -in {vp} {common} -color jet -legend 1 '
              f'-unit "Vp (m/s)" -out {outdir}/{shape}_vp.png >/dev/null 2>&1')
    os.system(f'x_showslice -in {im} {common} -color binary -legend 1 '
              f'-unit "Amplitude" -cmin {-clip:.6f} -cmax {clip:.6f} '
              f'-out {outdir}/{shape}_img.png >/dev/null 2>&1')

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 12.6))
    for ax, kind, lab in [(axes[0], 'vp', '(a)'), (axes[1], 'img', '(b)')]:
        ax.imshow(mpimg.imread(f'{outdir}/{shape}_{kind}.png'))
        ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=18,
                va='bottom')
        ax.axis('off')
    plt.subplots_adjust(hspace=0.06, left=0.01, right=0.99,
                        top=0.97, bottom=0.01)
    plt.savefig(f'../fig_{shape}_3d.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'saved fig_{shape}_3d.png')
