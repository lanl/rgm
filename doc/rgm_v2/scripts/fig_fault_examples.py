#
# The 3D fault example figures: strike-varying faults, spatially varying
# displacement, and displacement decay, plus the corresponding P-wave
# velocity models. Fault attributes are overlaid on the synthetic seismic
# image; the velocity models are montaged as panels (a) and (b).
#
# Data are the vstrike_1 / vdisp_1 / decay_1 blocks of
# example/example_rgm3_curved.f90.
#
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.image as mpimg

mpl.rcParams.update({'font.size': 16,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

n1, n2, n3 = 171, 251, 251
exdir = '../tmp'
outdir = '../tmp'
os.makedirs(outdir, exist_ok=True)

common = (f'-n1 {n1} -n2 {n2} -n3 {n3} -slice1 40 -slice2 125 -slice3 125 '
          f'-label1 "Z (Grid Number)" -label2 "Y (Grid Number)" '
          f'-label3 "X (Grid Number)"')

# Each figure pairs the velocity model (a) with the fault attributes
# overlaid on the synthetic seismic image (b).
cases = [
    ('fig_strike.png', 'vstrike_1', 'fstrike_vstrike_1', 'image_vstrike_1',
     'Fault Strike (degree)', '0:0,0.9:0,1:1'),
    ('fig_disp.png', 'vdisp_1', 'fdisp_vdisp_1', 'image_vdisp_1',
     'Fault Displacement (Grid Number)', '0:0,0.49:0,0.5:1'),
    ('fig_decay.png', 'decay_1', 'fdisp_decay_1', 'image_decay_1',
     'Fault Displacement (Grid Number)', '0:0,0.49:0,0.5:1'),
]

for out, tag, attr, back, unit, alphas in cases:
    os.system(f'x_showslice -in {exdir}/example_3d_vp_{tag}.bin {common} '
              f'-color jet -legend 1 -unit "Vp (m/s)" '
              f'-out {outdir}/{tag}_vp.png >/dev/null 2>&1')
    os.system(f'x_showslice -in {exdir}/example_3d_{attr}.bin '
              f'-background {exdir}/example_3d_{back}.bin {common} '
              f'-color jet -backcolor gray -alphas {alphas} '
              f'-legend 1 -unit "{unit}" '
              f'-out {outdir}/{tag}_img.png >/dev/null 2>&1')

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 12.6))
    for ax, kind, lab in [(axes[0], 'vp', '(a)'), (axes[1], 'img', '(b)')]:
        ax.imshow(mpimg.imread(f'{outdir}/{tag}_{kind}.png'))
        ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=18,
                va='bottom')
        ax.axis('off')
    plt.subplots_adjust(hspace=0.06, left=0.01, right=0.99,
                        top=0.97, bottom=0.01)
    plt.savefig('../' + out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('saved', out)
