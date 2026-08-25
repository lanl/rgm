#
# 2 x 2 galleries of 2D models with the four geomorphological unconformity
# types, from the v2.0 examples: the P-wave velocity sections and the
# corresponding synthetic seismic images.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size': 11,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

n1, n2 = 151, 301
exdir = '../tmp'
shapes = ['meander_channel', 'meander_canyon',
          'drainage_channel', 'drainage_canyon']
labels = ['(a) Meandering channels', '(b) Meandering canyon',
          '(c) Drainage channels', '(d) Drainage canyon']

clip = 0.06   # symmetric amplitude clip for the seismic images


def load(kind, shape):
    f = f'{exdir}/example_2d_{kind}_unconf_{shape}_151x301.bin'
    return np.fromfile(f, dtype=np.float32).reshape((n1, n2), order='F')


for kind, name in [('vp', 'fig_2d_sections.png'),
                   ('image', 'fig_2d_images.png')]:
    fig, axes = plt.subplots(2, 2, figsize=(13, 6.4), constrained_layout=True)
    for ax, shape, lab in zip(axes.flat, shapes, labels):
        d = load(kind, shape)
        if kind == 'vp':
            im = ax.imshow(d, cmap='jet', aspect='auto',
                           interpolation='bilinear')
            unit = '$V_p$ (m/s)'
        else:
            im = ax.imshow(d, cmap='gray_r', aspect='auto',
                           interpolation='bilinear', vmin=-clip, vmax=clip)
            unit = 'Amplitude'
        ax.set_title(lab, fontsize=12)
        ax.set_xlabel('X (Grid Number)')
        ax.set_ylabel('Z (Grid Number)')
        cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.015)
        cb.set_label(unit)
    plt.savefig('../' + name, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('saved', name)
