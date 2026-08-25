#
# Three-panel map-view figure of the meander-channel erosional depth map
# for unconf_channel_length = 12000, 25000 (default), and 100000,
# demonstrating the calibrated scale-safe length mapping.
#
# Input bins are produced by rgm/tmp/test_meander_length.f90 with the
# current library (seed 12345, width fraction 0.05).
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size': 11,
                     'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans',
                                         'DejaVu Sans']})

ny, nx = 201, 251
lens = [12000, 25000, 100000]
labels = ['(a) $L = 12000$', '(b) $L = 25000$ (default)', '(c) $L = 100000$']

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)

for ax, L, lab in zip(axes, lens, labels):
    d = np.fromfile(f'../tmp/mlen_channel_{L}.bin',
                    dtype=np.float32).reshape((ny, nx), order='F')
    im = ax.imshow(d, cmap='jet', aspect='equal', interpolation='bilinear')
    ax.set_title(lab, fontsize=12)
    ax.set_xlabel('X (Grid Number)')
    ax.set_ylabel('Y (Grid Number)')

cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.015)
cb.set_label('Erosional Depth (Grid Number)')

plt.savefig('../fig_channel_length.png', dpi=200, bbox_inches='tight')
print('saved fig_channel_length.png')
