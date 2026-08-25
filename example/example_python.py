#!/usr/bin/env python3
#
# © 2024-2026. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
# Security Administration. All rights in the program are reserved by
# Triad National Security, LLC, and the U.S. Department of Energy/National
# Nuclear Security Administration. The Government is granted for itself and
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display publicly,
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#
"""
Python examples for RGM v2.0.

One example per feature type of RGM v2.0 -- folded layers, faults (planar,
listric, strike-varying, slip-patch, drag-decayed), unconformities (random
topography, meandering channels/canyons, dendritic drainage networks/canyons),
salt bodies, karst cave systems, and elastic models.  Each example is
generated through the Python interface, exported to an HDF5 file, and plotted.

Calling RGM from Python is a three-liner: create an rgm2/rgm3 object, set the
parameters (identical names to the Fortran derived-type components, see
doc/README.md), and call generate().  Outputs are read back as numpy arrays:

    import rgm

    p = rgm.rgm3(n1=128, n2=160, n3=160, nf=4, seed=1234)
    p.delta_strike = [15.0, 25.0]
    p.yn_rgt = True
    p.generate()

    vp, image, rgt = p.vp, p.image, p.rgt      # numpy, shape (n1, n2, n3)

Requirements: the rgm Python package (pip install -e python/) with librgm.so
built (make -C src so), plus numpy, h5py, and matplotlib.

Usage:

    python3 example_python.py                     # all examples
    python3 example_python.py karst_3d salt_3d    # only these
    python3 example_python.py --list              # list available examples
    python3 example_python.py --outdir /tmp/rgm   # where to write (default ./rgm_python_examples)
    python3 example_python.py --no-plot           # HDF5 only

Set OMP_NUM_THREADS to control the number of threads used by the Fortran core.
"""

import argparse
import os
import sys
import time

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import rgm


# ----------------------------------------------------------------------------
# Example definitions
#
# Each example is a plain dict of RGM parameters, applied to an rgm2/rgm3
# object in order.  The names are exactly the Fortran derived-type component
# names, so any example here can be translated line by line to Fortran.
# ----------------------------------------------------------------------------

# Parameters shared by the 2D examples: a moderately folded, 25-layer medium
COMMON_2D = dict(
    n1=151,
    n2=301,
    nl=25,
    lwv=0.4,
    lwh=0.2,
    refl_shape='perlin',
    refl_shape_top='perlin',
    refl_smooth=3,
    refl_smooth_top=3,
    refl_height=[0.0, 15.0],
    refl_height_top=[0.0, 4.0],
    noise_level=0.01,
    psf_sigma=[10.0, 1.0],
)

# Parameters shared by the 3D examples
COMMON_3D = dict(
    n1=128,
    n2=160,
    n3=160,
    nl=25,
    lwv=0.4,
    lwh=0.2,
    refl_shape='perlin',
    refl_shape_top='perlin',
    refl_smooth=3,
    refl_smooth_top=3,
    refl_height=[0.0, 15.0],
    refl_height_top=[0.0, 4.0],
    noise_level=0.01,
    psf_sigma=[10.0, 1.0, 1.0],
)


def _case(dim, doc, base, x1=0.5, **params):
    # x1 is the relative depth of the map-view slice used when plotting a 3D
    # model; it only affects the figure, not the model, and is set per example
    # so that the slice cuts through the feature the example demonstrates
    return dict(dim=dim, doc=doc, x1=x1, params=dict(base, **params))


EXAMPLES = {

    # ---- layers and folds --------------------------------------------------

    'layers_2d': _case(
        2, 'Folded layers only (no faults), with RGT and facies labels',
        COMMON_2D,
        nf=0,                   # nf = 0 is what removes the faults; the
                                # yn_fault flag only controls whether the
                                # fault label arrays are output
        refl_height=[0.0, 25.0],
        yn_rgt=True,
        yn_facies=True,
        seed=101,
    ),

    'fold_3d': _case(
        3, 'Gaussian anticline/syncline folds (refl_shape = gaussian)',
        COMMON_3D,
        nf=0,
        refl_shape='gaussian',
        ng=3,
        rotate_fold=True,
        refl_sigma2=[30.0, 60.0],
        refl_sigma3=[30.0, 60.0],
        refl_height=[0.0, 50.0],
        refl_smooth=0,
        yn_rgt=True,
        seed=102,
    ),

    # ---- faults ------------------------------------------------------------

    'fault_2d': _case(
        2, 'Listric normal/reverse faults in 2D (dip, disp, delta_dip)',
        COMMON_2D,
        nf=4,
        dip=[60.0, 120.0],
        disp=[5.0, 12.0],
        delta_dip=[0.0, 20.0],
        yn_rgt=True,
        seed=201,
    ),

    'fault_3d': _case(
        3, 'Faults in 3D with prescribed dip, strike, and rake',
        COMMON_3D,
        nf=5,
        dip=[60.0, 120.0],
        strike=[20.0, 60.0],
        rake=[0.0, 30.0],
        disp=[8.0, 16.0],
        delta_dip=[0.0, 15.0],
        yn_rgt=True,
        seed=202,
    ),

    'fault_strike_3d': _case(
        3, 'Strike-varying (curved in map view) faults: delta_strike > 0 [v2.0]',
        COMMON_3D,
        nf=4,
        dip=[60.0, 120.0],
        disp=[10.0, 20.0],
        delta_dip=[0.0, 15.0],
        delta_strike=[15.0, 25.0],
        strike_nperiod=2,
        yn_rgt=True,
        seed=203,
    ),

    'fault_vary_disp_3d': _case(
        3, 'Elliptical slip patch: displacement dies out at the fault tips [v2.0]',
        COMMON_3D,
        nf=4,
        dip=[60.0, 120.0],
        disp=[10.0, 20.0],
        delta_strike=[10.0, 20.0],
        yn_vary_disp=True,
        disp_radius_strike=[0.4, 0.6],
        disp_radius_dip=[0.5, 0.8],
        disp_center_dip=[0.3, 0.6],
        yn_rgt=True,
        seed=204,
    ),

    'fault_decay_3d': _case(
        3, 'Displacement decay away from the fault: drag folds/rollover [v2.0]',
        COMMON_3D,
        nf=3,
        dip=[60.0, 120.0],
        disp=[12.0, 22.0],
        delta_strike=[10.0, 20.0],
        yn_vary_disp=True,
        disp_radius_strike=[0.4, 0.6],
        disp_radius_dip=[0.5, 0.8],
        yn_disp_decay=True,
        disp_decay_width=[0.3, 0.5],
        yn_rgt=True,
        seed=205,
    ),

    # ---- unconformities ----------------------------------------------------

    'unconf_2d': _case(
        2, 'Two unconformities with random (Perlin) erosional topography',
        COMMON_2D,
        nf=3,
        disp=[5.0, 10.0],
        unconf=2,
        unconf_z=[0.15, 0.4],
        unconf_height=[5.0, 15.0],
        unconf_nl=10,
        yn_rgt=True,
        seed=301,
    ),

    'meander_channel_3d': _case(
        3, 'Unconformity carved by meandering river channels [v2.0]',
        COMMON_3D,
        nf=2,
        disp=[5.0, 10.0],
        unconf=1,
        unconf_z=[0.2, 0.3],
        unconf_shape='meander_channel',
        x1=0.30,
        unconf_channel_width=[0.04, 0.08],
        unconf_channel_sinuosity=1.2,
        unconf_topo=0.25,
        yn_rgt=True,
        seed=302,
    ),

    'meander_canyon_3d': _case(
        3, 'Unconformity carved by a meandering incised canyon [v2.0]',
        COMMON_3D,
        nf=2,
        disp=[5.0, 10.0],
        unconf=1,
        unconf_z=[0.2, 0.3],
        unconf_shape='meander_canyon',
        x1=0.30,
        unconf_channel_width=[0.05, 0.10],
        unconf_height=[10.0, 20.0],
        yn_rgt=True,
        seed=303,
    ),

    'drainage_channel_3d': _case(
        3, 'Unconformity carved by a dendritic drainage network [v2.0]',
        COMMON_3D,
        nf=2,
        disp=[5.0, 10.0],
        unconf=1,
        unconf_z=[0.2, 0.3],
        unconf_shape='drainage_channel',
        x1=0.235,
        unconf_channel_density=[0.03, 0.08],
        yn_rgt=True,
        seed=304,
    ),

    'drainage_canyon_3d': _case(
        3, 'Unconformity carved by a dendritic drainage canyon system [v2.0]',
        COMMON_3D,
        nf=2,
        disp=[5.0, 10.0],
        unconf=1,
        unconf_z=[0.2, 0.3],
        unconf_shape='drainage_canyon',
        x1=0.36,
        unconf_channel_density=[0.03, 0.08],
        unconf_height=[10.0, 20.0],
        yn_rgt=True,
        seed=305,
    ),

    # ---- salt --------------------------------------------------------------

    'salt_3d': _case(
        3, 'Salt bodies (mask in %salt), inserted into the medium models',
        COMMON_3D,
        nf=2,
        disp=[5.0, 10.0],
        yn_salt=True,
        nsalt=2,
        salt_radius=[20.0, 30.0],
        salt_top_z=[0.4, 0.6],
        salt_nnode=8,
        salt_path_variation=6.0,
        seed=401,
    ),

    # ---- karst -------------------------------------------------------------

    'karst_3d': _case(
        3, 'Karst cave system: connected tube network (mask in %karst) [v2.0]',
        COMMON_3D,
        nf=3,
        dip=[60.0, 120.0],
        disp=[5.0, 10.0],
        yn_karst=True,
        karst_z=[0.45, 0.85],
        karst_npassage=25,
        karst_nctrl=18,
        karst_connect=0.4,
        karst_tortuosity=0.6,
        yn_rgt=True,
        seed=501,
    ),

    # ---- elastic -----------------------------------------------------------

    'elastic_3d': _case(
        3, 'Elastic model: Vp, Vs, density, and the four PP/PS/SP/SS images',
        COMMON_3D,
        nf=3,
        dip=[60.0, 120.0],
        disp=[8.0, 16.0],
        delta_strike=[10.0, 20.0],
        yn_elastic=True,
        vpvsratio=[1.6, 1.9],
        seed=601,
    ),
}


# ----------------------------------------------------------------------------
# Output arrays
# ----------------------------------------------------------------------------

# All output arrays RGM can produce; those not enabled for a given example are
# simply absent (the interface raises AttributeError, and we skip them).
OUTPUTS = [
    'vp', 'vs', 'rho',
    'image', 'image_pp', 'image_ps', 'image_sp', 'image_ss',
    'rgt', 'facies',
    'fault', 'fault_dip', 'fault_strike', 'fault_rake', 'fault_disp',
    'salt', 'karst',
]

# Plotting styles: seismic images are symmetric gray, labels that use zero as
# background are masked so the background stays white, masks are drawn on top
# of Vp.
SEISMIC = {'image', 'image_pp', 'image_ps', 'image_sp', 'image_ss'}
OVERLAY = {'fault', 'salt', 'karst'}
LABELS = {'fault_dip', 'fault_strike', 'fault_rake', 'fault_disp'}

TITLES = {
    'vp': 'Vp (m/s)',
    'vs': 'Vs (m/s)',
    'rho': 'Density (kg/m^3)',
    'image': 'Seismic image',
    'image_pp': 'PP image',
    'image_ps': 'PS image',
    'image_sp': 'SP image',
    'image_ss': 'SS image',
    'rgt': 'Relative geologic time',
    'facies': 'Facies',
    'fault': 'Fault index',
    'fault_dip': 'Fault dip (deg)',
    'fault_strike': 'Fault strike (deg)',
    'fault_rake': 'Fault rake (deg)',
    'fault_disp': 'Fault displacement (grid)',
    'salt': 'Salt mask',
    'karst': 'Karst mask',
}


def generate(name):
    """Build, generate, and read back one example. Returns (arrays, params)."""

    spec = EXAMPLES[name]
    params = spec['params']

    p = rgm.rgm3() if spec['dim'] == 3 else rgm.rgm2()
    for key, value in params.items():
        setattr(p, key, value)

    p.generate()

    arrays = {}
    for out in OUTPUTS:
        try:
            arrays[out] = getattr(p, out)
        except AttributeError:
            pass        # not enabled for this example

    p.free()
    return arrays, params


def export_hdf5(path, name, spec, arrays, params):
    """Write all output arrays, plus the parameters used, to an HDF5 file."""

    with h5py.File(path, 'w') as f:
        f.attrs['name'] = name
        f.attrs['description'] = spec['doc']
        f.attrs['ndim'] = spec['dim']
        f.attrs['generator'] = 'rgm2_curved' if spec['dim'] == 2 else 'rgm3_curved'
        # the parameters are stored so that an example is reproducible from
        # the file alone
        for key, value in params.items():
            f.attrs['param_' + key] = value
        for key, value in arrays.items():
            d = f.create_dataset(key, data=value, compression='gzip',
                                 compression_opts=4)
            d.attrs['long_name'] = TITLES.get(key, key)


def _slices(a, ndim, x1=0.5):
    """Panels to draw: the array itself in 2D, three orthogonal sections in 3D."""

    if ndim == 2:
        return [(a, '')]
    n1, n2, n3 = a.shape
    k = min(int(x1*n1), n1 - 1)
    return [(a[:, :, n3//2], 'x3 = %d' % (n3//2)),
            (a[:, n2//2, :], 'x2 = %d' % (n2//2)),
            (a[k, :, :], 'x1 = %d (map view)' % k)]


def _draw(ax, name, panel, background):
    """Draw one panel and return the mappable used for the colorbar."""

    if name in SEISMIC:
        c = np.percentile(np.abs(panel), 99.0)
        c = c if c > 0 else 1.0
        return ax.imshow(panel, cmap='gray', vmin=-c, vmax=c, aspect='auto')

    if name in OVERLAY:
        if background is not None:
            ax.imshow(background, cmap='binary', aspect='auto')
        m = np.ma.masked_where(panel == 0, panel)
        return ax.imshow(m, cmap='autumn', aspect='auto')

    if name in LABELS:
        m = np.ma.masked_where(panel == 0, panel)
        return ax.imshow(m, cmap='jet', aspect='auto')

    cmap = 'nipy_spectral' if name == 'facies' else 'jet'
    return ax.imshow(panel, cmap=cmap, aspect='auto')


def plot(path, name, spec, arrays):
    """One figure per example: one row per output array."""

    ndim = spec['dim']
    names = [k for k in OUTPUTS if k in arrays]
    ncol = 3 if ndim == 3 else 1
    nrow = len(names)

    fig, axes = plt.subplots(nrow, ncol, squeeze=False, layout='constrained',
                             figsize=(4.6*ncol + 1.2, 2.5*nrow))
    fig.suptitle('%s -- %s' % (name, spec['doc']), fontsize=11)

    x1 = spec['x1']
    vp = arrays.get('vp')
    backgrounds = _slices(vp, ndim, x1) if vp is not None else None
    for i, key in enumerate(names):
        im = None
        for j, (panel, title) in enumerate(_slices(arrays[key], ndim, x1)):
            ax = axes[i][j]
            bg = backgrounds[j][0] if backgrounds is not None else None
            im = _draw(ax, key, panel, bg)
            if i == 0:
                ax.set_title(title, fontsize=9)
            ax.tick_params(labelsize=7)
        axes[i][0].set_ylabel(TITLES.get(key, key), fontsize=9)
        if key not in OVERLAY:      # a 0/1 mask needs no colorbar
            cb = fig.colorbar(im, ax=list(axes[i]), fraction=0.03, pad=0.01)
            cb.ax.tick_params(labelsize=7)

    fig.savefig(path, dpi=110)
    plt.close(fig)


def main(argv=None):

    parser = argparse.ArgumentParser(
        description='Generate the RGM v2.0 feature examples through the '
                    'Python interface, export them as HDF5, and plot them.')
    parser.add_argument('cases', nargs='*',
                        help='examples to run (default: all)')
    parser.add_argument('--list', action='store_true',
                        help='list the available examples and exit')
    parser.add_argument('--outdir', default='./rgm_python_examples',
                        help='output directory (default: %(default)s)')
    parser.add_argument('--no-plot', action='store_true',
                        help='export HDF5 only, do not plot')
    args = parser.parse_args(argv)

    if args.list:
        for name, spec in EXAMPLES.items():
            print('%-22s %dD  %s' % (name, spec['dim'], spec['doc']))
        return 0

    cases = args.cases or list(EXAMPLES)
    unknown = [c for c in cases if c not in EXAMPLES]
    if unknown:
        parser.error('unknown example(s): %s (use --list)' % ', '.join(unknown))

    os.makedirs(args.outdir, exist_ok=True)

    for name in cases:
        spec = EXAMPLES[name]
        t = time.time()
        arrays, params = generate(name)

        h5 = os.path.join(args.outdir, name + '.h5')
        export_hdf5(h5, name, spec, arrays, params)

        png = ''
        if not args.no_plot:
            png = os.path.join(args.outdir, name + '.png')
            plot(png, name, spec, arrays)

        print('%-22s %5.1f s  %dD  %s  ->  %s %s'
              % (name, time.time() - t, spec['dim'],
                 ', '.join(sorted(arrays)), h5, png))
        sys.stdout.flush()

    print('\nDone. Output written to %s' % os.path.abspath(args.outdir))
    return 0


if __name__ == '__main__':
    sys.exit(main())
