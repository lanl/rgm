# RGM Python interface

A thin `ctypes` wrapper around the RGM Fortran library — no algorithm is
reimplemented in Python.

## Build and install

```bash
# 1. Build the shared library (requires libflit compiled with -fPIC)
cd rgm/src
make so

# 2. Install the Python package (editable)
pip install -e ../python
```

If `librgm.so` is moved, set `RGM_LIB=/path/to/librgm.so`.

## Usage

```python
import rgm

p = rgm.rgm3(n1=128, n2=128, n3=128, nf=4, seed=1234)
p.delta_strike = [15, 25]        # strike-varying faults
p.yn_vary_disp = True            # spatially varying displacement
p.unconf = 1                     # channel-type unconformity
p.unconf_shape = 'meander_channel'
p.yn_karst = True                # karst cave system
p.nthrust = 2                    # thrust faults, older than the faults
p.yn_rgt = True
p.generate()

vp = p.vp                        # numpy array, shape (n1, n2, n3)
image = p.image
rgt = p.rgt
karst = p.karst
```

The classes `rgm.rgm2` and `rgm.rgm3` are in fact the Fortran types
`rgm2_curved` and `rgm3_curved`, not the legacy Fortran types `rgm2` and
`rgm3`, which the Python interface does not provide. Parameter names are
identical to the Fortran derived-type components of
`rgm2_curved`/`rgm3_curved` (see `doc/README.md`). Scalars, booleans,
strings, and 1D numeric lists/arrays are supported; the custom input
arrays (`refl`, `refl_top`, custom `psf`) are not exposed through the
Python interface.

## The C-interoperable shim

The C-interoperable shim (`src/module_geological_model_c.f90`) is not part of
the repository. `make so` generates it from the Fortran type definitions with
`python3 python/generate_shim.py`, and regenerates it whenever they change, so
new Fortran parameters become available in Python after `make so`. The shim is
part of `librgm.so` only; the static library `librgm.a` does not contain it.
