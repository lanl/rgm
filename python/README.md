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
p.yn_rgt = True
p.generate()

vp = p.vp                        # numpy array, shape (n1, n2, n3)
image = p.image
rgt = p.rgt
karst = p.karst
```

Parameter names are identical to the Fortran derived-type components of
`rgm2_curved`/`rgm3_curved` (see `doc/README.md`). Scalars, booleans,
strings, and 1D numeric lists/arrays are supported; the custom input
arrays (`refl`, `refl_top`, custom `psf`) are not exposed through the
Python interface.

## Regenerating the shim after adding Fortran parameters

The C-interoperable shim (`src/module_geological_model_c.f90`) is
generated from the Fortran type definitions:

```bash
python3 python/generate_shim.py
cd src && make so
```

New parameters become available in Python automatically after this.
