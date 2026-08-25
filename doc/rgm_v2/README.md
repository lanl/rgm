# Figures of the RGM 2.0 paper (LA-UR-26-27587)

This directory holds the figures of the RGM 2.0 manuscript together with
everything needed to regenerate them. All generators and figure scripts use
fixed random seeds, so a rerun reproduces the published figures exactly.

```
rgm_v2/
  fig_*.png     the figures, as used in the manuscript
  scripts/      Fortran model generators (gen_*.f90) and
                Python figure scripts (fig_*.py), driven by run.sh
  tmp/          created by run.sh; intermediate models and images
```

## Requirements

* `librgm.a` built in the package `lib/` directory (`make -C src`), plus
  [FLIT](https://github.com/lanl/flit) and the Intel toolchain
  (`mpiifx`, MKL, HDF5) used by the package Makefiles.
* Python 3 with `numpy`, `scipy`, `matplotlib`, and `pyvista`.
* [pymplot](https://github.com/lanl/pymplot) — `x_showslice`, `x_showmatrix`,
  and `x_showcolorbar` must be on the `PATH` (used by `fig_gallery.py`).

## Reproducing everything

```bash
cd doc/rgm_v2/scripts
bash run.sh            # data + figures
bash run.sh data       # only the Fortran stage, writes ../tmp
bash run.sh figures    # only the Python stage, writes ../fig_*.png
```

The Fortran stage compiles and runs the six generators into `../tmp`; the
Python stage reads those binaries and writes each `fig_*.png` one level up,
overwriting the copies shipped here. Expect the data stage to take a while
and to produce a few GB in `tmp/` — the 3D gallery models dominate both.

Paths and tools are picked up from the environment, so override them if your
setup differs:

```bash
PYTHON=/path/to/python3 FC=ifx FLITDIR=$HOME/src/flit \
RGMROOT=/path/to/rgm bash run.sh
```

## Rebuilding a single figure

Each figure script is standalone. Run it from `scripts/` after the generator
that produces its input has run:

```bash
cd doc/rgm_v2/scripts
python3 fig_karst_pairs.py     # -> ../fig_karst_3d.png, ../fig_karst_2d.png
```

| Figure | Script | Data generator |
| --- | --- | --- |
| `fig_fault_schematic.png` | `fig_fault_schematic.py` | none (analytic) |
| `fig_disp_schematic.png` | `fig_disp_schematic.py` | none (analytic) |
| `fig_disp_3d.png` | `fig_disp_3d.py` | none (analytic) |
| `fig_geomorph_schematic.png` | `fig_geomorph_schematic.py` | none (analytic) |
| `fig_strike.png`, `fig_disp.png`, `fig_decay.png` | `fig_fault_examples.py` | `gen_fault_examples.f90` |
| `fig_fault_3d.png` | `fig_fault_3d.py` | `gen_fault_examples.f90` |
| `fig_channel_length.png` | `fig_channel_length.py` | `gen_meander_length.f90` |
| `fig_geomorph_topo.png` | `fig_geomorph_topo.py` | `gen_extras.f90` |
| `fig_karst_iso_connected.png` | `fig_karst_iso_connected.py` | `gen_extras.f90` |
| `fig_meander_channel_3d.png`, `fig_meander_canyon_3d.png`, `fig_drainage_channel_3d.png`, `fig_drainage_canyon_3d.png` | `fig_unconf_pairs.py` | `gen_unconf_models.f90` |
| `fig_2d_sections.png`, `fig_2d_images.png` | `fig_2d_sections.py` | `gen_unconf_models.f90` |
| `fig_karst_3d.png`, `fig_karst_2d.png` | `fig_karst_pairs.py` | `gen_karst_models.f90` |
| `fig_karst_iso.png` | `fig_karst_iso.py` | `gen_karst_models.f90` |
| `fig_gallery_vp.png`, `fig_gallery_image.png`, `fig_gallery_strike.png`, `fig_gallery_disp.png` | `fig_gallery.py` | `gen_gallery.f90` |
