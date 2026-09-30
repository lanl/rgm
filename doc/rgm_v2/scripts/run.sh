#!/bin/bash
#
# Reproduce all data and figures of the RGM 2.0 manuscript.
#
#   bash run.sh data      generate all model data into ../tmp (Fortran)
#   bash run.sh figures   render all figures into ../ (Python)
#   bash run.sh           both stages
#
# Requirements: librgm.a built in ../../lib (make -C ../../src), libflit,
# the Intel toolchain used by the package Makefiles, and a Python with
# numpy, scipy, matplotlib, and pyvista (plus pymplot's x_showslice,
# x_showmatrix, and x_showcolorbar on the PATH).
#
# All generators and figure scripts use fixed random seeds, so rerunning
# this script reproduces the figures of the manuscript exactly.
#
set -e

SCRIPTS=$(cd "$(dirname "$0")" && pwd)
ROOT=${RGMROOT:-$(cd "$SCRIPTS/../../.." && pwd)}
TMP="$SCRIPTS/../tmp"
mkdir -p "$TMP"

PYTHON=${PYTHON:-python3}
FC=${FC:-mpiifx}
FLITDIR=${FLITDIR:-$HOME/src/flit}

FFLAGS="-O3 -fpp -qopenmp -heap-arrays -xHost -fp-model=consistent
        -fpscomp logicals -I$FLITDIR/lib -I$ROOT/lib"
LIBS="$ROOT/lib/librgm.a $FLITDIR/lib/libflit.a -qmkl -lstdc++ -lstdc++fs
      -L$HOME/intel/lib -L$HOME/intel/mkl/lib
      $HOME/intel/mkl/lib/libmkl_core.a
      $HOME/intel/mkl/lib/libmkl_intel_thread.a
      $HOME/intel/mkl/lib/libmkl_intel_lp64.a
      $HOME/intel/mkl/lib/libmkl_blas95_lp64.a
      $HOME/intel/mkl/lib/libmkl_lapack95_lp64.a
      -L$HOME/intel/mpi/lib -L$HOME/intel/mpi/lib/release
      -L$HOME/tool/hdf5/lib -lhdf5_fortran -lhdf5 -lz -ldl"

GENERATORS="gen_fault_examples gen_meander_length gen_extras
            gen_unconf_models gen_karst_models gen_gallery"

FIGURES="fig_fault_schematic fig_disp_schematic fig_disp_3d
         fig_geomorph_schematic fig_channel_length fig_geomorph_topo
         fig_fault_3d fig_fault_examples fig_karst_iso
         fig_karst_iso_connected fig_unconf_pairs fig_karst_pairs
         fig_2d_sections fig_gallery"

stage_data() {
    for prog in $GENERATORS; do
        echo "=== $prog ==="
        $FC -o "$TMP/$prog" $FFLAGS "$SCRIPTS/$prog.f90" $LIBS
        (cd "$TMP" && "./$prog")
    done
}

stage_figures() {
    for fig in $FIGURES; do
        echo "=== $fig ==="
        (cd "$SCRIPTS" && "$PYTHON" "$fig.py")
    done
}

case "${1:-all}" in
    data)    stage_data ;;
    figures) stage_figures ;;
    all)     stage_data; stage_figures ;;
    *)       echo "usage: bash run.sh [data|figures|all]"; exit 1 ;;
esac
echo "DONE"
