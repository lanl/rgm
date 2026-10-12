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


# paths
bindir = $(PWD)/../lib
objdir = $(PWD)/../obj
moddir = $(PWD)/../lib

# Dependencies
common = \
    module_geological_model_utility.o \
    module_geological_model_realism.o \
    module_geological_model_thrust.o \
    module_geological_model_2d.o \
    module_geological_model_3d.o \
    module_geological_model_meander.o \
    module_geological_model_drainage.o \
    module_geological_model_karst.o \
    module_geological_model_2d_curved.o \
    module_geological_model_3d_curved.o \
    module_geological_model.o

objlib = $(addprefix $(objdir)/, $(common))

# targets
lib = $(moddir)/librgm.a

all: makedir $(lib)

# options
flitdir = $(HOME)/src/flit

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)
fflags = $(base_fflags) \
	# -O0 -g -check bounds

# compile
$(objdir)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link; the archive is rebuilt from scratch, so that it holds only these objects
$(lib): $(objlib)
	-@rm -f $@
	$(ll) rcs $@ $^

# shared library for the Python interface, which also holds the C-interoperable
# shim; the shim is regenerated from the derived types whenever they change
solib = $(moddir)/librgm.so
shim = $(PWD)/module_geological_model_c.f90

so: makedir $(solib)

$(solib): $(objlib) $(objdir)/module_geological_model_c.o
	$(fc) -shared -o $@ $^ $(lflags)

$(shim): $(PWD)/module_geological_model_2d_curved.f90 $(PWD)/module_geological_model_3d_curved.f90 \
		$(PWD)/../python/shim_template.f90 $(PWD)/../python/generate_shim.py
	python3 $(PWD)/../python/generate_shim.py

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod $(lib) $(shim)
