#
# © 2024 - 2025. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
# Security Administration. All rights in the program are reserved by
# Triad National Security, LLC, and the U.S. Department of Energy/National
# Nuclear Security Administration. The Government is granted for itself and
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare. derivative works,
# distribute copies to the public, perform publicly and display publicly,
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#


# paths
bindir = $(PWD)/bin
objdir = $(PWD)/obj
moddir = $(PWD)/obj

# Dependencies
object1 = example_mtl_2d.o
object2 = example_mtl_3d.o
object3 = example_rgm2_curved.o
object4 = example_rgm3_curved.o
object5 = example_fault.o

obj1 = $(addprefix $(objdir)/, $(object1))
obj2 = $(addprefix $(objdir)/, $(object2))
obj3 = $(addprefix $(objdir)/, $(object3))
obj4 = $(addprefix $(objdir)/, $(object4))
obj5 = $(addprefix $(objdir)/, $(object5))

# targets
exec1 = $(bindir)/x_generate_mtl_2d
exec2 = $(bindir)/x_generate_mtl_3d
exec3 = $(bindir)/x_generate_curved_2d
exec4 = $(bindir)/x_generate_curved_3d
exec5 = $(bindir)/x_generate_fault

all: makedir \
	$(exec1) $(exec2) \
	$(exec3) $(exec4) \
	$(exec5)

# options
flitdir = $(HOME)/src/flit
rgmdir = $(HOME)/src/rgm

include $(flitdir)/src/Makefile.in
inc = \
	$(base_inc) \
	-I$(flitdir)/lib \
	-I$(rgmdir)/lib
lflags = \
	$(rgmdir)/lib/librgm.a \
	$(flitdir)/lib/libflit.a \
	$(base_lflags)
fflags = $(base_fflags) \
	#	-O0 -g -check bounds

# compile
$(objdir)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec1) : $(obj1)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec2) : $(obj2)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec3) : $(obj3)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec4) : $(obj4)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec5) : $(obj5)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod \
		$(lib) \
		$(exec1) $(exec2) \
		$(exec3) $(exec4) \
		$(exec5)
