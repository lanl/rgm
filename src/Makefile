
# paths
bindir = $(PWD)/../bin
objdir = $(HOME)/bin/obj/geological_model
moddir = $(PWD)/../lib

# Dependencies
common = \
    module_geological_model_utility.o \
    module_geological_model_2d.o \
    module_geological_model_3d.o \
    module_geological_model_2d_elastic.o \
    module_geological_model_3d_elastic.o \
    module_geological_model.o

objlib = $(addprefix $(objdir)/, $(common))

# targets
lib = $(moddir)/librgm.a

all: makedir \
	$(lib)

# options
include $(HOME)/src/__libflit/src/Makefile.in
inc = \
	$(base_inc) \
	-I$(HOME)/src/__libflit/lib \
	-I$(HOME)/src/__geological_model/lib
lflags = \
	$(HOME)/src/__geological_model/lib/librgm.a \
	$(HOME)/src/__libflit/lib/libflit.a \
	$(base_lflags)
fflags = $(base_fflags) \
	#	-O0 -g -check bounds

# compile
$(objdir)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(lib): $(objlib)
	$(ll) rcs $@ $^

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod \
		$(lib)
