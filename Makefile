export
#
#  1. Define the compiler (gfortran, ifort, f90, xlf90,mpif90)
#
FC=mpif90
#
#  2. Code vesion
#
ver=1.0
EXE=P_ATLAS.$(ver).x
#
#  3. Location of the NetCDF library (MN,HOME)
#
LIB_NetCDF= -L$(shell nf-config --prefix)/lib -lnetcdff -ldl -lm -Wl,-R$(shell nf-config --prefix)/lib
INC_NetCDF= $(shell nf-config --fflags) 
dir := Sources
#
#----------------------------------------
# Do not change anything beyond this line
#----------------------------------------
#
# Compiler selection
#
ifeq ($(FC),mpif90)
#  FFLAGS= -ffixed-line-length-132 -O2
  FFLAGS= -fallow-invalid-boz
  LINKER= $(FC)
  LFLAGS=
endif
ifeq ($(FC),f90)
  FFLAGS= -extend_source -r8 -64
  LINKER= $(FC) -64
  LFLAGS=
endif
ifeq ($(FC),ifort)
  FFLAGS= -132 -r8 -O1 -traceback
  LINKER= $(FC)
  LFLAGS=
endif
ifeq ($(FC),xlf90)
  LINKER= xlf90
  FOPT= -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto -q64 -WF,-D__USE_LARGEFILE64
  FFLAGS= -qfree=f90 -qrealsize=8 -qextname=flush $(FOPT)
  LINKER= $(FC)
  LFLAGS= $(FCFLAGS)
endif

all:
	@$(MAKE) -C $(dir)
	@mv $(dir)/$(EXE) Runs/

clean:
	@$(MAKE) -C $(dir) clean
	@rm -f Runs/$(EXE)
