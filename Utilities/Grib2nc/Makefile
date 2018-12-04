#! /bin/sh -f
#
#  1. Define the compiler (gfortran, ifort, f90, xlf90)
#
FC= gfortran
#  
#  2. Code vesion
#
ver=1.0
EXE=Grib2nc
#
#  3. Location of the NetCDF library (MN,HOME)
#
MACHINE=HOME
#
ifeq ($(MACHINE),HOME)
   LIB_NetCDF= /usr/local/lib/libnetcdff.so
   INC_NetCDF= -I/usr/local/include/
endif
ifeq ($(MACHINE),ARNAU)
   LIB_NetCDF= /Users/arnaufolch/Documents/Software/Libraries/netcdf/lib/libnetcdf.a
   INC_NetCDF= -I/Users/arnaufolch/Documents/Software/Libraries/netcdf/include
endif
ifeq ($(MACHINE),MN)
   LIB_NetCDF= /apps/NETCDF/3.6.2/lib/libnetcdf.a
   INC_NetCDF= -I/apps/NETCDF/3.6.2/include/
endif
#
#----------------------------------------
# Do not change anything beyond this line
#----------------------------------------
#
# Compiler selection
#
ifeq ($(FC),gfortran)
  FFLAGS= -ffixed-line-length-132  -fdefault-real-8 -O3
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
  LFLAGS= $(FFLAGS) 
endif
#
# Lists
#
FMODS= KindType.o InpOut.o Master.o 
FOBJS= Grib2nc.o inival.o openinp.o runend.o write_nc_grid.o write_nc_var2d.o write_nc_var3d.o

.SUFFIXES:.o .f90
.f90.o:
	$(FC) -c $(FFLAGS) ${INC_NetCDF} $<
#
#----------------------------------------------
#  Tasks
#----------------------------------------------
#
ATLAS: $(FMODS) $(FOBJS)
	$(LINKER) -o $(EXE) $(FMODS) $(FOBJS) $(LIB_NetCDF) $(LIB_NetCDFa)
	@chmod 770 $(EXE)
	@echo '---------------------------->>> END OF COMPILATION'
#
new:
	@rm -rf *.o *.mod core*
	@make
#
clean:
	@rm -f core core.* *~ *.o *.mod
