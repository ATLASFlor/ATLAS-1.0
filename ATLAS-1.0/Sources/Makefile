#! /bin/sh -f
#
#  1. Define the compiler (gfortran, ifort, f90, xlf90)
#
FC= gfortran
#  
#  2. Code vesion
#
ver=1.0
EXE=ATLAS.$(ver).x
#
#  3. Location of the NetCDF library (MN,HOME)
#
MACHINE=HOME2
#
ifeq ($(MACHINE),HOME2)
   LIB_NetCDF= /usr/local/lib/libnetcdff.so
   INC_NetCDF= -I/usr/local/include/
endif
ifeq ($(MACHINE),HOME)
   LIB_NetCDF= /home/hp/local/lib/libnetcdff.so
   INC_NetCDF= -I/home/hp/local/include/
endif
#
#----------------------------------------
# Do not change anything beyond this line
#----------------------------------------
#
# Compiler selection
#
ifeq ($(FC),gfortran)
  FFLAGS= -ffixed-line-length-132 -O2
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
FMODS= KindType.o Numeric.o InpOut.o Master.o def_elsest.o TimeFun.o elsest_mod.o Wrf.o Gfs.o\
       Debug.o Sourcemod.o 
FOBJS= ATLAS.o openlog.o runend.o initialize.o readinp.o set_meteo_model.o set_back_mesh.o \
       outmeteo.o meteos.o readmeteo.o interpolate.o advection.o sedimentation.o checkdomain.o\
       outpart_nc.o setpsi.o outpart_kml.o readgrn.o grnfile.o source.o outrest.o readrestart.o\
       readnumpart.o readbpts.o init_random_seed.o calcdif.o diffusion.o

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
	@mv $(EXE) ../.
	@echo '---------------------------->>> END OF COMPILATION'
#
new:
	@rm -rf *.o *.mod core*
	@make
#
clean:
	@rm -f core core.* *~ *.o *.mod
