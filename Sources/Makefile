#
# Lists
#
FMODS= KindType.o Numeric.o InpOut.o Master.o def_elsest.o TimeFun.o elsest_mod.o Wrf.o Gfs.o\
       Debug.o Sourcemod.o EraInt.o Era5.o
FOBJS= ATLAS.o openlog.o runend.o initialize.o readinp.o set_meteo_model.o set_back_mesh.o \
       init_netcdf.o init_kml.o init_outvar.o writerest.o \
       outmeteo.o meteos.o readmeteo.o interpolate.o advection.o sedimentation.o checkdomain.o\
       make_part_nc.o reduce_mass.o outpart_nc.o setpsi.o outpart_kml.o readgrn.o grnfile.o source.o outrest.o readrestart.o\
       readnumpart.f90 readbpts.o init_random_seed.o calcdif.o diffusion.o splitting.o

.SUFFIXES:.o .f90
.f90.o:
	$(FC) -c $(FFLAGS) $(FCFLAGS) ${INC_NetCDF} $<
#
#----------------------------------------------
#  Tasks
#----------------------------------------------
#
ATLAS: $(FMODS) $(FOBJS)
	$(LINKER) -o $(EXE) $(FMODS) $(FOBJS) $(LIB_NetCDF) $(LIB_NetCDFa) -fcheck=all
	@chmod 770 $(EXE)	
	@echo '---------------------------->>> END OF COMPILATION'
#
new:
	@rm -rf *.o *.mod core*
	@make
#
clean:
	@rm -f core core.* *~ *.o *.mod
