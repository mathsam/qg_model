# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/

CPPDEFS = -O3


include /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/lib/intel.mk


.DEFAULT:
	-echo $@ does not exist.
all: qg_run.x
eig_pak.o: $(SRCROOT)eig_pak.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)eig_pak.f90
eigrg1.o: $(SRCROOT)eigrg1.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)eigrg1.f
fft_fftw.o: $(SRCROOT)fft_fftw.f90 par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fft_fftw.f90
io_tools.o: $(SRCROOT)io_tools.f90 par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)io_tools.f90
numerics_lib.o: $(SRCROOT)numerics_lib.f90 op_rules.o io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)numerics_lib.f90
op_rules.o: $(SRCROOT)op_rules.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)op_rules.f90
par_tools.o: $(SRCROOT)par_tools.f90  $(MPICH_DIR)/include/mpif.h 
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(MPICH_DIR)/include/mpif.h	$(SRCROOT)par_tools.f90
qg_arrays.o: $(SRCROOT)qg_arrays.f90 qg_params.o io_tools.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_arrays.f90
qg_diagnostics.o: $(SRCROOT)qg_diagnostics.f90 op_rules.o qg_params.o qg_arrays.o qg_strat_and_shear.o par_tools.o io_tools.o numerics_lib.o rmf_forcing.o qg_tracers.o qg_filter_tools.o transform_tools.o strat_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_diagnostics.f90
qg_driver.o: $(SRCROOT)qg_driver.f90 qg_arrays.o qg_params.o qg_strat_and_shear.o qg_tracers.o qg_topo.o qg_init_streamfunction.o qg_output.o qg_diagnostics.o qg_filter_tools.o rmf_forcing.o transform_tools.o numerics_lib.o io_tools.o par_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_driver.f90
qg_filter_tools.o: $(SRCROOT)qg_filter_tools.f90 io_tools.o qg_params.o qg_arrays.o par_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_filter_tools.f90
qg_init_streamfunction.o: $(SRCROOT)qg_init_streamfunction.f90 op_rules.o io_tools.o qg_params.o qg_arrays.o qg_strat_and_shear.o qg_diagnostics.o numerics_lib.o transform_tools.o strat_tools.o par_tools.o nc_io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_init_streamfunction.f90
qg_output.o: $(SRCROOT)qg_output.f90 io_tools.o qg_params.o qg_arrays.o qg_tracers.o par_tools.o rmf_forcing.o nc_io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_output.f90
nc_io_tools.o: $(SRCROOT)nc_io_tools.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c  -I$(NETCDF_DIR)/include      $(SRCROOT)nc_io_tools.f90
qg_params.o: $(SRCROOT)qg_params.f90 io_tools.o par_tools.o nc_io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_params.f90
qg_strat_and_shear.o: $(SRCROOT)qg_strat_and_shear.f90 io_tools.o qg_arrays.o qg_params.o numerics_lib.o par_tools.o strat_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_strat_and_shear.f90
qg_topo.o: $(SRCROOT)qg_topo.f90 io_tools.o qg_params.o qg_strat_and_shear.o qg_arrays.o transform_tools.o numerics_lib.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_topo.f90
qg_tracers.o: $(SRCROOT)qg_tracers.f90 op_rules.o io_tools.o qg_filter_tools.o qg_strat_and_shear.o qg_params.o par_tools.o qg_arrays.o transform_tools.o strat_tools.o numerics_lib.o nc_io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)qg_tracers.f90
rmf_forcing.o: $(SRCROOT)rmf_forcing.f90 io_tools.o qg_params.o par_tools.o op_rules.o qg_arrays.o qg_filter_tools.o numerics_lib.o nc_io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)rmf_forcing.f90
strat_tools.o: $(SRCROOT)strat_tools.f90 io_tools.o numerics_lib.o eig_pak.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)strat_tools.f90
transform_tools.o: $(SRCROOT)transform_tools.f90 fft_fftw.o par_tools.o io_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)transform_tools.f90
./nc_io_tools.f90: $(SRCROOT)nc_io_tools.f90
	cp $(SRCROOT)nc_io_tools.f90
./mpif.h: $(SRCROOT)mpif.h
	cp $(SRCROOT)mpif.h .
./fft_fftw.f90: $(SRCROOT)fft_fftw.f90
	cp $(SRCROOT)fft_fftw.f90 .
./transform_tools.f90: $(SRCROOT)transform_tools.f90
	cp $(SRCROOT)transform_tools.f90 .
./eig_pak.f90: $(SRCROOT)eig_pak.f90
	cp $(SRCROOT)eig_pak.f90 .
./qg_driver.f90: $(SRCROOT)qg_driver.f90
	cp $(SRCROOT)qg_driver.f90 .
./eigrg1.f: $(SRCROOT)eigrg1.f
	cp $(SRCROOT)eigrg1.f .
./qg_strat_and_shear.f90: $(SRCROOT)qg_strat_and_shear.f90
	cp $(SRCROOT)qg_strat_and_shear.f90 .
./qg_init_streamfunction.f90: $(SRCROOT)qg_init_streamfunction.f90
	cp $(SRCROOT)qg_init_streamfunction.f90 .
./io_tools.f90: $(SRCROOT)io_tools.f90
	cp $(SRCROOT)io_tools.f90 .
./rmf_forcing.f90: $(SRCROOT)rmf_forcing.f90
	cp $(SRCROOT)rmf_forcing.f90 .
./qg_diagnostics.f90: $(SRCROOT)qg_diagnostics.f90
	cp $(SRCROOT)qg_diagnostics.f90 .
./qg_arrays.f90: $(SRCROOT)qg_arrays.f90
	cp $(SRCROOT)qg_arrays.f90 .
./qg_output.f90: $(SRCROOT)qg_output.f90
	cp $(SRCROOT)qg_output.f90 .
./qg_topo.f90: $(SRCROOT)qg_topo.f90
	cp $(SRCROOT)qg_topo.f90 .
./qg_params.f90: $(SRCROOT)qg_params.f90
	cp $(SRCROOT)qg_params.f90 .
./numerics_lib.f90: $(SRCROOT)numerics_lib.f90
	cp $(SRCROOT)numerics_lib.f90 .
./par_tools.f90: $(SRCROOT)par_tools.f90
	cp $(SRCROOT)par_tools.f90 .
./qg_filter_tools.f90: $(SRCROOT)qg_filter_tools.f90
	cp $(SRCROOT)qg_filter_tools.f90 .
./op_rules.f90: $(SRCROOT)op_rules.f90
	cp $(SRCROOT)op_rules.f90 .
./strat_tools.f90: $(SRCROOT)strat_tools.f90
	cp $(SRCROOT)strat_tools.f90 .
./qg_tracers.f90: $(SRCROOT)qg_tracers.f90
	cp $(SRCROOT)qg_tracers.f90 .
SRC = $(SRCROOT)qg_params.f90 $(SRCROOT)par_tools.f90 $(SRCROOT)eig_pak.f90 $(SRCROOT)strat_tools.f90 $(SRCROOT)qg_output.f90 $(SRCROOT)transform_tools.f90 $(SRCROOT)eigrg1.f $(SRCROOT)qg_topo.f90 $(SRCROOT)qg_init_streamfunction.f90 $(SRCROOT)rmf_forcing.f90 $(SRCROOT)io_tools.f90 $(SRCROOT)qg_arrays.f90 $(SRCROOT)qg_driver.f90 $(SRCROOT)qg_strat_and_shear.f90 $(SRCROOT)qg_diagnostics.f90 $(SRCROOT)qg_tracers.f90 $(SRCROOT)op_rules.f90 $(SRCROOT)fft_fftw.f90 $(SRCROOT)qg_filter_tools.f90 $(SRCROOT)numerics_lib.f90 $(SRCROOT)mpif.h $(SRCROOT)nc_io_tools.f90
OBJ = qg_params.o par_tools.o eig_pak.o strat_tools.o qg_output.o transform_tools.o eigrg1.o qg_topo.o qg_init_streamfunction.o rmf_forcing.o io_tools.o qg_arrays.o qg_driver.o qg_strat_and_shear.o qg_diagnostics.o qg_tracers.o op_rules.o fft_fftw.o qg_filter_tools.o numerics_lib.o nc_io_tools.o
OFF = $(SRCROOT)mpif.h $(SRCROOT)fft_fftw.f90 $(SRCROOT)transform_tools.f90 $(SRCROOT)eig_pak.f90 $(SRCROOT)qg_driver.f90 $(SRCROOT)eigrg1.f $(SRCROOT)qg_strat_and_shear.f90 $(SRCROOT)qg_init_streamfunction.f90 $(SRCROOT)io_tools.f90 $(SRCROOT)rmf_forcing.f90 $(SRCROOT)qg_diagnostics.f90 $(SRCROOT)qg_arrays.f90 $(SRCROOT)qg_output.f90 $(SRCROOT)qg_topo.f90 $(SRCROOT)qg_params.f90 $(SRCROOT)numerics_lib.f90 $(SRCROOT)par_tools.f90 $(SRCROOT)qg_filter_tools.f90 $(SRCROOT)op_rules.f90 $(SRCROOT)strat_tools.f90 $(SRCROOT)qg_tracers.f90
clean: neat
	-rm -f .qg_run.x.cppdefs $(OBJ) qg_run.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
qg_run.x: $(OBJ) 
	$(LD) -ipo $(OBJ) -o qg_run.x  $(LDFLAGS)
