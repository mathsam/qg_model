# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/

CPPDEFS = -O3


include /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/lib/intel.mk


.DEFAULT:
	-echo $@ does not exist.
all: qg_run.x
eig_pak.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90
eigrg1.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f
fft_fftw.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90 par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90
io_tools.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90 par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90
numerics_lib.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90 op_rules.o io_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90
op_rules.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90
par_tools.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./mpif.h
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT).	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90
qg_arrays.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90 qg_params.o io_tools.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90
qg_diagnostics.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90 op_rules.o qg_params.o qg_arrays.o qg_strat_and_shear.o par_tools.o io_tools.o numerics_lib.o rmf_forcing.o qg_tracers.o qg_filter_tools.o transform_tools.o strat_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90
qg_driver.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90 qg_arrays.o qg_params.o qg_strat_and_shear.o qg_tracers.o qg_topo.o qg_init_streamfunction.o qg_output.o qg_diagnostics.o qg_filter_tools.o rmf_forcing.o transform_tools.o numerics_lib.o io_tools.o par_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90
qg_filter_tools.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90 io_tools.o qg_params.o qg_arrays.o par_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90
qg_init_streamfunction.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90 op_rules.o io_tools.o qg_params.o qg_arrays.o qg_strat_and_shear.o qg_diagnostics.o numerics_lib.o transform_tools.o strat_tools.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90
qg_output.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90 io_tools.o qg_params.o qg_arrays.o qg_tracers.o par_tools.o rmf_forcing.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90
qg_params.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90 io_tools.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90
qg_strat_and_shear.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90 io_tools.o qg_arrays.o qg_params.o numerics_lib.o par_tools.o strat_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90
qg_topo.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90 io_tools.o qg_params.o qg_strat_and_shear.o qg_arrays.o transform_tools.o numerics_lib.o par_tools.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90
qg_tracers.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90 op_rules.o io_tools.o qg_filter_tools.o qg_strat_and_shear.o qg_params.o par_tools.o qg_arrays.o transform_tools.o strat_tools.o numerics_lib.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90
rmf_forcing.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90 io_tools.o qg_params.o par_tools.o op_rules.o qg_arrays.o qg_filter_tools.o numerics_lib.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90
strat_tools.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90 io_tools.o numerics_lib.o eig_pak.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90
transform_tools.o: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90 fft_fftw.o par_tools.o io_tools.o op_rules.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90
./mpif.h: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./mpif.h
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./mpif.h .
./fft_fftw.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90 .
./transform_tools.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90 .
./eig_pak.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90 .
./qg_driver.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90 .
./eigrg1.f: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f .
./qg_strat_and_shear.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90 .
./qg_init_streamfunction.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90 .
./io_tools.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90 .
./rmf_forcing.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90 .
./qg_diagnostics.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90 .
./qg_arrays.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90 .
./qg_output.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90 .
./qg_topo.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90 .
./qg_params.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90 .
./numerics_lib.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90 .
./par_tools.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90 .
./qg_filter_tools.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90 .
./op_rules.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90 .
./strat_tools.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90 .
./qg_tracers.f90: /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90
	cp /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90 .
SRC = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./mpif.h
OBJ = qg_params.o par_tools.o eig_pak.o strat_tools.o qg_output.o transform_tools.o eigrg1.o qg_topo.o qg_init_streamfunction.o rmf_forcing.o io_tools.o qg_arrays.o qg_driver.o qg_strat_and_shear.o qg_diagnostics.o qg_tracers.o op_rules.o fft_fftw.o qg_filter_tools.o numerics_lib.o
OFF = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./mpif.h /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./fft_fftw.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./transform_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eig_pak.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_driver.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./eigrg1.f /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_strat_and_shear.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_init_streamfunction.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./io_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./rmf_forcing.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_diagnostics.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_arrays.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_output.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_topo.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_params.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./numerics_lib.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./par_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_filter_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./op_rules.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./strat_tools.f90 /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src/./qg_tracers.f90
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
	$(LD) $(OBJ) -o qg_run.x  $(LDFLAGS)
