#!/bin/csh -e

set echo

set cppDefs = (-O3)
set root_dir = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/
set src_dir = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/src
set bld_dir = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/bld
set run_dir = /lustre/f1/Junyi.Chai/testqg 
set template_dir = /lustre/f1/unswept/Junyi.Chai/qg_model_mpi/lib

mkdir -p $bld_dir $run_dir

module load fre

cd $src_dir

\rm -f pathnames_solo
list_paths -o pathnames_solo

cd $bld_dir

mkmf -a $src_dir -t $template_dir/intel.mk -c "$cppDefs" -p qg_run.x $src_dir/pathnames_solo

module unload PrgEnv-pgi
module load PrgEnv-intel
module load netcdf

make DEBUG=ON -j 32

cp -f $bld_dir/qg_run.x $run_dir/qg_run.x

