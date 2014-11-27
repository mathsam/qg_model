#!/bin/csh

set NUMPRO = 1
set i_run = 0

module unload PrgEnv-pgi
module load PrgEnv-intel

cd ../run

time -p aprun -n $NUMPRO ./qg_run.x | tee ${i_run}.${$}.log || exit 1
#aprun -n $NUMPRO ./qg_run.x | tee ${i_run}.${$}.log || exit 1
