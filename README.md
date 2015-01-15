
12/20/07: With J Marshall, discovered that tracer equation is stirred
only by initial (target) mean velocity.  Whe time-varying mean state
is allowed, tracer must also be stirred by adjusted mean.

So, as a cluge the current version (.43) includes time varying mean in
tracer, but cannot be run at higher vertical res in tracer (cant have
nzt?nz)

ALSO corrected:

1.  on restart, reads in current ubar and vbar, but reads in original
input file for u_target and v_target (just realized that if original
ubar is set in model, then ubar_in_file is not set - should change
this).

2.  now updates toposhift at each timestep if topo is used and mean is stepped

3.  includes hb in calculation of vq and uq:  (q -> q+hb)

THIS version includes forcing of the mean zonal shear via

dU/dt = <vq> - a*(U - Uo)

New parameter for namelist:  inv_restore_time (= a)

If this is nonzero, Program will read in or create ubar, then read
that value into ubar_target, which will stay constant, then step ubar
forward via leapfrog with same timestep as full model.  This is all done in
new subprogram in main qg_driver, called Step_U

A logical flag (not a namelist parameter) called time_varying_mean is
set to false by default, then set to true in check_parameters, at
model start, if inv_restore_time/=0


VERSION 4.41 does the same for V:

dV/dt = -<uq> - a*(V - Vo)



Some recent changes made to this version of code:

1.  Added default true flag 'normalize_ubar' that toggles removal of BT mode and normalization
of mean velocities by their maxima.

2.  Added vertical tracer diffusivity, with parameter 'kappa_v'.  Required alteration of
qg_tracers.f90 and numerics_lib.f90.

Version with these changes called sqgp4.1


11/14/06:  Added straight horizontal laplacian diffusion to tracer, with one new optional
parameter:  kappa_h

5/11/07: Version 4.32 has slightly cleaned up implementation of
thermal drag. No conditioning on vaue of Fe, which didn't make sense
anyway.

Version 4.43:  Includes:

diffz in numerics_lib
Get_corrs in qg_diagnostics
do_corrs in qg_params
