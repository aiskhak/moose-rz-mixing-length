# 03_tamu_demo

This folder contains the reactor-relevant axisymmetric demonstration case.

The case uses the final RZ mixing-length effective-viscosity setup:

`mu_t = rho*l_m^2*|S|_RZ`

`mu_eff = mu + mu_t`

The effective viscosity `mu_eff` is supplied to

`INSFVMomentumDiffusion`

`INSFVMomentumViscousSourceRZ`

`INSFVSymmetryVelocityBC`

## Case setup

The mesh is read from

`tamu.msh`

and scaled by the inlet diameter.

The inlet axial velocity uses a one-seventh power-law profile normalized by circular-area average velocity:

`uz_in = -u_avg*(60/49)*(1 - r/R)^(1/7)`

with

`u_avg = 1`

`R = 0.5`

The negative sign sets the inlet flow direction.

## Output

The input samples both

`mu_t_out`

and

`mu_eff_out`

for post-processing.

## Run

Edit the MOOSE executable path in the run script if needed, then run the case from this folder.

## Notes

This case is a demonstration of the final input-file pattern. It is not a manufactured-solution verification case.
