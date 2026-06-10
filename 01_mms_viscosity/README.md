# 01_mms_viscosity

This case verifies the RZ mixing-length turbulent-viscosity evaluation used by

`INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ`.

The test checks the viscosity evaluation only. It does not solve the momentum equations or verify boundary-condition treatment.

## Prescribed field

The velocity field is prescribed on

`x in [0, 1]`, `r in [0.5, 1.5]`

as

`uz(x,r) = sin(pi*x)*(1 + r^2)`

`ur(x,r) = cos(pi*x)*r^2`

with

`rho = 1`

`l_m = 0.5`

`mu = 3.333333e-3`

## Exact turbulent viscosity

The RZ strain-rate invariant is

`|S|_RZ^2 = 2*(duz/dx)^2 + 2*(dur/dr)^2 + (duz/dr + dur/dx)^2 + 2*(ur/r)^2`

and the exact turbulent viscosity is

`mu_t = rho*l_m^2*|S|_RZ`.

The effective viscosity is

`mu_eff = mu + mu_t`.

## Verification

The script compares the sampled MOOSE `mu_t` functor output against the analytical `mu_t`.

The reported table contains

`Grid, h, L2(mu_t), order, Linf(mu_t), order`.

Interior-cell norms are used to avoid boundary-reconstruction effects.

A sanity check verifies that `mu_eff = mu + mu_t` to roundoff level.

## Run

`./run_all.sh`

The paper-ready tables are

`mms_rz_viscosity_paper_table.csv`

`mms_rz_viscosity_paper_table_numeric.csv`

## Expected result

The `mu_t` error should decrease at approximately second order under mesh refinement.
