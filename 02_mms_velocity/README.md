# 02_mms_velocity

This case verifies the full RZ finite-volume solver setup using the final effective-viscosity formulation.

The test uses a manufactured solution for velocity and pressure. The exact velocity field is RZ-divergence-free by construction.

## Verified setup

The input uses

`mu_t = rho*l_m^2*|S|_RZ`

and

`mu_eff = mu + mu_t`.

The effective viscosity `mu_eff` is supplied to

`INSFVMomentumDiffusion`

`INSFVMomentumViscousSourceRZ`

`INSFVSymmetryVelocityBC`

The manufactured forcing is applied through

`INSFVBodyForce`.

The pressure is pinned using

`FVPointValueConstraint`.

## Verification

The script compares the computed solution against the analytical manufactured fields.

The reported table contains

`Grid, h, L2(uz), order, L2(ur), order, L2(p), order, L2(mu_t), order`.

Interior-cell norms are used to avoid boundary-reconstruction effects.

A sanity check verifies that `mu_eff = mu + mu_t` to roundoff level.

## Run

`./run_all_parallel.sh`

The paper-ready tables are

`velocity_mms_paper_table.csv`

`velocity_mms_convergence_numeric.csv`

## Expected result

The velocity, pressure, and turbulent-viscosity errors should decrease under mesh refinement. The velocity errors should approach second-order convergence.
