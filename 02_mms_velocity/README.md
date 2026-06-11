# 02_mms_velocity

This case verifies the full RZ INSFV solver setup using a manufactured solution. The test exercises the final RZ mixing-length effective-viscosity path:

`mu_t = rho*l_m^2*|S|_RZ`

`mu_eff = mu + mu_t`

The effective viscosity `mu_eff` is supplied to

`INSFVMomentumDiffusion`

`INSFVMomentumViscousSourceRZ`

`INSFVSymmetryVelocityBC`

The case also uses manufactured body forces, pressure coupling, Rhie--Chow interpolation, and a pressure pin through `FVPointValueConstraint`.

## Manufactured field

The MOOSE coordinate convention is

`x = z`

`y = r`

with `y = 0` as the RZ axis.

The velocity field is constructed from an axisymmetric streamfunction so that the RZ continuity equation is satisfied analytically. Let

`A(x) = x^2*(1 - x)^2`

The manufactured fields are

`uz(x,r) = A(x)*(2 - 8*r^2 + 6*r^4)`

`ur(x,r) = -A'(x)*(r - 2*r^3 + r^5)`

`p(x,r) = x*r^2`

This choice gives

`ur/r = -A'(x)*(1 - 2*r^2 + r^4)`

so that

`ur/r -> -A'(x) = dur/dr`

at the RZ axis. Therefore, this single solver-level MMS also exercises the regularity replacement `ur/r -> dur/dr` used by the RZ mixing-length functor.

## Verification

SymPy is used to generate the manufactured forcing terms. The generated MOOSE input files are stored directly in the grid folders and are the reproducible source for the reported calculations.

The script compares the MOOSE solution against the analytical `uz`, `ur`, `p`, and `mu_t` values. Interior-cell norms are used to avoid direct comparison of finite-volume boundary states, which can be affected by boundary-condition reconstruction.

A sanity check verifies that

`mu_eff = mu + mu_t`

to roundoff level.

## Run

Generate the cases and run the verification with

```bash
python3 make_cases.py
./run_all_parallel.sh
```

or run both steps with

```bash
python3 make_cases.py --run
```

The paper-ready tables are

`velocity_mms_paper_table.csv`

`velocity_mms_convergence_numeric.csv`

The script also writes pointwise diagnostic files in each grid folder.

## Expected result

The errors in `uz`, `ur`, `p`, and `mu_t` should decrease monotonically under mesh refinement. On the refined grids, the observed orders should approach near-second-order behavior for the full nonlinear RZ finite-volume solve.
