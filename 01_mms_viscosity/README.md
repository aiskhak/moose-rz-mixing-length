01_mms_viscosity
This case verifies the RZ mixing-length turbulent-viscosity evaluation used by the final effective-viscosity functor material:
INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ

The test is a grid-convergence check for the viscosity evaluation only. It does not solve the momentum equations and does not verify boundary-condition treatment.

*Prescribed velocity field*

The velocity field is prescribed analytically on the RZ domain:
x in [0, 1], r in [0.5, 1.5]
with
uz(x,r) = sin(pix)(1 + r^2)
ur(x,r) = cos(pi*x)*r^2

The constants are:
rho = 1
l_m = 0.5
mu = 3.333333e-3

*Exact RZ strain invariant*

The exact derivatives are:
duz/dx = picos(pix)*(1 + r^2)
duz/dr = 2rsin(pi*x)
dur/dx = -pisin(pix)*r^2
dur/dr = 2rcos(pi*x)
ur/r = rcos(pix)

The RZ strain-rate invariant is:
|S|_RZ^2 =
2*(duz/dx)^2
+ 2*(dur/dr)^2

(duz/dr + dur/dx)^2
+ 2*(ur/r)^2

The exact turbulent viscosity is:
mu_t = rhol_m^2|S|_RZ

The exact effective viscosity is:
mu_eff = mu + mu_t

*What is verified*

The primary verification compares the MOOSE sampled solve-path turbulent-viscosity functor output against the analytical value of mu_t.
The script reports the paper-ready table:
Grid, h, L2(mu_t), order, Linf(mu_t), order
Interior-cell norms are used as the primary results to avoid boundary-reconstruction effects.
As a sanity check, the script also verifies that the effective-viscosity output satisfies:
mu_eff = mu + mu_t
to roundoff level.

*How to run*

python3 make_cases.py
./run_all.sh

The convergence table is written to:
mms_rz_viscosity_paper_table.csv
mms_rz_viscosity_paper_table_numeric.csv

*Expected behavior*

The turbulent-viscosity error should decrease at approximately second order under mesh refinement. The effective-viscosity sanity check should remain near roundoff.