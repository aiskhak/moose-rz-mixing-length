#02_mms_velocity

This case verifies the full RZ finite-volume solver setup using the final effective-viscosity formulation.

The test uses a manufactured solution for velocity and pressure. The velocity field is constructed from an axisymmetric streamfunction, so the RZ incompressibility constraint is satisfied exactly.

##Manufactured solution

The RZ domain is

x in [0, 1], r in [0, 1]

with r = 0 as the symmetry axis.

Here, psi is the axisymmetric streamfunction. It is used only to construct a smooth manufactured velocity field that satisfies RZ continuity:

uz = (1/r)*dpsi/dr

ur = -(1/r)*dpsi/dx

The streamfunction is written in nonsingular polynomial form using

A(x) = x^2*(1 - x)^2

R(r)/r = r^2 - 2*r^3 + r^4

(dR/dr)/r = 3*r - 8*r^2 + 5*r^3

Therefore,

uz(x,r) = A(x)*(3*r - 8*r^2 + 5*r^3)

ur(x,r) = -dA/dx*(r^2 - 2*r^3 + r^4)

The manufactured pressure is

p(x,r) = x*r^2

The pressure is pinned at

(x,r) = (0.5, 0.5).

Effective viscosity

The input uses

mu_t = rho*l_m^2*|S|_RZ

and

mu_eff = mu + mu_t.

The constants are

rho = 1

mu = 1.0e-2

l_m = 0.20

Manufactured forcing

The analytical forcing terms were generated with SymPy from the manufactured velocity, pressure, and effective-viscosity fields.

The symbolic generation script is not included in this repository. The resulting forcing expressions are embedded directly in the generated MOOSE input files in each grid folder:

n16/velocity_mms_n16.i

n32/velocity_mms_n32.i

n64/velocity_mms_n64.i

n128/velocity_mms_n128.i

These input files are the reproducible source for the MMS calculations reported here.

Verified setup

The effective viscosity mu_eff is supplied to

INSFVMomentumDiffusion

INSFVMomentumViscousSourceRZ

INSFVSymmetryVelocityBC

The manufactured forcing is applied through

INSFVBodyForce.

The pressure is pinned using

FVPointValueConstraint.

Verification

The script compares the computed solution against the analytical manufactured fields.

The reported table contains

Grid, h, L2(uz), order, L2(ur), order, L2(p), order, L2(mu_t), order.

Interior-cell norms are used to avoid boundary-reconstruction effects.

A sanity check verifies that mu_eff = mu + mu_t to roundoff level.

Run

./run_all_parallel.sh

The paper-ready tables are

velocity_mms_paper_table.csv

velocity_mms_convergence_numeric.csv

Expected result

The velocity, pressure, and turbulent-viscosity errors should decrease under mesh refinement. The velocity errors should approach second-order convergence.
