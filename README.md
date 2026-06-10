# MOOSE RZ Mixing-Length Verification Cases

This repository contains verification and demonstration files supporting the technical note:

**A Note on Using Axisymmetric Mixing-Length Eddy Viscosity in MOOSE**

The cases support an RZ-consistent mixing-length eddy-viscosity implementation for MOOSE incompressible finite-volume Navier--Stokes calculations.

## Purpose

The repository documents and verifies the use of an axisymmetric mixing-length turbulent viscosity in MOOSE finite-volume RANS calculations. The implemented formulation computes the RZ-consistent turbulent viscosity

`mu_t = rho*l_m^2*|S|_RZ`

and the effective viscosity

`mu_eff = mu + mu_t`.

The effective viscosity is supplied consistently to the native MOOSE finite-volume viscous momentum operators, the RZ radial viscous source term, and the symmetry velocity boundary condition.

## Repository contents

### `01_mms_viscosity/`

Grid-convergence verification of the RZ mixing-length turbulent-viscosity evaluation.

This case prescribes a smooth manufactured velocity field and compares the sampled MOOSE functor output for `mu_t` against the analytical value. It verifies the RZ strain-rate invariant used in the closure, including the hoop-strain contribution.

Primary reported output:

* `mms_rz_viscosity_paper_table.csv`
* `mms_rz_viscosity_paper_table_numeric.csv`

### `02_mms_velocity/`

Full solver-level manufactured-solution verification using the final effective-viscosity setup.

This case verifies the complete recommended solve path:

* RZ mixing-length turbulent viscosity
* effective viscosity `mu_eff = mu + mu_t`
* native `INSFVMomentumDiffusion` with `mu_eff`
* native `INSFVMomentumViscousSourceRZ` with `mu_eff`
* `INSFVSymmetryVelocityBC` with `mu_eff`
* manufactured forcing through `INSFVBodyForce`
* pressure pinning through `FVPointValueConstraint`

Primary reported output:

* `velocity_mms_paper_table.csv`
* `velocity_mms_convergence_numeric.csv`

### `tamu/`

Reactor-relevant axisymmetric demonstration case based on the final recommended input-file pattern.

This case demonstrates practical use of the RZ mixing-length effective-viscosity formulation in an axisymmetric thermal-hydraulics configuration.

## Implementation files

The repository includes the MOOSE object used in the verification and demonstration cases:

* `INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.h`
* `INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.C`

The material exposes two functor properties:

* `mu_t`: turbulent dynamic viscosity
* `mu_eff`: effective dynamic viscosity

The intended input-file pattern is to use `mu_eff` in all viscous momentum operators and boundary viscous treatments, while sampling `mu_t` for verification and post-processing.

## Running the cases

The repository contains generated input files and run/post-processing scripts. The case-generation scripts are not included.

Before running, edit the MOOSE executable path in the run scripts if needed.

Example:

`APP=/path/to/fv_app-opt`

Run the viscosity verification:

`cd 01_mms_viscosity`

`./run_all.sh`

Run the full solver MMS verification:

`cd 02_mms_velocity`

`./run_all_parallel.sh`

The scripts execute the existing input files and then run the convergence post-processing scripts.

## Expected verification behavior

For `01_mms_viscosity`, the turbulent-viscosity error should decrease at approximately second order under mesh refinement.

For `02_mms_velocity`, the velocity, pressure, and turbulent-viscosity errors should decrease under mesh refinement, with approximately second-order behavior for the velocity fields.

## Notes

The verification cases are intended to support reproducibility of the paper results. They are not a general-purpose turbulence-model validation suite. The manufactured-solution cases verify implementation consistency and convergence behavior; the reactor-relevant case demonstrates practical use of the final RZ input-file pattern.
