# INSFV source files

This folder contains the MOOSE source files for the RZ mixing-length effective-viscosity functor material:

- `INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.h`
- `INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.C`

The material computes the RZ turbulent dynamic viscosity

`mu_t = rho*l_m^2*|S|_RZ`

and the effective dynamic viscosity

`mu_eff = mu + mu_t`.

The intended input-file pattern is to use `mu_eff` in:

- `INSFVMomentumDiffusion`
- `INSFVMomentumViscousSourceRZ`
- `INSFVSymmetryVelocityBC`

The turbulent viscosity `mu_t` can be sampled separately for verification and post-processing.

## Use

Copy the source files into the corresponding MOOSE application or module source/include directories and rebuild the application.

The calculations in this repository used the executable and MOOSE version listed in `../MOOSE_VERSION.txt`.

## Scope

This implementation is intended for the INSFV RZ use pattern demonstrated in this repository. It is not presented here as a general-purpose replacement for all finite-element, nodal, or non-RZ functor-material evaluations.
