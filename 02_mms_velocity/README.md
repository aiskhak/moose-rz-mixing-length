# 02_mms_velocity

Purpose: full solver-level MMS verification for the final RZ mixing-length effective-viscosity setup.

This verifies the complete recommended path:
- RZ mixing-length turbulent viscosity
- effective viscosity mu_eff = mu + mu_t
- native INSFVMomentumDiffusion with mu_eff
- native INSFVMomentumViscousSourceRZ with mu_eff
- INSFVSymmetryVelocityBC with mu_eff
- manufactured momentum forcing through INSFVBodyForce
- pressure pinning through FVPointValueConstraint
- axisymmetric finite-volume setup

The exact velocity field is RZ-divergence-free by construction.
