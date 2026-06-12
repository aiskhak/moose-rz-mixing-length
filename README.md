# MOOSE RZ Mixing-Length Verification Cases

This repository contains source files, MOOSE input files, and verification results supporting the technical note:

**A Note on Using Axisymmetric Mixing-Length Eddy Viscosity in MOOSE**
https://doi.org/10.13140/RG.2.2.25524.41607

The repository is organized as focused verification and demonstration cases. Details for each case are provided in the corresponding folder README.

## Contents

- `INSFV/`  
  MOOSE source files for the RZ mixing-length effective-viscosity functor material.

- `01_mms_viscosity/`  
  Verification of the RZ turbulent-viscosity evaluation.

- `02_mms_velocity/`  
  Solver-level manufactured-solution verification.

- `03_tamu_demo/`  
  Reactor-relevant axisymmetric demonstration case.

- `MOOSE_VERSION.txt`  
  MOOSE and dependency versions used for the reported calculations.

## Use

Each case folder contains ready-to-run input files, run scripts, post-processing scripts, and reported tables.

The stored MOOSE input files are the reproducible inputs used for the reported results.

Edit the MOOSE executable path in the run scripts before running on a different system.

## Citation

If this repository is used, please cite the accompanying technical note. The final citation will be added after publication.
