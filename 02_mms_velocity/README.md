# 02_mms_velocity

Full solver-level MMS verification for the final RZ mixing-length effective-viscosity setup.

The manufactured velocity field is RZ-divergence-free by construction. It also has a finite, generally nonzero axis limit for `ur/r`, so this single solver-level MMS exercises the regularity replacement `ur/r -> dur/dr` used by the RZ mixing-length functor at the RZ axis.

SymPy is used to generate the manufactured forcing terms. The generated MOOSE input files are stored directly in the grid folders and are the reproducible source for the reported calculations.

Run:

```bash
python3 make_cases.py
./run_all_parallel.sh
```

The reported convergence table is written by `check_convergence.py`.

Interior-cell norms are reported to avoid direct comparison of finite-volume boundary states, which can be affected by boundary-condition reconstruction.
