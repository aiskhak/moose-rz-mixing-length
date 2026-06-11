#!/usr/bin/env python3
from pathlib import Path
import subprocess
import argparse
import textwrap
import sympy as sp

# 02_mms_velocity
#
# Full RZ mixing-length solver MMS.
#
# This script creates:
#   n16/velocity_mms_n16.i
#   n32/velocity_mms_n32.i
#   n64/velocity_mms_n64.i
#   n128/velocity_mms_n128.i
#   n256/velocity_mms_n256.i
#   run_all_parallel.sh
#   check_convergence.py
#   manufactured_solution_functions.txt
#   README.md
#
# It verifies the complete recommended setup:
#   manufactured uz, ur, p
#   -> RZ mixing-length mu_t = rho*l_m^2*|S|_RZ
#   -> effective viscosity mu_eff = mu + mu_t
#   -> native INSFVMomentumDiffusion with mu_eff
#   -> native INSFVMomentumViscousSourceRZ with mu_eff
#   -> INSFVSymmetryVelocityBC with mu_eff
#   -> INSFVBodyForce manufactured forcing
#   -> FVPointValueConstraint pressure pin
#
# MOOSE coordinate convention in these inputs:
#   x = axial coordinate z
#   y = radial coordinate r
#   uz = axial velocity
#   ur = radial velocity

levels = {
    "n16": 16,
    "n32": 32,
    "n64": 64,
    "n128": 128,
    "n256": 256,
}

rho = 1.0
mu = 1.0e-2
lm = 0.20

x, y = sp.symbols("x y")

# Manufactured streamfunction construction.
#
# Domain:
#   x in [0,1]
#   y in [0,1], with y=0 as the RZ axis.
#
# Axisymmetric incompressible streamfunction:
#   uz = (1/y) dpsi/dy
#   ur = -(1/y) dpsi/dx
#
# We use:
#   psi = A(x) F(y)
#   F(y) = y^2 (1 - y^2)^2
#
# Then:
#   F(y)/y  = y - 2 y^3 + y^5
#   F'(y)/y = 2 - 8 y^2 + 6 y^4
#
# Therefore:
#   ur/y = -A'(x) * (1 - 2 y^2 + y^4)
#
# and at the axis:
#   ur/y -> -A'(x),
#   dur/dy -> -A'(x).
#
# This directly exercises the regularity replacement
#   ur/r -> dur/dr
# used by the RZ mixing-length functor at the RZ axis.

A = x**2 * (1 - x)**2
Ap = sp.diff(A, x)

R_over_y = y - 2*y**3 + y**5
Rp_over_y = 2 - 8*y**2 + 6*y**4

uz = A * Rp_over_y
ur = -Ap * R_over_y

p = x * y**2

uz_x = sp.diff(uz, x)
uz_y = sp.diff(uz, y)
ur_x = sp.diff(ur, x)
ur_y = sp.diff(ur, y)

# RZ strain invariant.
hoop = sp.simplify(ur / y)

s2 = (
    2 * uz_x**2
    + 2 * ur_y**2
    + (uz_y + ur_x)**2
    + 2 * hoop**2
)

mu_t = rho * lm**2 * sp.sqrt(s2)
mu_eff = mu + mu_t

# Effective viscous stresses for no-swirl axisymmetric flow.
tau_xx = 2 * mu_eff * uz_x
tau_rr = 2 * mu_eff * ur_y
tau_xr = mu_eff * (uz_y + ur_x)
tau_tt = 2 * mu_eff * hoop

# RZ divergence of viscous stress.
div_tau_x = (
    sp.diff(tau_xx, x)
    + sp.diff(tau_xr, y)
    + tau_xr / y
)

div_tau_r = (
    sp.diff(tau_xr, x)
    + sp.diff(tau_rr, y)
    + (tau_rr - tau_tt) / y
)

# Advective terms.
adv_x = rho * (uz * uz_x + ur * uz_y)
adv_r = rho * (uz * ur_x + ur * ur_y)

# Residual convention:
#   adv + grad(p) - div(tau) - forcing = 0
forcing_uz = adv_x + sp.diff(p, x) - div_tau_x
forcing_ur = adv_r + sp.diff(p, y) - div_tau_r

# Analytical checks.
div_check = sp.simplify(sp.diff(uz, x) + sp.diff(ur, y) + ur / y)
axis_limit = sp.simplify(sp.limit(ur / y, y, 0))
axis_gradient = sp.simplify(sp.limit(sp.diff(ur, y), y, 0))
axis_regular_check = sp.simplify(axis_limit - axis_gradient)

print("RZ divergence check:", div_check)
print("Axis limit ur/y:", axis_limit)
print("Axis gradient dur/dy:", axis_gradient)
print("Axis regularity check:", axis_regular_check)

pin_x = 0.5
pin_y = 0.5
pin_p = float(p.subs({x: pin_x, y: pin_y}))
print("Pressure pin value:", pin_p)


def moose_expr(expr):
    # MOOSE ParsedFunction uses ^ for powers.
    # Do not aggressively simplify the large forcing expressions.
    s = str(expr)
    s = s.replace("**", "^")
    return s


exprs = {
    "exact_uz": moose_expr(uz),
    "exact_ur": moose_expr(ur),
    "exact_p": moose_expr(p),
    "forcing_uz": moose_expr(forcing_uz),
    "forcing_ur": moose_expr(forcing_ur),
}


Path("manufactured_solution_functions.txt").write_text(
    "\n".join(
        [
            "Manufactured solution and forcing functions for 02_mms_velocity",
            "",
            "MOOSE coordinates:",
            "  x = z",
            "  y = r",
            "",
            f"RZ divergence check: {div_check}",
            f"Axis limit ur/y: {axis_limit}",
            f"Axis gradient dur/dy: {axis_gradient}",
            f"Axis regularity check: {axis_regular_check}",
            f"Pressure pin point: ({pin_x}, {pin_y})",
            f"Pressure pin value: {pin_p:.16e}",
            "",
            "[exact_uz]",
            exprs["exact_uz"],
            "",
            "[exact_ur]",
            exprs["exact_ur"],
            "",
            "[exact_p]",
            exprs["exact_p"],
            "",
            "[forcing_uz]",
            exprs["forcing_uz"],
            "",
            "[forcing_ur]",
            exprs["forcing_ur"],
            "",
        ]
    )
)


Path("README.md").write_text(
    """# 02_mms_velocity

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
"""
)


for name, n in levels.items():
    Path(name).mkdir(exist_ok=True)

    text = f"""# 02_mms_velocity / {name}
# Full RZ mixing-length solver MMS.
#
# MOOSE coordinates:
#   x = axial coordinate z
#   y = radial coordinate r

mu = {mu:.16e}
rho = {rho:.16e}
lm = {lm:.16e}

advected_interp_method = 'average'
velocity_interp_method = 'rc'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = {n}
    ny = {n}
  []
  coord_type = 'RZ'
  rz_coord_axis = x
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[GlobalParams]
  rhie_chow_user_object = rc
[]

[Variables]
  [uz]
    type = INSFVVelocityVariable
    initial_condition = 0
  []

  [ur]
    type = INSFVVelocityVariable
    initial_condition = 0
  []

  [p]
    type = INSFVPressureVariable
    initial_condition = 0
  []

  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[AuxVariables]
  [lm]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []

  [mu_t_out]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []

  [mu_eff_out]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[ICs]
  [lm_ic]
    type = ConstantIC
    variable = lm
    value = ${{lm}}
  []
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = uz
    v = ur
    pressure = p
  []
[]

[FunctorMaterials]
  [mixing_length_viscosity]
    type = INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ
    property_name = mu_eff
    turbulent_viscosity_property_name = mu_t
    molecular_viscosity = ${{mu}}
    rho = ${{rho}}
    mixing_length = lm
    u = uz
    v = ur
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = p
    rho = ${{rho}}
    advected_interp_method = ${{advected_interp_method}}
    velocity_interp_method = ${{velocity_interp_method}}
  []

  [pressure_pin]
    type = FVPointValueConstraint
    lambda = lambda
    variable = p
    point = '{pin_x} {pin_y} 0'
    phi0 = {pin_p:.16e}
  []

  [uz_advection]
    type = INSFVMomentumAdvection
    variable = uz
    advected_interp_method = ${{advected_interp_method}}
    velocity_interp_method = ${{velocity_interp_method}}
    rho = ${{rho}}
    momentum_component = 'x'
  []

  [uz_diffusion]
    type = INSFVMomentumDiffusion
    variable = uz
    mu = mu_eff
    momentum_component = 'x'
    complete_expansion = true
    u = uz
    v = ur
  []

  [uz_pressure]
    type = INSFVMomentumPressure
    variable = uz
    momentum_component = 'x'
    pressure = p
  []

  [uz_body]
    type = INSFVBodyForce
    variable = uz
    functor = forcing_uz
    momentum_component = 'x'
  []

  [ur_advection]
    type = INSFVMomentumAdvection
    variable = ur
    advected_interp_method = ${{advected_interp_method}}
    velocity_interp_method = ${{velocity_interp_method}}
    rho = ${{rho}}
    momentum_component = 'y'
  []

  [ur_diffusion]
    type = INSFVMomentumDiffusion
    variable = ur
    mu = mu_eff
    momentum_component = 'y'
    complete_expansion = true
    u = uz
    v = ur
  []

  [ur_diffusion_rz]
    type = INSFVMomentumViscousSourceRZ
    variable = ur
    mu = mu_eff
    momentum_component = 'y'
    complete_expansion = true
  []

  [ur_pressure]
    type = INSFVMomentumPressure
    variable = ur
    momentum_component = 'y'
    pressure = p
  []

  [ur_body]
    type = INSFVBodyForce
    variable = ur
    functor = forcing_ur
    momentum_component = 'y'
  []
[]

[AuxKernels]
  [mu_t_out_aux]
    type = FunctorAux
    variable = mu_t_out
    functor = mu_t
    execute_on = 'initial final'
  []

  [mu_eff_out_aux]
    type = FunctorAux
    variable = mu_eff_out
    functor = mu_eff
    execute_on = 'initial final'
  []
[]

[FVBCs]
  [axis_uz]
    type = INSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = uz
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = x
  []

  [axis_ur]
    type = INSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = ur
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = y
  []

  [axis_p]
    type = INSFVSymmetryPressureBC
    boundary = 'bottom'
    variable = p
  []

  [dirichlet_uz]
    type = INSFVNoSlipWallBC
    boundary = 'left right top'
    variable = uz
    function = exact_uz
  []

  [dirichlet_ur]
    type = INSFVNoSlipWallBC
    boundary = 'left right top'
    variable = ur
    function = exact_ur
  []
[]

[Functions]
  [exact_uz]
    type = ParsedFunction
    expression = '{exprs["exact_uz"]}'
  []

  [exact_ur]
    type = ParsedFunction
    expression = '{exprs["exact_ur"]}'
  []

  [exact_p]
    type = ParsedFunction
    expression = '{exprs["exact_p"]}'
  []

  [forcing_uz]
    type = ParsedFunction
    expression = '{exprs["forcing_uz"]}'
  []

  [forcing_ur]
    type = ParsedFunction
    expression = '{exprs["forcing_ur"]}'
  []
[]

[VectorPostprocessors]
  [vpp]
    type = ElementValueSampler
    variable = 'uz ur p mu_t_out mu_eff_out'
    sort_by = x
    execute_on = 'final'
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'lu NONZERO 200'
  []
[]

[Executioner]
  type = Steady
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 50
  l_tol = 1e-8
  l_max_its = 200
[]

[Outputs]
  print_linear_residuals = false

  [exodus]
    type = Exodus
    execute_on = FINAL
    file_base = velocity_mms_{name}_out
  []

  [csv]
    type = CSV
    execute_on = FINAL
    file_base = velocity_mms_{name}_csv
  []
[]
"""
    (Path(name) / f"velocity_mms_{name}.i").write_text(text)


Path("run_all_parallel.sh").write_text(
    """#!/bin/bash
set -e

APP=${APP:-/homes/aiskhak/projects/fv_app/fv_app-opt}

run_case () {
  level=$1
  np=$2

  echo "Running $level on $np MPI ranks"
  cd $level
  rm -f *.csv *.e *.log

  mpiexec -n $np $APP -i velocity_mms_${level}.i > run_np${np}.log 2>&1

  tail -25 run_np${np}.log
  cd ..
}

run_case n16 2
run_case n32 4
run_case n64 8
run_case n128 8
run_case n256 16

python3 check_convergence.py
"""
)
Path("run_all_parallel.sh").chmod(0o755)


Path("check_convergence.py").write_text(
    r'''#!/usr/bin/env python3
from pathlib import Path
import math
import pandas as pd
import numpy as np

levels = [
    ("n16", 16),
    ("n32", 32),
    ("n64", 64),
    ("n128", 128),
    ("n256", 256),
]

rho = 1.0
mu = 1.0e-2
lm = 0.20

xmin, xmax = 0.0, 1.0
rmin, rmax = 0.0, 1.0


def exact_fields(z, r):
    A = z**2 * (1.0 - z)**2
    Ap = 2.0 * z * (1.0 - z)**2 - 2.0 * z**2 * (1.0 - z)
    App = 2.0 - 12.0*z + 12.0*z**2

    R_over_r = r - 2.0*r**3 + r**5
    Rp_over_r = 2.0 - 8.0*r**2 + 6.0*r**4

    d_Rp_over_r_dr = -16.0*r + 24.0*r**3
    d_R_over_r_dr = 1.0 - 6.0*r**2 + 5.0*r**4

    uz = A * Rp_over_r
    ur = -Ap * R_over_r
    p = z * r**2

    uz_z = Ap * Rp_over_r
    uz_r = A * d_Rp_over_r_dr

    ur_z = -App * R_over_r
    ur_r = -Ap * d_R_over_r_dr

    # Nonsingular expression for ur/r.
    # At r = 0, this equals dur/dr = -A'(z).
    hoop = -Ap * (1.0 - 2.0*r**2 + r**4)

    s2 = (
        2.0 * uz_z**2
        + 2.0 * ur_r**2
        + (uz_r + ur_z)**2
        + 2.0 * hoop**2
    )

    mu_t = rho * lm**2 * np.sqrt(np.maximum(s2, 0.0))
    mu_eff = mu + mu_t

    return uz, ur, p, mu_t, mu_eff


def find_vpp_file(level):
    candidates = sorted(Path(level).glob(f"velocity_mms_{level}_csv_vpp*.csv"))
    if candidates:
        return candidates[-1]

    candidates = sorted(Path(level).glob("*vpp*.csv"))
    if candidates:
        return candidates[-1]

    raise FileNotFoundError(f"No VPP CSV file found in {level}")


def get_column(df, names):
    lower = {c.lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower:
            return df[lower[name.lower()]].to_numpy(dtype=float)
    raise KeyError(f"Could not find any of columns {names}. Available columns: {list(df.columns)}")


def l2(err):
    return float(np.sqrt(np.mean(err**2)))


def linf(err):
    return float(np.max(np.abs(err)))


def observed_order(e_coarse, e_fine):
    return math.log(e_coarse / e_fine) / math.log(2.0)


def offset_free_pressure_error(p_num, p_exact):
    e = p_num - p_exact
    return e - np.mean(e)


rows = []
pointwise_outputs = []

for level, n in levels:
    h = (xmax - xmin) / n
    f = find_vpp_file(level)
    df = pd.read_csv(f)

    z = get_column(df, ["x"])
    r = get_column(df, ["y"])

    uz_num = get_column(df, ["uz"])
    ur_num = get_column(df, ["ur"])
    p_num = get_column(df, ["p"])
    mu_t_num = get_column(df, ["mu_t_out"])
    mu_eff_num = get_column(df, ["mu_eff_out"])

    uz_ex, ur_ex, p_ex, mu_t_ex, mu_eff_ex = exact_fields(z, r)

    # Exclude two layers of boundary cells, consistent with the verification
    # approach used to avoid boundary reconstruction effects.
    interior = (
        (z > xmin + 2*h)
        & (z < xmax - 2*h)
        & (r > rmin + 2*h)
        & (r < rmax - 2*h)
    )

    if not np.any(interior):
        raise RuntimeError(f"No interior cells selected for {level}")

    uz_err = uz_num[interior] - uz_ex[interior]
    ur_err = ur_num[interior] - ur_ex[interior]
    p_err = offset_free_pressure_error(p_num[interior], p_ex[interior])
    mu_t_err = mu_t_num[interior] - mu_t_ex[interior]

    sanity = linf((mu_eff_num - mu_eff_ex) - (mu_t_num - mu_t_ex))

    rows.append(
        {
            "grid": f"{n}x{n}",
            "level": level,
            "n": n,
            "h": h,
            "cells": int(np.sum(interior)),
            "L2_uz": l2(uz_err),
            "L2_ur": l2(ur_err),
            "L2_p": l2(p_err),
            "L2_mu_t": l2(mu_t_err),
            "Linf_uz": linf(uz_err),
            "Linf_ur": linf(ur_err),
            "Linf_p": linf(p_err),
            "Linf_mu_t": linf(mu_t_err),
            "mu_eff_minus_mu_t_error_linf": sanity,
            "file": str(f),
        }
    )

    pointwise = pd.DataFrame(
        {
            "x": z,
            "y": r,
            "interior": interior,
            "uz": uz_num,
            "uz_exact": uz_ex,
            "uz_error": uz_num - uz_ex,
            "ur": ur_num,
            "ur_exact": ur_ex,
            "ur_error": ur_num - ur_ex,
            "p": p_num,
            "p_exact": p_ex,
            "p_error": p_num - p_ex,
            "mu_t": mu_t_num,
            "mu_t_exact": mu_t_ex,
            "mu_t_error": mu_t_num - mu_t_ex,
            "mu_eff": mu_eff_num,
            "mu_eff_exact": mu_eff_ex,
            "mu_eff_error": mu_eff_num - mu_eff_ex,
        }
    )
    pointwise_path = Path(level) / f"velocity_mms_{level}_pointwise.csv"
    pointwise.to_csv(pointwise_path, index=False)
    pointwise_outputs.append(str(pointwise_path))


for i, row in enumerate(rows):
    if i == 0:
        row["order_uz"] = "--"
        row["order_ur"] = "--"
        row["order_p"] = "--"
        row["order_mu_t"] = "--"
    else:
        prev = rows[i - 1]
        row["order_uz"] = observed_order(prev["L2_uz"], row["L2_uz"])
        row["order_ur"] = observed_order(prev["L2_ur"], row["L2_ur"])
        row["order_p"] = observed_order(prev["L2_p"], row["L2_p"])
        row["order_mu_t"] = observed_order(prev["L2_mu_t"], row["L2_mu_t"])


table = pd.DataFrame(rows)

paper = table[
    [
        "grid",
        "h",
        "L2_uz",
        "order_uz",
        "L2_ur",
        "order_ur",
        "L2_p",
        "order_p",
        "L2_mu_t",
        "order_mu_t",
    ]
].copy()

paper_display = paper.copy()

for col in ["h", "L2_uz", "L2_ur", "L2_p", "L2_mu_t"]:
    paper_display[col] = paper_display[col].map(lambda v: f"{v:.4e}")

for col in ["order_uz", "order_ur", "order_p", "order_mu_t"]:
    paper_display[col] = paper_display[col].map(
        lambda v: "--" if v == "--" else f"{float(v):.2f}"
    )

paper_display.columns = [
    "grid",
    "h",
    "L2(uz)",
    "order_L2_uz",
    "L2(ur)",
    "order_L2_ur",
    "L2(p)",
    "order_L2_p",
    "L2(mu_t)",
    "order_L2_mu_t",
]

print()
print("Table: full RZ INSFV velocity MMS with effective viscosity")
print()
print(paper_display.to_string(index=False))

sanity_max = max(row["mu_eff_minus_mu_t_error_linf"] for row in rows)

print()
print("Sanity check:")
print(f"  max |(mu_eff - mu_eff_exact) - (mu_t - mu_t_exact)| = {sanity_max:.3e}")

paper_display.to_csv("velocity_mms_paper_table.csv", index=False)
table.to_csv("velocity_mms_convergence_numeric.csv", index=False)

print()
print("Wrote:")
print("  velocity_mms_paper_table.csv")
print("  velocity_mms_convergence_numeric.csv")
for pth in pointwise_outputs:
    print(f"  {pth}")

# Simple pass/fail check.
passed = True

for key in ["L2_uz", "L2_ur", "L2_p", "L2_mu_t"]:
    vals = [row[key] for row in rows]
    if not all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)):
        passed = False
        print(f"  ERROR: {key} does not monotonically decrease")

if sanity_max > 1e-10:
    passed = False
    print("  ERROR: mu_eff consistency sanity check failed")

print()
print("Verification result:")
print("  PASS" if passed else "  FAIL")

if not passed:
    raise SystemExit(1)
'''
)
Path("check_convergence.py").chmod(0o755)


print("Wrote:")
for name in levels:
    print(f"  {name}/velocity_mms_{name}.i")
print("  run_all_parallel.sh")
print("  check_convergence.py")
print("  README.md")
print("  manufactured_solution_functions.txt")


parser = argparse.ArgumentParser()
parser.add_argument(
    "--run",
    action="store_true",
    help="After generating files, run ./run_all_parallel.sh",
)
args = parser.parse_args()

if args.run:
    subprocess.run(["./run_all_parallel.sh"], check=True)
