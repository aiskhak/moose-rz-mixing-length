#!/usr/bin/env python3
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
