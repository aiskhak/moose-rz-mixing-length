import numpy as np
import pandas as pd
from pathlib import Path

levels = [
    ("n16", 16),
    ("n32", 32),
    ("n64", 64),
    ("n128", 128),
]

rho = 1.0
mu = 1.0e-2
lm = 0.20

xmin, xmax = 0.0, 1.0
rmin, rmax = 0.0, 1.0


def exact_fields(x, r):
    A = x**2 * (1.0 - x)**2
    Ap = 2.0*x*(1.0 - x)**2 - 2.0*x**2*(1.0 - x)

    R_over_r = r**2 - 2.0*r**3 + r**4
    Rp_over_r = 3.0*r - 8.0*r**2 + 5.0*r**3

    uz = A * Rp_over_r
    ur = -Ap * R_over_r
    p = x * r**2

    return uz, ur, p


def exact_mu_t(x, r):
    A = x**2 * (1.0 - x)**2
    Ap = 2.0*x*(1.0 - x)**2 - 2.0*x**2*(1.0 - x)
    App = 2.0 - 12.0*x + 12.0*x**2

    R_over_r = r**2 - 2.0*r**3 + r**4
    Rp_over_r = 3.0*r - 8.0*r**2 + 5.0*r**3

    d_Rp_over_r_dr = 3.0 - 16.0*r + 15.0*r**2
    d_R_over_r_dr = 2.0*r - 6.0*r**2 + 4.0*r**3

    uz_x = Ap * Rp_over_r
    uz_r = A * d_Rp_over_r_dr

    ur_x = -App * R_over_r
    ur_r = -Ap * d_R_over_r_dr

    # This manufactured field is regular at r=0.
    # ur/r = -Ap*(r - 2r^2 + r^3), including the r=0 limit.
    hoop = -Ap * (r - 2.0*r**2 + r**3)

    s2 = (
        2.0 * uz_x**2
        + 2.0 * ur_r**2
        + (uz_r + ur_x)**2
        + 2.0 * hoop**2
    )

    return rho * lm**2 * np.sqrt(s2)


def l2_linf(error):
    return np.sqrt(np.mean(error**2)), np.max(np.abs(error))


def offset_free_pressure_error(p_num, p_exact):
    # Pressure is pinned, but this makes the norm robust to tiny constant offsets.
    offset = np.mean(p_num - p_exact)
    return (p_num - offset) - p_exact


rows = []

for name, n in levels:
    files = sorted(Path(name).glob("*vpp*.csv"))
    if not files:
        raise RuntimeError(f"No VPP file found for {name}")

    df = pd.read_csv(files[-1])
    h = (xmax - xmin) / n

    # Interior-cell norms: exclude boundary/axis-adjacent cells.
    interior = (
        (df["x"] > xmin + 2*h)
        & (df["x"] < xmax - 2*h)
        & (df["y"] > rmin + 2*h)
        & (df["y"] < rmax - 2*h)
    )

    d = df[interior].copy()
    xq = d["x"].to_numpy()
    rq = d["y"].to_numpy()

    uz_exact, ur_exact, p_exact = exact_fields(xq, rq)
    mu_t_exact = exact_mu_t(xq, rq)
    mu_eff_exact = mu + mu_t_exact

    e_uz = d["uz"].to_numpy() - uz_exact
    e_ur = d["ur"].to_numpy() - ur_exact
    e_p = offset_free_pressure_error(d["p"].to_numpy(), p_exact)
    e_mu_t = d["mu_t_out"].to_numpy() - mu_t_exact
    e_mu_eff = d["mu_eff_out"].to_numpy() - mu_eff_exact

    uz_l2, uz_linf = l2_linf(e_uz)
    ur_l2, ur_linf = l2_linf(e_ur)
    p_l2, p_linf = l2_linf(e_p)
    mu_t_l2, mu_t_linf = l2_linf(e_mu_t)

    rows.append({
        "grid": f"{n}x{n}",
        "n": n,
        "h": h,
        "cells": len(d),

        "L2(uz)": uz_l2,
        "order_L2_uz": np.nan,

        "L2(ur)": ur_l2,
        "order_L2_ur": np.nan,

        "L2(p)": p_l2,
        "order_L2_p": np.nan,

        "L2(mu_t)": mu_t_l2,
        "order_L2_mu_t": np.nan,

        "Linf(uz)": uz_linf,
        "Linf(ur)": ur_linf,
        "Linf(p)": p_linf,
        "Linf(mu_t)": mu_t_linf,

        "mu_eff_check": np.max(np.abs(e_mu_eff - e_mu_t)),
        "file": str(files[-1]),
    })

res = pd.DataFrame(rows)

for i in range(1, len(res)):
    h1, h2 = res.loc[i - 1, "h"], res.loc[i, "h"]

    for field, order_col in [
        ("L2(uz)", "order_L2_uz"),
        ("L2(ur)", "order_L2_ur"),
        ("L2(p)", "order_L2_p"),
        ("L2(mu_t)", "order_L2_mu_t"),
    ]:
        e1, e2 = res.loc[i - 1, field], res.loc[i, field]
        res.loc[i, order_col] = np.log(e1 / e2) / np.log(h1 / h2)

res.to_csv("velocity_mms_convergence_numeric.csv", index=False)

paper = res[
    [
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
].copy()

for col in ["h", "L2(uz)", "L2(ur)", "L2(p)", "L2(mu_t)"]:
    paper[col] = paper[col].map(lambda x: f"{x:.4e}")

for col in ["order_L2_uz", "order_L2_ur", "order_L2_p", "order_L2_mu_t"]:
    paper[col] = paper[col].map(lambda x: "--" if pd.isna(x) else f"{x:.2f}")

paper.to_csv("velocity_mms_paper_table.csv", index=False)

print()
print("Table: full RZ INSFV velocity MMS with effective viscosity")
print()
print(paper.to_string(index=False))
print()
print("Sanity check:")
print(f"  max |(mu_eff - mu_eff_exact) - (mu_t - mu_t_exact)| = {res['mu_eff_check'].max():.3e}")
print()
print("Wrote:")
print("  velocity_mms_paper_table.csv")
print("  velocity_mms_convergence_numeric.csv")
print()

ok = (
    res.loc[len(res)-1, "L2(uz)"] < res.loc[0, "L2(uz)"]
    and res.loc[len(res)-1, "L2(ur)"] < res.loc[0, "L2(ur)"]
    and res.loc[len(res)-1, "L2(p)"] < res.loc[0, "L2(p)"]
    and res.loc[len(res)-1, "L2(mu_t)"] < res.loc[0, "L2(mu_t)"]
    and res["mu_eff_check"].max() < 1.0e-10
)

print("Verification result:")
if ok:
    print("  PASS")
else:
    print("  NOT PASS")
