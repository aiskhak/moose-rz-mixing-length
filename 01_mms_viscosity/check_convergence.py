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
lm = 0.5
mu = 3.333333e-3

xmin, xmax = 0.0, 1.0
rmin, rmax = 0.5, 1.5


def exact_mu_t(x, r):
    duz_dx = np.pi * np.cos(np.pi * x) * (1.0 + r**2)
    duz_dr = 2.0 * r * np.sin(np.pi * x)

    dur_dx = -np.pi * np.sin(np.pi * x) * r**2
    dur_dr = 2.0 * r * np.cos(np.pi * x)

    hoop = r * np.cos(np.pi * x)

    s2 = (
        2.0 * duz_dx**2
        + 2.0 * dur_dr**2
        + (duz_dr + dur_dx)**2
        + 2.0 * hoop**2
    )

    return rho * lm**2 * np.sqrt(s2)


def norm_report(error):
    return np.sqrt(np.mean(error**2)), np.max(np.abs(error))


rows = []

for name, n in levels:
    files = sorted(Path(name).glob("*vpp*.csv"))
    if not files:
        raise RuntimeError(f"No VPP file found for {name}")

    df = pd.read_csv(files[-1])
    h = (xmax - xmin) / n

    interior = (
        (df["x"] > xmin + 2*h)
        & (df["x"] < xmax - 2*h)
        & (df["y"] > rmin + 2*h)
        & (df["y"] < rmax - 2*h)
    )

    d = df[interior].copy()

    mu_t_exact = exact_mu_t(d["x"].to_numpy(), d["y"].to_numpy())
    mu_eff_exact = mu + mu_t_exact

    mu_t_error = d["mu_t_out"].to_numpy() - mu_t_exact
    mu_eff_error = d["mu_eff_out"].to_numpy() - mu_eff_exact

    mu_t_l2, mu_t_linf = norm_report(mu_t_error)

    rows.append({
        "grid": f"{n}x{n}",
        "n": n,
        "h": h,
        "cells": len(d),
        "L2(mu_t)": mu_t_l2,
        "order_L2": np.nan,
        "Linf(mu_t)": mu_t_linf,
        "order_Linf": np.nan,
        "mu_eff_check": np.max(np.abs(mu_eff_error - mu_t_error)),
        "file": str(files[-1]),
    })

res = pd.DataFrame(rows)

for i in range(1, len(res)):
    h1, h2 = res.loc[i - 1, "h"], res.loc[i, "h"]

    e1, e2 = res.loc[i - 1, "L2(mu_t)"], res.loc[i, "L2(mu_t)"]
    res.loc[i, "order_L2"] = np.log(e1 / e2) / np.log(h1 / h2)

    e1, e2 = res.loc[i - 1, "Linf(mu_t)"], res.loc[i, "Linf(mu_t)"]
    res.loc[i, "order_Linf"] = np.log(e1 / e2) / np.log(h1 / h2)

# Machine-readable full summary
res.to_csv("mms_rz_viscosity_paper_table_numeric.csv", index=False)

# Paper-ready table
paper = res[["grid", "h", "L2(mu_t)", "order_L2", "Linf(mu_t)", "order_Linf"]].copy()

for col in ["h", "L2(mu_t)", "Linf(mu_t)"]:
    paper[col] = paper[col].map(lambda x: f"{x:.4e}")

for col in ["order_L2", "order_Linf"]:
    paper[col] = paper[col].map(lambda x: "--" if pd.isna(x) else f"{x:.2f}")

paper.to_csv("mms_rz_viscosity_paper_table.csv", index=False)

print()
print("Table: RZ mixing-length turbulent-viscosity verification")
print()
print(paper.to_string(index=False))
print()
print("Sanity check:")
print(f"  max |(mu_eff - mu_eff_exact) - (mu_t - mu_t_exact)| = {res['mu_eff_check'].max():.3e}")
print()
print("Wrote:")
print("  mms_rz_viscosity_paper_table.csv")
print("  mms_rz_viscosity_paper_table_numeric.csv")
print()

ok = (
    res.loc[len(res) - 1, "L2(mu_t)"] < res.loc[0, "L2(mu_t)"]
    and res.loc[len(res) - 1, "Linf(mu_t)"] < res.loc[0, "Linf(mu_t)"]
    and res["mu_eff_check"].max() < 1.0e-10
)

print("Verification result:")
if ok:
    print("  PASS")
else:
    print("  NOT PASS")
