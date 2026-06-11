from pathlib import Path
import sympy as sp

# Generate manufactured solution expressions for 02_mms_velocity.
#
# This script prints the exact fields
# and forcing functions that can be pasted manually into the [Functions] block
# of the MOOSE input files.
#
# Coordinate convention:
#   x = axial coordinate z
#   y = radial coordinate r

rho = 1.0
mu = 1.0e-2
lm = 0.20
FORCE_SIGN = 1.0

x, y = sp.symbols("x y")
pi = sp.pi

# Streamfunction construction:
#   uz = (1/r) dpsi/dr
#   ur = -(1/r) dpsi/dx
#
# We write the result in nonsingular polynomial form.
A = x**2 * (1 - x)**2
Ap = sp.diff(A, x)

R_over_y = y**2 - 2*y**3 + y**4
Rp_over_y = 3*y - 8*y**2 + 5*y**3

uz = A * Rp_over_y
ur = -Ap * R_over_y

p = x * y**2

uz_x = sp.diff(uz, x)
uz_y = sp.diff(uz, y)
ur_x = sp.diff(ur, x)
ur_y = sp.diff(ur, y)

# RZ strain invariant for no-swirl axisymmetric flow.
hoop = ur / y
s2 = (
    2 * uz_x**2
    + 2 * ur_y**2
    + (uz_y + ur_x)**2
    + 2 * hoop**2
)

mu_t = rho * lm**2 * sp.sqrt(s2)
mu_eff = mu + mu_t

# RZ viscous stress components.
tau_xx = 2 * mu_eff * uz_x
tau_rr = 2 * mu_eff * ur_y
tau_xr = mu_eff * (uz_y + ur_x)
tau_tt = 2 * mu_eff * hoop

# Divergence of viscous stress in RZ.
div_tau_x = sp.diff(tau_xx, x) + sp.diff(tau_xr, y) + tau_xr / y
div_tau_r = sp.diff(tau_xr, x) + sp.diff(tau_rr, y) + (tau_rr - tau_tt) / y

# Advective terms.
adv_x = rho * (uz * uz_x + ur * uz_y)
adv_r = rho * (uz * ur_x + ur * ur_y)

# MOOSE residual convention used in the MMS input:
#   adv + grad(p) - div(tau) - forcing = 0
forcing_uz = FORCE_SIGN * (adv_x + sp.diff(p, x) - div_tau_x)
forcing_ur = FORCE_SIGN * (adv_r + sp.diff(p, y) - div_tau_r)

div_check = sp.simplify(sp.diff(uz, x) + sp.diff(ur, y) + ur / y)

pin_x = 0.5
pin_y = 0.5
pin_p = float(p.subs({x: pin_x, y: pin_y}))


def moose_expr(expr):
    # Do not simplify large forcing expressions. This keeps generation fast.
    return str(expr).replace("**", "^")


expressions = {
    "exact_uz": moose_expr(uz),
    "exact_ur": moose_expr(ur),
    "exact_p": moose_expr(p),
    "forcing_uz": moose_expr(forcing_uz),
    "forcing_ur": moose_expr(forcing_ur),
}


print("RZ divergence check:")
print(div_check)
print()

print("Pressure pin:")
print(f"point = '{pin_x} {pin_y} 0'")
print(f"phi0  = {pin_p:.16e}")
print()

print("MOOSE [Functions] block:")
print()
print("[Functions]")
for name, expr in expressions.items():
    print(f"  [{name}]")
    print("    type = ParsedFunction")
    print(f"    expression = '{expr}'")
    print("  []")
    print("[]")
    print()
print("[]")

Path("manufactured_solution_functions.txt").write_text(
    "RZ divergence check:\n"
    f"{div_check}\n\n"
    "Pressure pin:\n"
    f"point = '{pin_x} {pin_y} 0'\n"
    f"phi0  = {pin_p:.16e}\n\n"
    "MOOSE [Functions] block:\n\n"
    + "\n".join(
        [
            "[Functions]",
            *sum(
                [
                    [
                        f"  [{name}]",
                        "    type = ParsedFunction",
                        f"    expression = '{expr}'",
                        "  []",
                        "",
                    ]
                    for name, expr in expressions.items()
                ],
                [],
            ),
            "[]",
            "",
        ]
    )
)

print("Wrote manufactured_solution_functions.txt")