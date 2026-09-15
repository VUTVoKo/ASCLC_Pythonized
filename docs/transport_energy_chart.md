# Voltage and current density versus Fermi energy

Status: complete and accepted by the user on 2026-09-15, including explicit
axis labels, voltage minimum 1 V and current-density minimum 1 A/m².

Implemented 2026-09-15 at the user's request. All parameters and calculated
arrays are unchanged; `plot_model_transport_energy` uses `ModelCurrent`
and is registered as `model_transport_energy` in the notebook export map.

## Source mapping

Inspected `data/SCLCKopecky.xlsx`, `xl/charts/chart10.xml` and `j(U)`
formulas. All four series use `j(U)!A9:A3009`, referencing `n(E)!A`
(absolute Fermi energy in eV). No energy offset is added.

| Legend | Y range | Calculated array | Axis |
|---|---|---|---|
| U (nt) | B9:B3009 | U_n | Left, logarithmic |
| U (pt) | C9:C3009 | U_p | Left, logarithmic |
| j (nf) | D9:D3009 | J_n | Right, logarithmic |
| j (pf) | E9:E3009 | J_p | Right, logarithmic |

The voltage series use axis IDs 432344064/432345856; current series use
1290794863/1290796783. Both X axes are linear with the same energy quantity,
so the implementation shares X. There are no source axis titles. At the user's
subsequent request, the chart explicitly labels X as Fermi energy, E_F (eV),
left Y as Voltage, U (V), and right Y as Current density, j (A/m²).
No cached values or fixed source bounds are used. With no measured series,
the shared horizontal limits and vertical upper limits are automatic.
Following the user's screenshot and request to correct the display, the
notebook sets voltage minimum to 1 V and current-density minimum to 1 A/m².
This supersedes the earlier 10 V interpretation and automatic current minimum.
These are display bounds, not changes to the arrays; automatic upper limits
remain independent, so relative curve heights depend on the two axis scales.

## Governing equations and limitations

The [M4 derivation](M4_model_current.md) starts with steady one-dimensional
Poisson and drift equations, dF/dx = e q/epsilon and J = e mu_0 f F,
constant current, ideal injection F(0)=0, and the assumed profile
F(x)=F_L(x/L)^(1-gamma). Integration and the collecting-boundary density give

    U = e L² q / [epsilon (1-gamma)(2-gamma)]
    J = e mu_0 f (2-gamma) U/L.

Here epsilon=epsilon_0 epsilon_r, q and f are number densities in m⁻³,
U is in V, and J is in A/m². The established calculation assigns
q=n_t or p_t and f=n_f or p_f. As documented in M4, these labelled trapped
densities are not generally total Poisson charge; plotting them does not
resolve that physical limitation. This change adds no calculation.
Signed values remain intact; nonpositive/nonfinite ordinates appear as
gaps on logarithmic axes, without absolute values or branch reversal.

## Validation

Executed all notebook calculation and plotting cells with their existing
parameters, skipping the export cell: passed. Checked all four series against
the calculated arrays, axis assignments, linear/log/log scales, gaps at
nonpositive values and unchanged input arrays: passed. Inspected the rendered
preview; all four branches and both vertical scales are visible.
