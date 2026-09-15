# Drift mobility versus voltage

Requested on 2026-09-15 as the next of the four remaining source charts.
Implemented with existing physical and numerical parameters unchanged.
Complete and accepted by the user on 2026-09-15, including the expanded
Y range and model hole-mobility series.

## Source mapping

Inspected `data/SCLCKopecky.xlsx`, `xl/charts/chart1.xml`, and the referenced
formulas. The three series are:

- Measured: `MODEL!E9:E1003` versus `I9:I1003`, calculated here as
  `mobility.U` versus `mobility.mu_eff` from the measurement input.
- Model holes: `MODEL!Z9:Z3011` versus `AF9:AF3008`, labelled μdm (p).
  Z refers to `j(U)!C`, the hole voltage; AF refers to `j(U)!K`, whose
  formula is the hole fraction `n(E)!N` times microscopic mobility.
  Use `model_current.U_p` versus `model_mobility.mu_p` on their shared
  energy sweep. Unequal source range lengths do not justify changing grids.
- Microscopic mobility: the constant reference from `MODEL!Q2:Q3` versus
  `R2:R3`, calculated here from the existing `params['mu_0']`.

Axes are Voltage, U (V), linear, and Drift mobility, μ (m²/V/s), logarithmic.
The model legend uses μd,m (p); the measured markers use μd. Both display
limits are captured from finite measured voltage and positive finite mobility
before adding either model or reference. Negative voltage remains eligible
on the linear axis. Invalid model points leave gaps, with every energy row
retained and no clipping, sorting, resampling or recalculation.

The user subsequently requested an expanded range on 2026-09-15. The
notebook overrides the Y maximum to 10^-2 m²/V/s, showing μ0=0.0027
and the model curve across the measured voltage window. Lower Y and
X limits remain measured-based. The source model series was rechecked:
Z versus AF is the model hole series already present in this comparison.

## Physical derivation and limits

For steady single-carrier drift with immobile traps and negligible diffusion
and recombination, J = e μ0 pf F. Defining J = e μd ps F gives
μd = μ0 pf/ps; the density ratio is dimensionless, so mobility has units
m² V⁻¹ s⁻¹. This local identity needs no contact boundary condition.
See [M5](M5_model_mobility.md) for its derivation and analytical limits.

The horizontal coordinate uses the existing power-field closure with ideal
injection F(0)=0, uniform permittivity ε and thickness L:
U = e L² q / [ε(1−γ)(2−γ)]. Its derivation from Poisson's equation,
current conservation and the voltage integral is recorded in
[M4](M4_model_current.md). No new transport calculation is introduced.

The retained populations use absolute pf and excess ps, with Θ=|pf/ps|;
the voltage uses q=pt=ps−pf0. These definitions are not generally the
consistent total densities needed for a spatial transport solution.
Pairing the arrays reproduces the intended quantities but does not resolve
that existing physical limitation or force μd to agree with an inversion
of the M4 current. Values above μ0 and singularities stay intact.

## Notebook and validation

The completed comparison follows the M5 calculation and replaces the
`effective_mobility` entry used by export. The earlier measured view remains
available. The user removed the supporting energy-domain mobility figure on
2026-09-15; the M5 calculation it drew stays, feeding the model series here.

Validation executes all notebook calculation and plotting cells at their
unchanged settings, checks all three series, shared energy-row pairing,
linear/log scales and measured-only limits before the notebook override.
The expanded view also checks that μ0 lies inside the displayed Y range.
The existing 72 tests pass, including the independent Mott–Gurney mobility
inversion. This chart change adds no formulas or parameter adjustments.
