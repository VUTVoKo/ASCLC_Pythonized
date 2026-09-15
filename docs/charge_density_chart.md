# Charge density versus voltage, linear/log

Implemented on 2026-09-15 as the next chart after the accepted linear/log
hole-current chart. The notebook calls `plot_carrier_densities` with model
occupations and current outputs and registers the comparison under
`carrier_densities` for export. The earlier measured-only diagnostic remains.

Complete and accepted by the user on 2026-09-15 with X limits 0–4 V and
Y limits 10^10–10^18 m⁻³, and the measured-data display-limit rule recorded
prominently in `AGENTS.md`.

## Source mapping

Inspected `xl/charts/chart4.xml` and referenced formulas in
`data/SCLCKopecky.xlsx`:

| Series | Source X / Y | Calculated X / Y |
| --- | --- | --- |
| pt (m⁻³) | MODEL E / M | carriers.U / carriers.p_t |
| pf (m⁻³) | MODEL E / K | carriers.U / carriers.p_f |
| ptm | MODEL Z / AJ | model_current.U_p / model.p_t |
| pfm | MODEL Z / AH | model_current.U_p / model.p_f |

Z9 references j(U)!C9; AJ9 and AH9 reference n(E)!O9 and L9.
The axis says `charge density, nf, nt (m⁻³)` while the series contain holes.
The requested source labels are retained; this does not change the densities
into electron populations. The unequal source pfm range lengths are not
propagated. Model arrays must share the same Fermi-energy rows.
No source cached values are used for plotting.

## Physical definitions and limits

This chart adds no calculation or parameter changes. Under the steady
one-dimensional drift and power-field closure derived in
[M4](M4_model_current.md), with ideal injecting boundary F(0)=0,
F_L=(2-gamma)U/L. Drift J=e mu_0 p_f F_L gives
p_f=L J/[e mu_0 (2-gamma)U]. Poisson gives the total boundary density
p_s=epsilon(1-gamma)(2-gamma)U/(e L²). The existing measurement convention
uses voltage/current magnitudes and p_t=p_s-p_f. Each density has units m⁻³;
the measured decomposition conserves p_s=p_f+p_t.

The source M column instead labels the Poisson expression itself pt.
The comparison retains the project's measured decomposition, as documented
in the backend, rather than replacing it with that source expression.
Model p_f is the absolute valence-band hole integral, while
p_t=A_p(E_F)-A_p(E_F0)-p_f0, with DOS occupation definitions in
[M2](M2_excel_conventions.md). Its voltage uses this labelled p_t as q in
U_p=e L² q/[epsilon(1-gamma)(2-gamma)]. As already documented in M2/M4,
that model q is not generally a consistent total space-charge density.
The overlay preserves this unresolved convention; it does not establish
physical consistency between measured and model trapped populations.

Voltage is linear and retains finite signed values; density is logarithmic.
Nonpositive or nonfinite densities cannot be displayed. Invalid model rows
leave gaps, with no absolute values, clipping, interpolation, or resampling.
The plotting function determines default limits from measured data before
adding model curves. User correction, 2026-09-15: the notebook explicitly
sets X to 0–4 V and Y to 10^10–10^18 m⁻³. Model data must not expand
comparison-chart limits. This corrects the initial full-model automatic view.

## Validation

Executed every notebook calculation and plotting cell with the existing
settings, excluding export. Verified all four plotted series against their
calculated arrays, their row pairing, linear/log scales, measured-only default limits and the explicit notebook bounds.
Inspected the rendered comparison and saved its output in the new notebook
cell. Existing physical and numerical parameters remain unchanged.

## Log/log view

Complete and accepted by the user on 2026-09-15.

Implemented on 2026-09-15 in the next notebook section using
`plot_carrier_densities(..., xscale="log")`, registered for export as
`carrier_densities_loglog`. Inspected `xl/charts/chart8.xml`: its four
series use the same MODEL E/M, E/K, Z/AJ and Z/AH mapping above, with
base-10 logarithmic axes and the same axis and legend labels. Its unequal
model source ranges are not propagated, and cached values are not used.

The physical definitions, derivation and unresolved model density convention
above also apply here. This view adds no calculations or parameter changes.
Both coordinates must be finite and positive to appear on logarithmic axes;
invalid model rows leave gaps. Limits are captured from the two measured
series before adding model curves. The linear/log chart's explicit bounds
apply only to that chart. Automatic log/log bounds with the current inputs
are approximately 0.06449–3.51576 V and 3.098e11–3.304e16 m⁻³.

Validation: executed all notebook calculation and plotting cells, excluding
export, with unchanged settings. Verified all four lines against their
calculated arrays and paired rows, both logarithmic scales, and both limits
against an independent measured-only plot. Confirmed the existing linear/log
view retains its explicit bounds. Inspected the rendered chart and saved its
output in the new notebook cell.
