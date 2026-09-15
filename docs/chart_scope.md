# Chart reproduction scope

User scope decision, 2026-09-15: remove the total density-of-states series
from the density-of-states chart. The valence-band, conduction-band and trap
components remain, as do the energy markers and the chart's vertical range,
which is now taken from the largest plotted component.

User scope decision, 2026-09-15: remove the trap density-of-states close-up
view from the notebook. The full density-of-states chart remains. The close-up
energy grid stays in the notebook for the trap saturation check; parameters
and grids are unchanged.

User scope decision, 2026-09-15: exclude the unlabelled n(E) chart on
g(E), source S1304:S1404 (#13). It is not needed and must not be
implemented as part of this chart reproduction work. The remaining scope
contains 14 of the 15 native source charts.

User override, 2026-09-15: label the Fermi–Dirac occupation chart axes
as relative energy E−E_F (eV) and occupation probability (dimensionless).
Set X display limits to −1 to 1 eV; Y remains automatic and both axes
remain linear. Calculation parameters and grids are unchanged.

User correction, 2026-09-15: for U(nt), U(pt), j(nf), j(pf) versus
Fermi energy, the screenshot and request to fix the display supersede the
previous `10e0` interpretation. Set voltage minimum to 1 V. Set current-density
minimum to 1 A/m² to remove the tiny-current tail from the displayed range.
Upper and horizontal limits remain automatic; calculated arrays are unchanged.

User override, 2026-09-15: label the axes of that transport-versus-energy
chart explicitly: Fermi energy, E_F (eV) on X, Voltage, U (V) on left Y,
and Current density, j (A/m²) on right Y, despite absent source axis titles.

User clarification, 2026-09-14: reproduce the charts in
`data/SCLCKopecky.xlsx` with respect to labels, plotted content, and each
axis's linear or logarithmic scale. Content includes the actual series,
axis quantities, units, and energy references used by each chart.

Reproduce that content using this project's calculations and measurement
inputs. Do not copy or hard-code plotted data from the workbook, including
cached chart values. The source charts define which quantities and series
to show; they do not supply the implementation's plotted numerical data.

For comparison charts, determine automatic axis limits from the measured
data only; model curves must not widen the view. Explicit user-specified
display bounds take precedence. For charts without measured data, use
automatic limits. Do not copy source fixed bounds. Styling need not match
the source charts.

Inspect the actual chart definitions and series references before choosing
an energy offset; do not infer a chart's reference from a related formula
or column alone. Keep existing physical and numerical parameters unchanged.
Source provenance belongs in docs; code and displayed labels describe the
mathematical quantities directly.

Identify charts in discussion and task lists by their visible axis or legend
labels. Keep chart XML numbers only as secondary references. For unlabelled
charts, identify the sheet and source range explicitly.

For Trapped charge versus Free charge, determine display limits from the
measured pt and pf series only. Model curves must not widen these limits.
Align the hidden secondary Y limits with the primary Y limits and fit the
secondary X range using model rows within the measured free-charge window.

User override, 2026-09-15: for Current density versus voltage, linear/log
(measured and model hole current), set X limits to 0–4 V and Y limits to
10^-8–10^0 A/m² in the notebook. This supersedes automatic limits for this
chart only; calculation parameters remain unchanged.

User clarification, 2026-09-15: measured data determine the display limits
for comparison charts generally, not just Trapped charge versus Free charge.
For Charge density versus voltage, linear/log, explicitly set X to 0–4 V
and Y to 10^10–10^18 m⁻³. These are display limits only.

User override, 2026-09-15: expand the Drift mobility chart range to show
the microscopic mobility reference and model curve. Implemented with a
logarithmic Y maximum of 10^-2 m²/V/s in the notebook; the measured lower
Y limit and voltage limits remain. This is a display-only override.
The model hole series remains calculated on the existing energy sweep.
