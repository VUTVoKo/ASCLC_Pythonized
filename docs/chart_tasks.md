# Chart inventory and reproduction tasks

Inspected 2026-09-14: chart XML, sheet/drawing relationships and selected
source-cell formulas in `data/SCLCKopecky.xlsx`; current `asclc_plotting.py`
and `notebooks/asclc.ipynb`. IDs below are the stable chart XML numbers,
not an inferred visual reading order. There are 15 native charts. The notebook
inventory below lists the figures present at that inspection, less the model
effective mobility view and the trap DOS close-up, both removed at the user's
request on 2026-09-15; figures added by later chart tasks are recorded in the
task list. Embedded images are not counted as charts.

Scope: 14 native source charts; the unlabelled n(E) chart on g(E),
S1304:S1404 (#13), was excluded by the user on 2026-09-15.
Match labels, series/quantities and linear/log scales. Calculate all
plotted values from project measurement inputs and implementations, never
from copied source chart data. Use measured-data limits for comparisons and independent styling.
Keep all existing parameters fixed. This is a task list, not implementation.

## Source chart inventory

Use visible axis or legend labels as chart names in discussion and tasks;
numeric IDs are secondary references only. Scales are X / Y. The charts have
no chart titles. Names below use Y — X axis labels, or series labels where
axis titles are absent; mathematical typography and spacing are normalized.
The unlabelled chart is identified by its sheet and source range.
Source label inconsistencies are preserved here and recorded below.

| Chart labels | ID | Sheet | Content | Scales |
|---|---|---|---|---|
| Drift mobility, μ (m²/V/s) — Voltage, U (V) | 1 | MODEL | Measured mobility, model hole mobility, microscopic mobility reference versus voltage | linear / log |
| Current density, j (A/m²) — Voltage, U (V), log/log | 2 | MODEL | Measured current density, hole and negative electron model branches, slope guides m=1 and m=2 versus voltage | log / log |
| Current density, j (A/m²) — Voltage, U (V), linear/log | 3 | MODEL | Measured current density and model hole current versus voltage | linear / log |
| charge density, nf, nt (m⁻³) — Voltage, U (V), linear/log | 4 | MODEL | Measured free/trapped holes and model free/trapped holes versus voltage | linear / log |
| Trapped charge, nt (m⁻³) — Free charge, nf (m⁻³) | 5 | MODEL | Measured trapped versus free holes, identity series, model trapped holes/electrons and free-hole identity series; additional secondary-axis model series | log / log; secondary X linear, Y log |
| g(E), dn/dEF (m⁻³eV⁻¹) — Energy, E, ΔEF (eV) | 6 | MODEL | DOS, measured/model density derivatives, model free/trapped electrons/holes and measured free/trapped holes versus energy | linear / log; secondary Y log |
| Energy, EF (eV) — Voltage, U (V); legend EF - Ec | 7 | MODEL | Measured and model Fermi energy versus voltage, equilibrium Fermi-energy reference | linear / linear |
| charge density, nf, nt (m⁻³) — Voltage, U (V), log/log | 8 | MODEL | Same four carrier series as charge density, nf, nt (m⁻³) — Voltage, U (V), linear/log (#4) | log / log |
| Theta, Θ (-) — Voltage, U (V) | 9 | MODEL | Measured and model free-charge fraction Theta versus voltage | linear / log |
| U (nt), U (pt), j (nf), j (pf) | 10 | j(U) | U(nt), U(pt), j(nf), j(pf) versus energy | linear / log; secondary Y log |
| pt (EF), nt (EF), pf (EF), nf (EF) | 11 | n(E) | Model pt(EF), nt(EF), pf(EF), nf(EF) versus Fermi energy | linear / log |
| Valenční pás, Vodivostní pás, Monoenergetická past | 12 | g(E) | Valence-band, conduction-band and localized-trap DOS versus energy | linear / log |
| Unlabelled chart on g(E), source S1304:S1404 | 13 | g(E) | Unnamed series from S1304:S1404, column headed n(E); no explicit X values | linear / log |
| fFD(E), 1-fFD(E) | 14 | FD-funkce | Fermi-Dirac occupation fFD(E) and its complement 1-fFD(E) versus energy | linear / linear |
| j (A/m²) — U (V) | 15 | Data-calculations | Measured current density versus voltage | log / log |

## Current notebook chart inventory

| Figure | Function | Content / scales | Source relationship |
|---|---|---|---|
| Measured current density | plot_measured_jv | j versus U; log/log | j (A/m²) — U (V) (#15) counterpart; labels need alignment |
| Smoothed current density | plot_smoothed_jv | Raw j and centered mean versus U; log/log | Additional diagnostic |
| Local slope | plot_local_gamma | gamma versus U, 1 and 1/2 guides; linear/linear | Additional diagnostic; fixed Y limits |
| Drift mobility | plot_effective_mobility | Measured mobility, model hole mobility and mu0 versus U; linear/log | Completed comparison (#1); expanded Y maximum by user request |
| Carrier concentrations | plot_carrier_densities | Measured pt and pf plus model pt and pf versus hole voltage; linear/log | Reproduced comparison (#4); earlier measured-only diagnostic retained |
| Carrier fraction | plot_carrier_fraction | Measured Theta versus U, 1 and 1/2 guides; linear/log | Partial Theta, Θ (-) — Voltage, U (V) (#9); model overlay absent |
| Full DOS | plot_dos | Valence-band, conduction-band and trap DOS plus energy markers; linear/log | Valenční pás, Vodivostní pás, Monoenergetická past (#12) counterpart with extra energy markers |
| Model occupations | plot_model_occupations | pt, nt, pf, nf versus EF plus EF0 marker; linear/log | pt (EF), nt (EF), pf (EF), nf (EF) (#11) counterpart with extra reference |
| Trapped charge versus free charge | plot_trapped_charge | Five primary log/log series plus independent reversed-linear/log secondary series | Reproduced source mapping; see [details](trapped_charge_chart.md) |
| Model Fermi level | plot_model_fermi_level | EF-Ev versus pf, zero reference; log/linear | Does not reproduce Energy, EF (eV) — Voltage, U (V); legend EF - Ec (#7) |
| Model current density | plot_model_current | Measured j, hole current and negated electron current versus hole voltage, m=1 and m=2 references; log/log | Implemented comparison (#2); measured-only limits |

The fixed gamma range needs review under the automatic-limit requirement.
Extra figures are listed, not authorized for deletion, except two the user
removed on 2026-09-15: the energy-domain model effective mobility view, whose
M5 calculation is retained for the drift mobility chart, and the trap DOS
close-up, whose close-up grid is retained for the trap saturation check.

## Task list

- [x] Inventory all 15 native source charts and 12 notebook figures.
- [x] Record content/scales and identify partial, missing and additional views.
- [ ] **j (A/m²) — U (V)** (#15): align axis labels (`U (V)`, `j (A/m²)`) and unnamed source series behavior; retain calculated measurement data and log/log scales.
- [x] **Drift mobility, μ (m²/V/s) — Voltage, U (V)** (#1): complete and accepted by the user on 2026-09-15. Measured mobility, model hole mobility versus hole voltage and mu0 on linear/log axes; Y maximum expanded to 10^-2 m²/V/s at user request. See [mapping and validation](drift_mobility_chart.md).
- [x] **Current density, j (A/m²) — Voltage, U (V), log/log** (#2): implemented both model series against hole voltage, the electron sign reversal, and physical m=1/m=2 reference curves with measured-only limits. Awaiting user review. See [mapping and validation](current_loglog_chart.md).
- [x] **Current density, j (A/m²) — Voltage, U (V), linear/log** (#3): complete and accepted by the user on 2026-09-15. Measured-plus-hole-current view with X limits 0–4 V and logarithmic Y limits 10^-8–10^0 A/m². See [mapping and validation](hole_current_chart.md).
- [x] **charge density, nf, nt (m⁻³) — Voltage, U (V), linear/log** (#4): complete and accepted by the user on 2026-09-15. Measured and model hole densities with user bounds 0–4 V and 10^10–10^18 m⁻³; source n/p label mismatch retained and explained. See [mapping and validation](charge_density_chart.md).
- [x] **charge density, nf, nt (m⁻³) — Voltage, U (V), log/log** (#8): complete and accepted by the user on 2026-09-15. Same four carrier series on log/log axes, with measured-only limits. See [mapping and validation](charge_density_chart.md#loglog-view).
- [x] **Trapped charge, nt (m⁻³) — Free charge, nf (m⁻³)** (#5): reproduced all six series and independent axes using calculated arrays; retained the source label mismatch. See [mapping and validation](trapped_charge_chart.md).
- [x] **g(E), dn/dEF (m⁻³eV⁻¹) — Energy, E, ΔEF (eV)** (#6): implemented measured secants, model three-point slopes and all nine series. See [mapping and validation](density_energy_chart.md).
- [x] **Energy, EF (eV) — Voltage, U (V); legend EF - Ec** (#7): implemented the absolute-shift measured formula, model hole-branch energy and EF0 reference on linear axes. See [implementation and validation](fermi_energy_chart.md).
- [x] **Theta, Θ (-) — Voltage, U (V)** (#9): added model fraction on the hole-voltage branch, with measured limits and linear/log scales. See [mapping and validation](theta_chart.md).
- [x] **U (nt), U (pt), j (nf), j (pf)** (#10): complete and accepted by the user on 2026-09-15. Labelled absolute Fermi-energy axis, voltage branches on the left logarithmic axis (minimum 1 V), and current branches on the right logarithmic axis (minimum 1 A/m²); upper and horizontal limits remain automatic. See [mapping and validation](transport_energy_chart.md).
- [x] **pt (EF), nt (EF), pf (EF), nf (EF)** (#11): aligned four occupation labels, removed additional EF0 marker and axis titles; preserved existing notebook bounds. Implemented, awaiting user review. See [mapping and validation](occupation_chart.md).
- [ ] **Valenční pás, Vodivostní pás, Monoenergetická past** (#12): align the three DOS components and labels; review the additional energy markers against source content. The total DOS series was removed at the user's request on 2026-09-15.
- **Excluded by the user, 2026-09-15: Unlabelled chart on g(E), source S1304:S1404** (#13). Not needed; removed from remaining work.
- [x] **fFD(E), 1-fFD(E)** (#14): complete and accepted by the user on 2026-09-15. Calculated occupation and complement versus labelled relative energy on linear/linear axes, with labelled probability axis and X limits −1 to 1 eV; existing parameters and grid preserved. See [mapping and validation](fermi_dirac_chart.md).
- [ ] Review extra diagnostics separately; do not count them as source-chart matches or remove them implicitly.
- [ ] Apply measured-data limits to comparison charts, honoring explicit user bounds; preserve numerical grids and other user settings.
- [ ] Validate each resulting chart's quantities, labels, series references and scales; validate any new calculations independently of source cached values.

## Source ambiguities to resolve, not silently reproduce

- Energy, EF (eV) — Voltage, U (V); legend EF - Ec (#7)'s vertical title is `Energy, EF (eV)`, but its series names say
  `EF - Ec`. `MODEL!AO9 = ('n(E)'!B9)+$B$13` restores Ec to the relative
  coordinate; `MODEL!S9` adds Ev to the measured inversion. These are
  absolute-energy formulas. The earlier suggestion to simply replace EF-Ev
  by EF-Ec is therefore insufficient for reproducing this actual chart.
- Charts 4, 5, 6 and 8 use n notation in some axis titles while several
  plotted columns contain p quantities. Trapped charge, nt (m⁻³) — Free charge, nf (m⁻³) (#5)'s secondary horizontal
  title says energy, but one series uses `MODEL!AG` (headed nfm) as X and
  `MODEL!AH` as Y. The plotted quantity and its title need reconciliation.
- g(E), dn/dEF (m⁻³eV⁻¹) — Energy, E, ΔEF (eV) (#6) has a `#REF!` series-name reference; its X/Y references still
  point to `MODEL!AQ/AR`. Several charts have unequal X/Y range lengths.
  Do not reproduce those defects by extending grids or copying caches.
- Excluded chart, retained here for provenance only: Unlabelled chart on g(E), source S1304:S1404 (#13) has no explicit X series or axis titles; its Y range is
  `'g(E)'!S1304:S1404`, with column heading n(E). It cannot be identified as
  the existing DOS close-up from this evidence.
- Charts 10–14 lack axis titles in the chart definitions. Their quantities
  above are inferred from series references and sheet headers. Rendered
  label handling must distinguish absent titles from inferred descriptions.
- Rich text uses Symbol-font glyphs and superscripts; decode these for
  mathematical labels rather than reproducing raw m/Q or plain exponents.
