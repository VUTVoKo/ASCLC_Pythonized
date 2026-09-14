# Chart inventory and reproduction tasks

Inspected 2026-09-14: chart XML, sheet/drawing relationships and selected
source-cell formulas in `data/SCLCKopecky.xlsx`; current `asclc_plotting.py`
and `notebooks/asclc.ipynb`. IDs below are the stable chart XML numbers,
not an inferred visual reading order. There are 15 native charts and now 13
current notebook figures (12 at the initial inventory). Embedded images are not counted as charts.

Scope: match labels, series/quantities and linear/log scales. Calculate all
plotted values from project measurement inputs and implementations, never
from copied source chart data. Use automatic limits and independent styling.
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
| Effective mobility | plot_effective_mobility | Measured mobility and mu0 versus U; linear/log | Partial Drift mobility, μ (m²/V/s) — Voltage, U (V) (#1); model overlay absent |
| Carrier concentrations | plot_carrier_densities | Measured pt and pf versus U, crossing annotations; linear/log | Partial charge density, nf, nt (m⁻³) — Voltage, U (V), linear/log (#4); model overlays absent |
| Carrier fraction | plot_carrier_fraction | Measured Theta versus U, 1 and 1/2 guides; linear/log | Partial Theta, Θ (-) — Voltage, U (V) (#9); model overlay absent |
| Full DOS | plot_dos | Band/trap DOS plus total and energy markers; linear/log | Valenční pás, Vodivostní pás, Monoenergetická past (#12) counterpart with extra content |
| Trap DOS close-up | plot_dos | Same quantities on close-up energy grid; linear/log | Additional view; not Unlabelled chart on g(E), source S1304:S1404 (#13) |
| Model occupations | plot_model_occupations | pt, nt, pf, nf versus EF plus EF0 marker; linear/log | pt (EF), nt (EF), pf (EF), nf (EF) (#11) counterpart with extra reference |
| Trapped charge versus free charge | plot_trapped_charge | Five primary log/log series plus independent reversed-linear/log secondary series | Reproduced source mapping; see [details](trapped_charge_chart.md) |
| Model Fermi level | plot_model_fermi_level | EF-Ev versus pf, zero reference; log/linear | Does not reproduce Energy, EF (eV) — Voltage, U (V); legend EF - Ec (#7) |
| Model current density | plot_model_current | Measured j and electron/hole currents versus respective voltages; log/log | Partial Current density, j (A/m²) — Voltage, U (V), log/log (#2); branch mapping and guides need review |
| Model effective mobility | plot_model_mobility | Electron/hole mobilities versus EF and mu0; linear/log | Supporting calculation view; does not reproduce Drift mobility, μ (m²/V/s) — Voltage, U (V) (#1) |

DOS views currently impose energy-domain X limits and a Y range relative to
the maximum DOS. These and the fixed gamma range need review under the
automatic-limit requirement. Extra figures are listed, not authorized for deletion.

## Task list

- [x] Inventory all 15 native source charts and 12 notebook figures.
- [x] Record content/scales and identify partial, missing and additional views.
- [ ] **j (A/m²) — U (V)** (#15): align axis labels (`U (V)`, `j (A/m²)`) and unnamed source series behavior; retain calculated measurement data and log/log scales.
- [ ] **Drift mobility, μ (m²/V/s) — Voltage, U (V)** (#1): combine measured mobility, calculated model hole mobility versus voltage and mu0; match mobility notation and linear/log scales.
- [ ] **Current density, j (A/m²) — Voltage, U (V), log/log** (#2): verify model voltage and current signs/branches against the referenced quantities; add both slope guides and align labels on log/log axes.
- [ ] **Current density, j (A/m²) — Voltage, U (V), linear/log** (#3): provide measured-plus-hole-current view on linear/log axes.
- [ ] **charge density, nf, nt (m⁻³) — Voltage, U (V), linear/log** (#4): add calculated model pt and pf to measured concentrations on linear/log axes; resolve source n/p label conflict.
- [ ] **charge density, nf, nt (m⁻³) — Voltage, U (V), log/log** (#8): provide the same four carrier series on log/log axes.
- [x] **Trapped charge, nt (m⁻³) — Free charge, nf (m⁻³)** (#5): reproduced all six series and independent axes using calculated arrays; retained the source label mismatch. See [mapping and validation](trapped_charge_chart.md).
- [x] **g(E), dn/dEF (m⁻³eV⁻¹) — Energy, E, ΔEF (eV)** (#6): implemented measured secants, model three-point slopes and all nine series. See [mapping and validation](density_energy_chart.md).
- [x] **Energy, EF (eV) — Voltage, U (V); legend EF - Ec** (#7): implemented the absolute-shift measured formula, model hole-branch energy and EF0 reference on linear axes. See [implementation and validation](fermi_energy_chart.md).
- [x] **Theta, Θ (-) — Voltage, U (V)** (#9): added model fraction on the hole-voltage branch, with measured limits and linear/log scales. See [mapping and validation](theta_chart.md).
- [ ] **U (nt), U (pt), j (nf), j (pf)** (#10): provide calculated voltage and current branches versus energy with separate logarithmic vertical axes; identify series-to-axis assignment.
- [ ] **pt (EF), nt (EF), pf (EF), nf (EF)** (#11): align four occupation labels and content; review additional EF0 marker against source content.
- [ ] **Valenční pás, Vodivostní pás, Monoenergetická past** (#12): align the three DOS components and labels; review additional total DOS and energy markers against source content.
- [ ] **Unlabelled chart on g(E), source S1304:S1404** (#13): establish the mathematical definition of n(E) and intended X coordinate before reproducing this incomplete source chart; do not copy its numerical series.
- [ ] **fFD(E), 1-fFD(E)** (#14): add calculated Fermi-Dirac occupation and complement on linear/linear axes using existing parameters.
- [ ] Review extra diagnostics separately; do not count them as source-chart matches or remove them implicitly.
- [ ] Apply automatic limits to reproduction charts and review existing forced limits; preserve numerical grids and other user settings.
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
- Unlabelled chart on g(E), source S1304:S1404 (#13) has no explicit X series or axis titles; its Y range is
  `'g(E)'!S1304:S1404`, with column heading n(E). It cannot be identified as
  the existing DOS close-up from this evidence.
- Charts 10–14 lack axis titles in the chart definitions. Their quantities
  above are inferred from series references and sheet headers. Rendered
  label handling must distinguish absent titles from inferred descriptions.
- Rich text uses Symbol-font glyphs and superscripts; decode these for
  mathematical labels rather than reproducing raw m/Q or plain exponents.
