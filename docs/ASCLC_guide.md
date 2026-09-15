# A-SCLC: current density, mobility, carrier densities, Fermi level shift

Model step M4 is available through `run_model_current` and `plot_model_current`.
It retains gamma = T_t/T and the established M2 density mapping; see
[M4 derivation and limitations](M4_model_current.md).

Model step M5 is available as `run_model_mobility(model, mu_0=params['mu_0'])`.
It evaluates both M2 carrier ratios at the unchanged microscopic mobility and
supplies the model hole mobility of the drift mobility chart. See
[M5 derivation and limits](M5_model_mobility.md).

## Scope

This repository reads a standard measurement CSV, plots the measured current
density against voltage, and carries out guide steps A2 to A5: the local slope
parameter `gamma`, the effective mobility of equation (4), the free and
trapped carrier concentrations of equations (5) and (6) with the parameter
`theta`, and the Fermi level shift of equation (7). It also carries the three
density-of-states distributions of the workbook's `g(E)` sheet — the two
parabolic bands and the trap — and their switched sum. M2 evaluates the
occupation sums and theta using the workbook's definitions. The spatial
transport prediction remains removed.

M2 follows the workbook's occupation sums, equilibrium subtraction, and
theta formulas. Its density definitions and known inconsistencies are recorded
in [M2: Excel conventions and known inconsistencies](M2_excel_conventions.md).
Reproduce those conventions without silently correcting them or changing
parameters; distinguish workbook agreement from physical validation.

The trap profile is **the one the workbook implements**, reproduced exactly,
including the Excel artefact that shapes it. `g(E)!M13` writes equation (S5) as

```
IF(E > E_t, M7*EXP(-(E-E_t)/kTt)/(1+EXP(-(E-E_t)/kTt)^2),
            M7*EXP( (E-E_t)/kTt)/(1+EXP( (E-E_t)/kTt)^2))     M7 = N_t/(2 k_B T_t)
```

and Excel binds `^` tighter than `+`, so `1+EXP(x)^2` evaluates as `1+e^(2x)`,
not the intended `(1+e^x)^2`. With `u = (E − E_t)/k_B T_t` both branches
therefore collapse to one symmetric profile,

```
g_t(E) = N_t/(2 k_B T_t) · e^(−|u|)/(1 + e^(−2|u|)) = N_t / (4 k_B T_t cosh u)
```

a **sech**, where equation (S5) as printed is `sech²(u/2)/4`. The `IF` is what
keeps it finite and is part of the form to port: it always evaluates the
decaying exponential, so `cosh` never overflows on the sheet's energy grid
(0 to −9.000 eV in 0.003 eV steps, 3001 points), where `|u|` reaches 745.
Evaluated this way the expression reproduces the cached `g(E)!M` column to
2.2e−16 on all 1284 rows that did not underflow to zero.

The normalization follows from the prefactor, since `∫ sech(u) du = π`:

```
∫ g_t(E) dE = (π/4) N_t = 3.69e16 m⁻³ at the shipped N_t = 4.7e16 m⁻³
```

— neither the `N_t` of equation (S5) as printed nor the `N_t/2` that this
prefactor would give with the intended denominator. The sheet's own Riemann sum
over column M is 3.6899e16, i.e. 0.78508 `N_t` against `π/4 = 0.785398`. So
`N_t` here is a scale factor rather than the trapped concentration at
saturation, and step A6 would read that saturation off as 3.69e16 m⁻³. The
parameter stays at its user-controlled 4.7e16: the `π/4` belongs to the profile,
not to the parameter.

The profile is symmetric about `E_t`, which therefore keeps the meaning
Table S2 gives it — the position of the distribution maximum — and no choice of
energy support arises. Because `k_B T_t = 2.59 meV` sits against
`k_B T = 25.8 meV` of thermal smearing, shape is nearly invisible under the
convolution of equation (S8): this sech, equation (S5) and equation (S6) agree
to within 0.4 % once normalized alike, so the `π/4` is the whole of the
difference the workbook's form makes.

## Parameters are user-controlled

Never change the sample, device, or numerical parameters supplied in the
notebook or in `examples/`. Only an explicit user instruction changing a
specific parameter overrides this. See [AGENTS.md](../AGENTS.md).

## Sources

| Source | Role |
| --- | --- |
| `data/s42005025022021.pdf` | Gavranovic, Zmeskal, Weiter & Pospisil, *Communications Physics* **8**, 280 (2025) |
| `data/42005_2025_2202_MOESM2_ESM.pdf` | Supplementary Information |
| `data/SCLCKopecky.xlsx` | The user's workbook: the smoothing, `gamma` and `E_a` columns, with cached results |
| `examples/` | Earlier exported data and calculations, kept for reference |

The workbook may be opened and cited directly; this guide quotes its cells where
they fix a parameter or define a column (see *Local-slope smoothing window*
below). Nothing in the package reads it at runtime — the modules take paths,
arrays, and parameters only.

## Frontend and backend

`notebooks/asclc.ipynb` is the **frontend**: an editable input path and
`device` dictionary, backend calls, and the plot. It contains no file parsing
or unit conversion.

The **backend** modules take paths, arrays, and parameters explicitly and
return data for the frontend to display:

- `asclc_backend.py`: `load_measurements` (standard CSV), `current_density`,
  `trailing_mean` (workbook columns E/G/N), `centered_mean` (an aligned variant
  for plotting), and `local_gamma` (column K).
  It also provides `effective_mobility` (equation (4)), `carrier_densities`
  (equations (5) and (6), with `theta`), `equal_density_voltage` (the
  `p_t = p_f` crossings), `fermi_level` (equation (7)) with `effective_dos`
  (equation (S4)), the `g(E)` sheet's three distributions — `valence_band_dos`
  (column F), `conduction_band_dos` (column I) and `trap_dos` (column M, see
  *Scope*) — with `total_dos` for their switched sum (column Q), and the
  constants `EPSILON_0`, `E_CHARGE`, `K_B`, `H`, `M_E`, `NONDEGENERATE_KT` and
  `DOS_FLOOR`.
- `asclc_workflow.py`: `run_calculation` (load + `V`, `J`),
  `run_effective_mobility` (steps A2 and A3, returning a `Mobility` with
  `gamma`, `mu_eff`, the `j` both were computed from, and `window`),
  `run_carrier_densities` (step A4 on those same rows, returning a `Carriers`
  with `p_f`, `p_t`, `p_s`, `theta`, the `mu_0` and `theta_model` used), and
  `export_results` (figures as PNG/SVG plus `measured.csv`).
- `asclc_plotting.py`: `plot_measured_jv`, `plot_smoothed_jv(result,
  window=...)` (a centred window-point mean of `J` over the raw markers, at the
  same voltages), `plot_local_gamma`, `plot_effective_mobility`,
  `plot_carrier_densities` (the axes of Fig. 4a, with any `p_t = p_f` crossing
  marked), `plot_carrier_fraction` (`theta`, the axes of Fig. S19) and
  `plot_dos(energy, ...)` (the three `g(E)` distributions against energy, the
  axes of Fig. S6), each returning `(figure, axes)` in plain publication
  styling with the notation of the published equations. `plot_dos`
  draws each band only where it exists, marks only the band edges and `E_t` that
  fall inside the energy grid it is given, and holds the view to that grid, so
  passing a narrow grid gives a close-up of the trap.

## CSV format

Standard CSV with the exact header `U(V),I(A)` and two numeric columns:
volts, amperes. A UTF-8 BOM is accepted. Incorrect headers or column
counts are rejected. Signs, row order, and nonfinite measurements are retained.

`run_calculation(path, device=device)` returns `V = U + device['voltage_offset']`
and `J = I / device['area']` (A/m²). The `device` dict is not mutated.

## Constants and units

Voltage in V; current density in A/m²; thickness in m; mobility in
m² V⁻¹ s⁻¹; carrier concentrations in m⁻³. The physical constants are the
vacuum permittivity `asclc_backend.EPSILON_0 = 8.8541878128e-12` F/m and the
elementary charge `asclc_backend.E_CHARGE = 1.602176634e-19` C (CODATA 2018),
used by equations (4) to (6). The workbook rounds both (`1.602e-19`,
`8.854e-12`), which shifts `p_f` and `p_t` by about 0.01 %.

Equations (S4) and (7) add `K_B = 1.380649e-23` J/K, `H = 6.62607015e-34` J s
(both exact) and `M_E = 9.1093837015e-31` kg. Energies are in eV on the vacuum
scale of the parameter files, so `E_c > E_v` numerically (−3.36 eV against
−5.58 eV for MAPbBr3); `k_B T` is divided by `E_CHARGE` to reach eV.
Temperatures are in K — the workbook's column N, `T(°C) + 273.15`.

## Notebook sequence

1. Imports and the editable input path and `device` dictionary.
2. `run_calculation(...)` to load the measurements.
3. `plot_measured_jv(result)` — measured J–V on logarithmic axes; nonpositive
   readings stay in the data but cannot appear on a logarithmic axis.
4. `plot_smoothed_jv(result, window=numerics["window"])` — the same data with a
   centred window-point mean of `J` overlaid at the same voltages.
5. `run_effective_mobility(result, material=material, window=numerics["window"])`
   — steps A2 and A3; then `plot_local_gamma` and `plot_effective_mobility`.
6. `run_carrier_densities(mobility, material=material, mu_0=params["mu_0"])`
   — step A4 on the same rows; then `plot_carrier_densities`,
   `plot_carrier_fraction` and `equal_density_voltage`.
7. `plot_dos(...)` — the `g(E)` sheet's three distributions over the band gap.
   The trap parameters `N_t`, `E_t` and `T_t` sit in the notebook's `params`
   dict beside `mu_0`, as `examples/MAPbBr3_S2.params` and `MODEL!B29/B25/B30`
   hold them.
8. `export_results(result, figures, directory)` to write the figures and
   `measured.csv` (columns `U_V,I_A,V_V,J_Am2`).

Changing the dataset or the `device` values is a configuration edit, not a code
change. Any measurement row count is accepted.

## Local-slope smoothing window

The workbook derives the logarithmic slope `m`, the parameter `gamma = 1/m`
(`Data-calculations` column K, headed `g (-)`), and the Arrhenius activation
energy `E_a` (column O) by a rolling straight-line fit through `ln j` versus
`ln U`. One cell sets how wide that fit is:

- **`Data-calculations!I1`** is the only control. Each fit spans rows
  `r … r + I1`, so the **window is `I1 + 1` points**. The same cell also sets
  the span of the trailing averages in columns E (`U`), G (`j`), and N (`T`).
- The window is trailing (forward-looking), not centred: `gamma` at a row
  reflects the slope just above that point, and the last `I1` rows get no
  value.
- **`I1 = 0` disables smoothing.** Every fit then sees a single point, so
  `LINEST` returns slope 0: `gamma ≡ 0`, `E_a ≡ 0`, and the "averaged" columns
  pass the raw values through. An early copy of the workbook shipped with
  `I1 = 0`; that is what made `gamma` come out zero.
- Current value: **`I1 = 10`, i.e. `window = 11`** (`numerics["window"]`).
  `trailing_mean(I, window=11) / area` reproduces the workbook's smoothed `j`
  column to three significant figures. The MAPbBr3 S2 dark scan is noisy enough
  that a 5-point window (the old `examples/MAPbBr3_S2.params` value) still looks
  like scatter; 11 points settles the sub-threshold decade into a line. Widen it
  if the curve is still jagged, narrow it if it smears the ohmic → trap-filling
  → Mott–Gurney transitions together.
- Even at `window = 11` the very lowest bias is a smoothed noise floor, not a
  real ohmic slope — the raw current there sign-flips within a few pA. Judge the
  physics only where `j` is clearly above noise.

Two backend functions port these columns, both taking `window` (the notebook
`numerics["window"]` / `examples/MAPbBr3_S2.params` `window`):

- `asclc_backend.trailing_mean(values, window=...)` — the column E/G/N average:
  point `i` is the mean of `values[i : i+window]`. The smoothed current density
  is `current_density(trailing_mean(I, window=w), area)`. Because the window is
  trailing, not centred, the smoothed series is shifted toward higher index
  (higher voltage) by about half a window.
- `asclc_backend.local_gamma(voltage, current, window=..., centered=False)` —
  column K: `gamma = d ln|U| / d ln|j|`, the `LINEST` of `ln|U|` against `ln|j|`
  over `window` points. `centered=False` is the sheet's trailing anchor;
  `centered=True` aligns the result with `voltage`. Not called by the notebook.

Both return an array the length of the input, with `nan` at the ends without a
full window and wherever a window covers a zero or nonfinite reading (the sheet
instead shrinks the window there and, for `local_gamma`, returns `0`).

### Which variable is fitted against which

The paper defines `gamma = 1/m = dlnV/dlnJ`, and the analysis takes the second
form directly: `local_gamma` fits `ln|U|` against `ln|j|` and keeps the slope.
That is column K, `LINEST(ln|U| over the window, ln|j| over the window)`.

The direction is part of the definition. Fitting `ln|j|` against `ln|U|` and
inverting is a different number on scattered data — the two agree only on a
noise-free power law — and it is not what the published workflow computes.

The sheet's `ln` columns are `LN(ABS(...))` of each **raw** row, so the slope is
fitted on the raw measurements; the trailing means in columns E and G are the
series the sheet plots, not the series it fits.

### Row layout

At each row the sheet holds three quantities taken over the same trailing
`window` starting at that row: the trailing means `U` and `j` (columns E, G),
`gamma` from the raw rows (column K), and the step-A3 mobility built from all
three. `run_effective_mobility` returns exactly that as `Mobility(U, j, gamma,
mu_eff, window)`, so `mu_eff[i]` uses `U[i]`, `j[i]` and `gamma[i]`. Plot
against `mobility.U`, the voltage each row describes — not `result.V`, which is
the raw row voltage.

For the MAPbBr3 S2 dark scan at `window = 11`, the first six rows reproduce
`Data-calculations` rows 13–18:

| sheet row | U (V) | j (A/m²) | gamma |
| --- | --- | --- | --- |
| 13 | 0.08 | 1.57e-07 | -1.45e-01 |
| 14 | 0.09 | 1.17e-07 | +1.65e-01 |
| 15 | 0.11 | 2.09e-07 | -9.66e-03 |
| 16 | 0.12 | 2.14e-07 | -2.01e-02 |
| 17 | 0.13 | 1.44e-07 | -2.79e-02 |
| 18 | 0.15 | 1.66e-07 | +5.32e-02 |

The scan's first decade of bias sits at the instrument's noise floor — ±5 pA
with sign changes inside the first window — so `gamma` there is fitted through
noise and lands near zero with either sign. That is the sheet's own result for
these rows, not a defect in the port; the physics starts where the current
clears the noise floor, around 0.3 V.

`asclc_backend.centered_mean(values, window=...)` is a non-workbook helper: the
same average centred on `i` instead of trailing, `nan` for the points at each
end where the window runs off the data. `plot_smoothed_jv` uses it so the
smoothed line stays aligned with the raw markers on the voltage axis.

## Free and trapped carriers, and theta (step A4)

Step A4 takes the rows of step A3 — the same trailing window, so `p_f[i]`,
`p_t[i]` and `theta[i]` describe `mobility.U[i]`, `mobility.j[i]` and
`mobility.gamma[i]` — and the hand-selected microscopic mobility `mu_0` read off
the step-A3 plot.

### The two equations

Equation (5) is Ohm's law in the form the SI uses for these devices
(equation (S14), `j = e mu_0 p_f (2 - gamma) U / L`) solved for `p_f`:

```
p_f = L j / (e mu_0 (2 - gamma) U)
```

Equation (6) is what remains when equation (2) is divided by that:

```
eps_0 eps_r Theta (1 - gamma)(2 - gamma) U / L^2 = e p_f
    =>  p_f / Theta = eps_0 eps_r (1 - gamma)(2 - gamma) U / (e L^2)
```

Both sides carry m⁻³: `[F/m][V]/([C][m²]) = m⁻³` and
`[m][A/m²]/([C][m²V⁻¹s⁻¹][V]) = m⁻³`. `mu_0` enters `p_f` only, so it slides the
free-carrier curve up and down and leaves equation (6) alone.

### Which concentration equation (6) computes

The derivation above returns `p_f / Theta`, that is the **total** concentration
`p_s`, with the trapped part `p_t = p_s - p_f`. The article labels equation (6)
`p_t` and then defines `Theta = p_f/(p_f + p_t)` (A4 and equation (S12)); those
two statements cannot both hold with equations (2) and (4), because reading
equation (6) as `p_t` makes `mu_eff/mu_0 = p_f/p_t` identically. Exact
Mott–Gurney data (`gamma = 1/2`, `mu_eff = mu_0`) would then give `Theta = 1/2`,
where the article says `Theta = 1` for that region.

The workbook carries both readings at once: `MODEL!O = K + M` adds equation (5)
and equation (6) as if equation (6) were `p_t`, while the `Theta` that feeds its
plots, `MODEL!Q = I/B21`, is `mu_eff/mu_0`.

`theta_model` selects the reading; all three give the same `p_f`:

| `theta_model` | `p_t` | `p_s` | `Theta` |
| --- | --- | --- | --- |
| `free_over_total` (default) | Eq. (6) − `p_f` | Eq. (6) | `p_f/p_s` = `mu_eff/mu_0` |
| `absolute_over_total` | Eq. (6) | `abs(p_f) + abs(p_t)` | `abs(p_f)/p_s` (SI, Eq. (S12)) |
| `mobility_ratio` | Eq. (6) | `p_f + p_t` | `mu_eff/mu_0` = `p_f/p_t` (workbook `MODEL!Q`) |

The default reading makes equations (2), (4), (5), (6) and (S12) one system with
a single `Theta`: 1 in the ohmic and Mott–Gurney regions, 1/2 where
`p_t = p_f` — and there, by equation (4), `mu_eff = mu_0/2`, which is the test of
`mu_0` the article proposes. `absolute_over_total` is the SI read literally, its
absolute values keeping the sign-flipping noise floor finite;
`mobility_ratio` reproduces the workbook column. The three differ only in
`p_t` and `Theta`, and only materially once `Theta` is not small: on the
MAPbBr3 S2 dark scan the largest `p_t` moves by 18 % between the default and the
other two.

### Conventions and what is left unclipped

Equations (5) and (6) are written for the forward branch, so `carrier_densities`
takes `abs(j)` and `abs(U)`, the convention of equation (4) in
`effective_mobility`; the reverse branch of a symmetric device maps onto the same
positive concentrations. A row is analysed or it is not: all four outputs are
`nan` where `U`, `j` or `gamma` is nonfinite (including the ends without a full
window), where `U = 0`, and where `gamma = 2` (equation (5) singular). `theta` is
additionally `nan` where its denominator vanishes — `gamma = 1`, the ohmic point,
where equation (6) leaves no space charge at all.

Two signals are deliberately left in the output rather than clipped: `1 < gamma <
2` makes equation (6) negative, the same report equation (4) makes with a
negative mobility, and `Theta > 1` (so a negative `p_t` under the default
reading) means the row needs more mobility than the `mu_0` supplied. Since `mu_0`
is hand-selected, that is how it gets judged.

### The MAPbBr3 S2 dark scan

At `window = 11` and the shipped `mu_0 = 2.7e-3 m² V⁻¹ s⁻¹`, 204 of the 214 rows
are analysed. `p_t` runs from 7.4e14 m⁻³ at the noise floor to 1.95e16 m⁻³
(1.95e10 cm⁻³) near 3 V, `p_f` from about 1e12 m⁻³ in the sub-threshold decade to
4.3e15 m⁻³, rising three decades through the trap-filling step just below 2 V.
`Theta` peaks at 0.222 near 2.5 V, so `p_t` stays between 3.5 and 3900 times
`p_f` and `equal_density_voltage` returns nothing: this scan never reaches the
`p_t = p_f` point, and the Mott–Gurney law never applies to it at this `mu_0`.
That is the same statement as `mu_eff` reaching only 0.22 × `mu_0` in step A3.

## Fermi level shift (step A5)

Step A5 turns the free-carrier concentration of step A4 into the position of the
quasi-Fermi level. Equation (7),

```
p_f = N_v exp(-dE_F / (k_B T))
```

defines the Fermi level shift `dE_F = E_F - E_v` against the valence band edge,
and `fermi_level` inverts it:

```
dE_F = k_B T ln(N_v / p_f)
E_F  = E_v + dE_F
```

`fermi_level` is a backend function; the notebook sequence above stops at step
A4, and there is as yet no `run_fermi_level` workflow step or `E_F`-versus-`V`
plot.

`N_v`, the concentration of delocalized valence states, is `effective_dos`, the
mono-energetic level of equation (S4):

```
N_v = 4 pi m_h* k_B T / h^3 sqrt(2 pi m_h* k_B T) = 2 (2 pi m_h* k_B T / h^2)^1.5
```

The two forms are the same number (`4 pi sqrt(2 pi) = 2 (2 pi)^1.5`); the second
is the textbook form and the workbook's `g(E)!F7`. It carries m⁻³:
`[kg][J/K][K]/[J s]^2 = kg/(J s^2) = m^-2`, and `(m^-2)^1.5 = m^-3`.

### `E_F0` stays hand-selected

Equation (7) fixes `E_F` against `E_v` alone. The thermodynamic Fermi level
`E_F0` — `E_F` at 0 V, one of the five Table S2 `Variable` parameters — does not
enter the calculation and is not fitted by it. The two are independent
statements, which is what makes the comparison in Fig. 5 informative: the
low-voltage end of the computed curve is the test of the hand-selected `E_F0`.
On the MAPbBr3 S2 dark scan `E_F` runs from −5.08 eV at the noise floor to
−4.74 eV at the top of the sweep, straddling the shipped `E_F0 = −4.84 eV`.

### Direction, degeneracy, and the workbook's `ABS`

On the vacuum scale `p_f < N_v` puts `E_F` above `E_v`, inside the gap, and
injecting free holes moves it back down toward the valence band — the direction
the article reports for MAPbBr3, whose Fermi level approaches the transport band
as the voltage rises.

`nondegenerate` flags `dE_F >= 3 k_B T` (equivalently `p_f <= N_v e^-3`, about
`0.0498 N_v`), the condition of Supplementary Note 1 under which equation (S4)
may replace the parabolic band and Boltzmann statistics replace Fermi–Dirac.
The SI states it as a *maximum* distance of `3 k_B T` from the band edge; that
is the wrong way round for the approximation it introduces, which needs the
Fermi level at least `3 k_B T` inside the gap, and `3 k_B T` is a convention
either way. Rows that fail it are returned, not dropped: there equation (7)
understates how far the level has moved. The MAPbBr3 S2 dark scan never comes
close — `p_f / N_v` peaks at 1.4e-9, so every analysed row is non-degenerate.

The workbook's `MODEL!S` column is this equation with an `ABS` around the shift,
`ABS(k_B T/e ln|p_f/N_v|) + E_v`, which would fold a degenerate row back into
the gap. `fermi_level` leaves the sign alone and reports the case through
`nondegenerate` instead.

A row is analysed or it is not: `E_F` and `N_v` are `nan` and `nondegenerate` is
`False` wherever `p_f` is nonpositive or nonfinite — the noise floor, the sign
flips, the ends without a full window — or the temperature is nonpositive or
nonfinite. `temperature` is either one value for the scan or one per point.

### Temperature in two places

`fermi_level` evaluates `N_v` and the Boltzmann exponent at the same
temperature. The workbook does not: `g(E)!F7` fixes `N_v` at the single sheet
temperature `Data-calculations!B6`, while `MODEL!S` feeds the per-row column N
into the exponent. Over this scan's ±16 K spread that choice, together with the
sheet's rounded constants, is the whole difference between the two: replaying
the sheet's own arithmetic reproduces its cached `MODEL!S` column to 9e-16 eV,
and `fermi_level` on the same `p_f` differs from it by at most 2.4 meV.

## Validation

`test_asclc_csv.py` checks standard CSV loading: unit/column enforcement and the
preservation of signs, order, and nonfinite rows. `test_asclc_smoothing.py`
checks `trailing_mean` and `centered_mean` against explicit windows, alignment,
`window = 1` identity, and `nan` handling. `test_asclc_gamma.py` checks `local_gamma` against analytic
power laws (`j proportional to V**k` gives `gamma = 1/k`), current scale/sign
invariance, the centred variant, and the `nan` handling. Because both fit
directions agree on a noise-free power law, it also pins the direction on
scattered data — `gamma` must equal `1/(slope of ln j on ln U)` and must differ
from the reversed regression — and checks the ohmic and Mott–Gurney reference
slopes and the `m = 0` case. `test_asclc_mobility.py`
checks `effective_mobility` against the independently written Mott–Gurney law
(exact round trip at `gamma = 1/2`, for several mobilities), the shape-factor
scaling away from `gamma = 1/2`, the `L^3` and `1/eps_r` dependences, the unit
scale, sign invariance and forward/reverse symmetry, the negative result for
`gamma > 1`, and the `nan` and validation cases. `test_asclc_carriers.py` checks
step A4 against independently written forms of equations (3) and (S14):
equation (5) returns the `p_f` that equation (S14) was built from at six slopes,
equation (6) matches a hand-checked magnitude, exact Mott-Gurney data give
`p_t = 0` and `Theta = 1`, `Theta` equals `mu_eff/mu_0` on scattered data,
`p_t = p_f` coincides with `mu_eff = mu_0/2`, and the three `theta_model`
readings are pinned against each other. It also checks the scaling in `mu_0`,
`L` and `eps_r`, forward/reverse symmetry, the sign-flipping noise floor, the
singular slopes, and the `nan` and validation cases. `equal_density_voltage` is
checked against a crossing placed at a known voltage. `test_asclc_fermi.py`
checks step A5: equation (7) is written out forward and `fermi_level` must
return the `E_F` it was built from, for a scalar and a per-point temperature;
`effective_dos` is checked against its own definition and its `T^1.5` scaling;
per-point temperatures must agree with the one-point calls they are made of and
must move `E_F` where a scalar temperature cannot; and the degeneracy flag, the
`nan` rows and the validation cases are pinned. Independently of the tests, the
equation is checked against the workbook as described above.
`test_asclc_trap.py` checks the trap profile: `trap_dos` must reproduce seven
cached `g(E)!M` values spanning 63 decades of the profile (passing the `T_t`
that matches the sheet's rounded `k_B T_t`, so the form is tested rather than
the constants), and must match an independently written `N_t/(4 k_B T_t cosh u)`
away from the tails. Its integral is checked against `(pi/4) N_t` and against
the 3.69e16 m⁻³ saturation, its peak against `N_t/(4 k_B T_t)`, and its symmetry
about `E_t` exactly. The far tail is checked to underflow to zero, monotonically
and without overflow, at the `|u| = 745` the sheet's own grid reaches — where
`cosh u` alone returns `inf`. The scalings in `N_t`, `T_t` and `E_t`, the
nonfinite energies and the validation cases are pinned. `test_asclc_dos.py`
checks the two bands and the sum: `valence_band_dos` and `conduction_band_dos`
must reproduce the cached `g(E)!F` and `g(E)!I` columns, and must do so by a
single constant factor across every band row — the sheet's rounded `h`, `e` and
`m_0` scale its `C` by 2.3e-4 and change nothing else — with the `DOS_FLOOR`
rows matched exactly. Both are checked against an independently written
`(1/2 pi^2)(2m/hbar^2)^1.5 sqrt(E)` in joules, against each other under mirrored
energies and a mass ratio, and, decisively, against `effective_dos`: the
Boltzmann integral of each band must return the `N_v` (`N_c`) of equation (S4),
which ties the energy-resolved profile to the integrated one. `total_dos` is
checked to be the exact sum of the three, to drop each term entirely when its
switch is off — floor included, as the sheet's `0*F13` does — and to reproduce
the cached `g(E)!Q` column where a band dominates it. Nonfinite energies,
shapes and the validation cases are pinned.
