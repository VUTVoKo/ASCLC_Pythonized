# A-SCLC: a guide to a reusable analysis and modeling notebook

The project is a reusable A-SCLC tool with a frontend/backend structure.
`notebooks/asclc.ipynb` is the frontend for user-supplied measurements and
material/device parameters; reusable functions are defined in separate Python
modules and called from the notebook. Supplied measurements and parameters are
examples, not fixed application inputs. Calculations use full-precision physical
constants and explicitly selected physical and numerical conventions.

The XLSX is used only for comparison tests. Application code and notebook cells
must neither mention nor depend on it. Exact agreement is not a requirement;
small differences are expected, and the user decides what difference is too much.
Never change constants, formulas, grids, or parameters merely to match it.

Build and inspect one stage at a time together. This document is a reference for
discussion, not an instruction to implement the entire model at once.

## 0. Sources, scope, and evidence

| Source | Role |
| --- | --- |
| `data/s42005025022021.pdf` | Gavranovic, Zmeskal, Weiter & Pospisil, *Communications Physics* **8**, 280 (2025); main Eqs. (1)-(7) |
| `data/42005_2025_2202_MOESM2_ESM.pdf` | Supplementary Information; DOS and occupation equations on pp. 3-4, parameter roles on p. 5, A1-A7 / M1-M5 on pp. 6-7 |
| `examples/` | Earlier exported data and calculations; check conventions before comparison |

Distinguish published relations from unresolved physical interpretations.
The supplied PDFs define the scientific sources; comparison-test data does not
define application behavior. See section 7 for the test-only comparison policy.

Two calculation paths share material and device parameters:

1. **Analysis (A):** measured voltage and current -> local slope, effective
   mobility, carrier-density estimates, and Fermi level.
2. **Modeling (M):** DOS and an energy sweep -> occupied densities -> voltage and
   current. Parameters are adjusted by hand and results compared with analysis.

Reusability means the user can change the dataset, material/device parameters,
and declared calculation conventions without rewriting backend functions.
It does not imply that one transport model applies to every material or regime.
The hole-injection equations below describe the first calculation path to develop;
other carrier types and extraction require their corresponding equations, not
just a different material name. Automated fitting is not required for the current
hand-fitting workflow.

### Frontend and backend responsibilities

The notebook is the **frontend**: editable configuration, backend function calls,
result inspection, tables, and plots. It should not contain implementations of
file parsing, validation, unit conversion, or transport calculations.

Separate Python modules form the **backend**. Functions receive paths, arrays,
parameters, and conventions explicitly, and return data/results for the frontend
to display. They must not depend on notebook globals, example filenames, or
interactive display calls. Backend validation raises informative exceptions or
returns flags; the notebook decides how to present them.

### Current calculator interface

The notebook accepts standard CSV files with the exact header `U(V),I(A),T(C)`
and three numeric columns: volts, amperes, and Celsius. Delimiter, header count,
column mapping, and unit selection are no longer notebook settings. The backend
rejects incorrect headers or column counts and retains signs, order, and
nonfinite measurements. A UTF-8 BOM is accepted.

- `asclc_backend.py`: standard CSV loading, current density, adjacent log slopes.
- `asclc_model.py`: DOS, occupations, injected densities, transport and slope guides.
- `asclc_workflow.py`: grid construction, complete calculation, and file exports.
- `asclc_plotting.py`: measured, slope, mobility, density, and combined model
  plots in plain publication styling; each function returns `(figure, axes)`.
  Axis titles and series names follow the reference workbook's chart labelling
  (its module docstring maps each function to a workbook `Graf`) so the figures
  can be compared against those charts; the workbook's colours and fonts are not
  reproduced.

The notebook keeps the input path and `material`, `device`, `model`, and
`numerics` dictionaries visible, then calls `run_calculation(...)`. The returned
result exposes raw `measurements`, corrected `V`, current density `J`, measured
`T_K`, interval midpoints `V_mid`, slopes `m`, model arrays, and plot series.
Model calculations use `device['temperature']` in kelvin, independently of the
measured temperature array. Both measured plots and slopes use the configured
voltage offset. No spreadsheet or precomputed model is loaded.

Energy and Fermi ranges and their common step are editable. The energy stop is
exclusive (rectangular integration); the Fermi stop is inclusive. Exporting is
an explicit notebook call to `export_results(result, figures, directory)`, where
`figures` maps a base filename to a Matplotlib figure; it writes each figure as
PNG/PDF/SVG alongside the measured and calculated CSVs.

The model conventions below remain documented alternatives and limitations;
the example selects `logistic` (SI S5), `voltage_density="reference"`, and
`gamma=T_t/T`. Parameters and guide multipliers remain explicit settings.

## 1. Constants, units, and example parameters

Physical constants are imported from `scipy.constants` in `asclc_model.py` and
shared by the model and analysis backend. Do not duplicate their numeric values
in application code. Material and device parameters are supplied per run.
The installed SciPy 1.18.1 supplies all published digits from the
[NIST CODATA 2022 listing](https://physics.nist.gov/cuu/Constants/Table/allascii.txt):

```python
from scipy.constants import (
    elementary_charge as E_CHARGE,  # C (exact SI definition)
    Boltzmann as K_B,               # J/K (exact SI definition)
    epsilon_0 as EPS0,              # F/m (measured)
    Planck as H,                    # J s (exact SI definition)
    electron_mass as M_E,           # kg (measured)
)
T0 = 273.15               # Celsius-to-kelvin offset
```

These are stored as ordinary double-precision floating-point values. Exact SI
values still have floating-point representation limits; EPS0 and M_E additionally
have measurement uncertainty. Do not truncate constants or offer a rounded
compatibility set. The equations below use `e`, `k_B`, `eps0`, `h`, and `m_e` as
notation for the corresponding constants above.

Energies use eV on the vacuum scale; densities use m^-3; DOS uses m^-3 eV^-1;
lengths use m; mobility uses m^2 V^-1 s^-1; current density uses A m^-2.
Define `kT_eV = k_B*T/e` and `kTt_eV = k_B*T_t/e` explicitly. In exponential
expressions, `exp(...)` denotes the exponential function; `e` denotes charge.

The example run uses the following configurable parameters; these are not
universal material defaults.

| Parameter | Example value |
| --- | --- |
| E_c, E_v | -3.36, -5.58 eV |
| Relative hole/electron masses | 0.305, 0.32 |
| Relative permittivity | 25.5 |
| Thickness L | 6.0e-4 m |
| Area S | 7.7e-6 m^2 |
| Temperature T | 299 K |
| Mobility mu0 | 2.7e-3 m^2 V^-1 s^-1 |
| Equilibrium Fermi level E_F0 | -4.84 eV |
| Trap position E_t | -4.82 eV |
| Trap concentration N_t | 4.7e16 m^-3 |
| Trap temperature T_t | 30 K |

Record parameters and modeling conventions with every run. Device dimensions
must describe the configured sample; Supplementary Table S1 lists 0.80 mm for
MAPbBr3, while this example uses 0.60 mm.

## 2. Equations and density definitions

### Transport equations

Write `eps = eps0*eps_r` and use `gamma` for gamma:

| Source | Relation |
| --- | --- |
| Main Eq. (1) | `J = e*mu0*p_f*V/L` |
| Main Eq. (2) | `J = eps*mu0*Theta*(1-gamma)*(2-gamma)**2*V**2/L**3` |
| Main Eq. (3) | `J = (9/8)*eps*mu0*Theta*V**2/L**3` |
| Main Eq. (4) | `mu_eff = L**3*J/(eps*(1-gamma)*(2-gamma)**2*V**2)` |
| Main Eq. (5) | `p_f = L*J/(e*mu0*(2-gamma)*V)` |
| Main Eq. (6), labeled p_t in the article | `q_eq6 = eps*(1-gamma)*(2-gamma)*V/(e*L**2)` |
| Main Eq. (7) | `p_f = N_v*exp(-(E_F-E_v)/kT_eV)` |
| SI Eq. (S12) | `Theta = p_f/(p_f+p_t)` |
| SI Eq. (S14), hole drift term used here | `J = e*mu0*p_f*(2-gamma)*V/L` |

`q_eq6` is a temporary neutral name for a **number density**, not charge in
coulombs. It avoids deciding prematurely whether Eq. (6) represents trapped,
total, or injected density. The factor `(2-gamma)` distinguishes Eq. (5) from
a direct inversion of Eq. (1).

From Eqs. (4) and (5), `mu_eff/mu0 = p_f/q_eq6`. Consequently, treating Eq. (6)
as total density makes these transport relations algebraically consistent with
S12. Treating it as trapped density and adding `p_f` gives a different concentration ratio. This is an unresolved physical
interpretation; retain explicit column names instead of silently switching them.

### Absolute densities versus changes from equilibrium

For an occupation integral, distinguish:

```text
P_abs(E_F)       = integral of g(E)*hole_occupation(E, E_F)
Delta_P(E_F)     = P_abs(E_F) - P_abs(E_F0)
p_f_abs(E_F)     = valence-band occupation integral
Delta_p_f(E_F)   = p_f_abs(E_F) - p_f_abs(E_F0)
```

For matching integration domains, an absolute trapped density is
`P_abs - p_f_abs`; its change is `Delta_P - Delta_p_f`. Neither should be
confused with `Delta_P - p_f_abs`: at equilibrium the latter is negative.
The full model DOS integral also includes the conduction band, so its absolute
hole count should not automatically be interpreted as a physical hole total.

SI S8/S11 restrict integration to parts of the gap using E_F0 as a boundary.
That is **not mathematically equivalent** to subtracting equilibrium occupation
over the full DOS. A restricted integral can be nonzero at equilibrium, whereas
an equilibrium-subtracted integral is zero there by construction. Keep these
as distinct model choices until their physical interpretation is resolved.

### Density of states and occupations

For an eV energy grid, the band DOS can be written as:

```python
C_v = 4*pi*(2*m_e*m_eff_h*e)**1.5/h**3
C_c = 4*pi*(2*m_e*m_eff_e*e)**1.5/h**3
g_v = C_v*sqrt(max(E_v-E, 0))
g_c = C_c*sqrt(max(E-E_c, 0))
N_v = 2*(2*pi*m_e*m_eff_h*k_B*T/h**2)**1.5
N_c = 2*(2*pi*m_e*m_eff_e*k_B*T/h**2)**1.5
```

Here `max` and `sqrt` are elementwise for arrays. Band DOS is zero outside
the corresponding band.

The calculator supports only the published trap profiles: `logistic` for the
biexponential distribution in SI S5, and `gaussian` for SI S6:

```text
S5: u = (E-E_t)/kTt_eV
    g_t = (N_t/kTt_eV)*exp(u)/(1+exp(u))**2
S6: sigma_eV = 2*k_B*T_t/e
    g_t = N_t/(sigma_eV*sqrt(2*pi))*exp(-(E-E_t)**2/(2*sigma_eV**2))
```

Both integrate to N_t over an infinite energy interval. Finite-grid integrals
have truncation and quadrature errors. For numerical evaluation of S5, use
`z = exp(-abs(u))` and `g_t = (N_t/kTt_eV)*z/(1+z)**2` to avoid overflow.

The notebook uses S5, so N_t is the integrated trap concentration.

Evaluate occupations directly, for example with `scipy.special.expit`:

```python
f_e = expit((E_F-E)/kT_eV)
f_h = expit((E-E_F)/kT_eV)
```

Hole occupation increases **above** the Fermi level. Computing `f_h` directly
avoids losing tiny hole populations through subtraction of `1-f_e`.

## 3. Measured-data analysis (A1-A7)

Apply this path to the configured measurement file and parameters.

1. Read U, I, and T from the standard CSV. Calculate `V = U + V_min`,
   `J = I/S`, and `T_K = T_C + 273.15`.
2. Calculate adjacent log differences for m and retain invalid pairs as NaN.
   Measured gamma is `1/m`; windowed slope fitting is not implemented.
3. Calculate `mu_eff`, `p_f`, and `q_eq6` from the transport equations.
4. State the density and Theta conventions used for the run; their unresolved
   interpretations are described in section 2.
5. When implementing Fermi-level inversion, use `E_v-kT_eV*log(p_f/N_v)`
   for positive p_f and check the nondegenerate-carrier regime of Eq. (7).
6. Plot densities against voltage and, when available, Fermi energy; inspect
   them before using plateaus to estimate parameters.
7. Interpret carrier type using the band/contact context.

For temperature-dependent analysis, explicitly choose measured per-point or
configured scalar temperature, and evaluate N_v at the same temperature as kT.
Model calculations currently use the configured scalar device temperature.

Keep original rows. Separate numerical finiteness from physical screening.
A useful SCLC screen requires positive V and J, finite outputs,
`0 < gamma < 1`, and `0 < Theta <= 1` with its definition stated. These are
necessary conditions, not proof of validity. Eq. (4) is singular at gamma = 1.

### Implemented A3 effective mobility

`asclc_backend.effective_mobility` evaluates main Eq. (4) with measured
`gamma = 1/m`, independently of the prescribed model gamma. Both V and J
are evaluated at geometric interval midpoints, consistent with the adjacent
logarithmic slopes. Returned arrays retain all intervals in acquisition order.
Invalid endpoint pairs and gamma outside `(0, 1)` produce NaN mobility;
gamma within an absolute `1e-12` of one is treated as the singular ohmic limit.
No smoothing is applied. Near-ohmic noise can strongly amplify mobility.

`asclc_plotting.plot_effective_mobility` overlays the measured effective
mobility, the modeled effective mobility `mu0*abs(p_f/delta_p)` over the measured
voltage interval, and a `mu0` reference line. Voltage is linear and mobility
logarithmic. The reference line defaults to the model `mobility`; an optional
`mu0_estimate` in m^2/(V s) overrides its value without modifying the model.
Estimate microscopic mobility only from a supported transport regime
(for example, a high-voltage Mott-Gurney plateau); a U-shaped curve's minimum
alone does not identify it. Ohmic estimation requires additional carrier-density
information rather than substituting gamma = 1 into Eq. (4).

`result.J_mid`, `result.gamma`, and `result.mu_eff` accompany `result.V_mid`
and `result.m`. Export writes `measured_mobility.csv` with these five columns and
saves every figure collected in the notebook in PNG/PDF/SVG formats.

### Implemented A4 carrier concentrations

`asclc_backend.carrier_densities` evaluates Eq. (5) for `p_f` and Eq. (6)
for `p_t`, following the trapped-density convention stated in SI A4. It returns
`p_s = p_f + p_t` and `theta = p_f/p_s`. `mu_eff/mu0` instead equals `p_f/p_t` under this convention. The unresolved
physical interpretation of Eq. (6) remains as described in section 2.

The workflow uses A3's geometric V/J midpoints and measured gamma, with the
same nonsingular `(0, 1)` gamma screen. Invalid intervals remain NaN in acquisition
order. Densities are in m^-3. `run_calculation(..., analysis_mobility=...)` accepts
a positive finite microscopic mobility in m^2/(V s); omission uses the model's
mobility. The notebook exposes this input explicitly. Updating the A3 reference
line alone does not update A4; change `analysis_mobility` and rerun the calculation.

Results expose `p_f`, `p_t`, `p_s`, `theta`, and `analysis_mobility`, separately
from model densities in `result.model`. `asclc_plotting` draws `p_f`, `p_t`,
`p_s` and Theta against voltage on logarithmic voltage axes. At `p_f = p_t`, this
Theta is 0.5; equality alone is not an automatic identification of a transport
regime. Export adds `measured_carriers.csv`, including the mobility input and
explicitly labeled Eq. (6)/Theta columns, plus the density and fraction figures.

## 4. Model calculation (M1-M5)

### Configurable grids and voltage-density conventions

Energy limits, step size, Fermi sweep, and DOS contributions are run settings.
Choose them for the supplied band edges and temperatures and check convergence.
The calculator evaluates occupations directly for each Fermi level and integrates
with rectangular weights of `energy_step`. The workflow exposes a common step,
with an exclusive energy stop and inclusive Fermi stop.

For the hole branch:

```text
pf0     = sum(g_v*f_h(E, E_F0))*energy_step
pf      = sum(g_v*f_h(E, E_F))*energy_step
Delta_P = sum(g_total*(f_h(E, E_F)-f_h(E, E_F0)))*energy_step
q       = Delta_P - pf0  (voltage_density="reference")
q       = Delta_P        (voltage_density="injected")
V       = e*L**2*q/(eps*(1-gamma)*(2-gamma))
J       = e*mu0*(2-gamma)*V*pf/L
```

The `reference` option subtracts equilibrium free carriers from the injected
total density. At E_F0, Delta_P is zero, so `injected` gives zero voltage while
`reference` gives a negative hole-branch voltage. These are explicit model
conventions with different boundary behavior. Their difference is relatively
small only when `abs(Delta_P)` is much greater than pf0. Neither choice resolves
the physical density interpretation automatically.

### Numerical stability

The integrated baseline can be of order 5e27 m^-3, whereas injected changes of interest
can be around 1e16 m^-3 or smaller. Adjacent float64 values at the baseline are
about 1.1e12 apart. Subtracting independently integrated totals loses precision,
especially close to equilibrium. FFT convolution of these totals can also lose
small carrier contributions and should not be assumed numerically equivalent.

For a stable calculation of the same mathematical difference, integrate
`g_total*(f_h(E,E_F)-f_h(E,E_F0))` before summing. Near saturated occupations,
evaluate that difference using the small complementary electron occupations
where appropriate, avoiding subtraction of two numbers near one. Validate such
rearrangements using equilibrium checks and quadrature convergence.

### Gamma and the displayed branch

The SI M4 prescription on p. 7 is `T/(T+T_t)` for T_t >= T and 0.5 otherwise.
The example configuration uses T_t/T. Choose gamma explicitly based on the
intended transport model and record the choice. A comparison at fixed parameters
cannot settle the validity of competing prescriptions.

Require `0 < gamma < 1` for an injection-model parameter choice. T_t/T fails
this condition at T_t >= T; it is not a general replacement for SI M4.
The model's prescribed gamma need not equal the inverse slope of the generated
J-V curve. That discrepancy should be inspected, not hidden by a round-trip test.

The hole-injection branch moves E_F downward from E_F0 toward E_v. Retain the
full grid for inspection and select the relevant positive-V, positive-J
branch for plots. Inspect density and Theta flags separately. M5 is manual
parameter adjustment; changing mu0 scales current for otherwise fixed inputs.

## 5. Notebook sequence

Work through these frontend stages together, checking each result before adding
the next. Implement the required reusable function in the backend and call it
from a short notebook cell; do not implement its algorithm in the notebook:

1. Imports and editable run configuration: input path, column/unit mapping,
   material/device parameters, model parameters, and numerical settings.
2. Call the backend loader for U, I, and T. Add backend conversion functions
   for V, J, and kelvin temperature when needed; plot the returned measurements
   in the notebook while preserving signs.
3. Calculate the analysis path using declared slope and density conventions.
4. Calculate and plot the band and trap DOS for the selected material; check
   N_v, N_c, units, and grid convergence.
5. Calculate occupations and the selected model branch from the run settings.
6. Overlay modeled and measured J-V curves, then density and mobility plots.
7. Adjust parameters by hand or compare a convention, with the effect visible.

Use the supplied example as the first worked case. Changing a dataset and its
parameters should require configuration edits, not changes to calculation code.
Accept arbitrary measurement row counts. The equations still need to match the selected transport regime.

The frontend/backend separation applies now; a larger package layout, CLI, and
automated fitting can wait. Record conventions beside the results.

## 6. Differences and open questions

Resolve these questions through physical interpretation and numerical checks:

| Topic | Question |
| --- | --- |
| Analysis gamma | Sensitivity of adjacent slopes to noise and signs |
| Model gamma | SI M4 versus other prescribed choices and generated curve slopes |
| Trap DOS | S5 or S6, with both normalized to N_t |
| Eq. (6) density | Absolute, trapped, and injected density interpretation |
| Model voltage | Density convention and equilibrium boundary behavior |
| Analysis E_F | Consistent temperature and nondegenerate validity regime |
| Energy grid | Convergence of integrals and interpolation of J(V) |

Published parameters also differ: Table 1 / Table S4 give the MAPbBr3 dark
values mu0 = 5.0e-3 m^2/V/s, N_t = 3.03e16 m^-3, E_t = -4.768 eV,
T_t = 8 K, E_F0 = -4.828 eV, and T = 300 K. Do not substitute these into a
comparison without labeling the change, including device dimensions.

Table S4's stated N_c of 4.52e16 cm^-3 is inconsistent with S3 evaluated with
m_eff_e = 0.32 and T = 299 K (approximately 4.52e18 cm^-3). Compute N_c and N_v from
their equations and report the discrepancy rather than copying that table value.

Earlier versions of this guide reported specific current ratios, grid errors,
and an unconverged gamma iteration without retaining a reproducible calculation.
Treat those numbers as unverified historical observations. Recompute them in the
notebook before using them to choose physics or impose pass/fail thresholds.
A better fit with one fixed parameter set supports that particular combination;
it does not prove a general equation or explain why earlier work was abandoned.

## 7. Validation

Separate general checks of the reusable tool from example regression checks.
Check backend functions independently of Jupyter, and verify that notebook
configuration and function calls work together.
For general usage, verify that changing the input path and parameter configuration
requires no code changes, that column/unit handling is explicit, and that the
calculation accepts different row counts.

The XLSX is exclusively an optional test input. Cached values represent its
saved state, and reading them does not recalculate formulas. Keep any reader,
cell references, or comparison-specific logic in test code, outside application
modules and notebook cells.

Report calculated values, cached values, and absolute/relative differences,
with parameter settings and known formula differences. Use the same full-precision
application constants during these tests. Small differences are expected;
exact matching and bitwise equality are not objectives. The user decides what
is too much. Do not invent an acceptance threshold, loosen a threshold, or
change application behavior to force agreement. Any comparison pass/fail limit
must be supplied by the user.

The cached model uses a different trap profile from the supported S5/S6 profiles,
so comparisons also need to distinguish formula differences from numerical error.
Comparison results do not establish physical validity.

Additional checks when the corresponding stage exists:

- Eq. (2)'s coefficient at gamma = 0.5 equals 9/8. Gamma = 0.5 alone does not
  imply mu_eff = mu0; that additionally requires Theta approaching one.
- Both supported trap integrals converge toward N_t (S5/S6).
- An equilibrium-subtracted density is zero at E_F0 within numerical tolerance.
- An algebraic forward/inverse check uses the **same gamma and density
  definitions** in both directions, away from singular points. It cannot be
  assumed across different definitions or an analysis that re-estimates gamma.
- Compare model and data over the measured voltage interval, declaring branch,
  interpolation, exclusions, and parameter settings. Agreement with measurements
  is a separate assessment from comparison-test agreement.

Numerical tests of mathematical properties and convergence remain separate from
XLSX comparisons. Their numerical tolerances do not set the acceptable discrepancy
from the XLSX; that decision belongs to the user.


