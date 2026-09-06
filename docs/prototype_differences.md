# Differences from the prototype spreadsheet

The model was first written as an Excel workbook, `SCLCKopecky.xlsx`. That
workbook was a prototype: it established the method and produced the article's
figures, but it is not the specification. **The published equations are** — the
article and its Supplementary Information — and this port follows them.

The prototype remains in `data/` for one purpose: it is the project's only
independent numerical oracle. Every other check in the suite is an analytic
identity or an internal consistency relation, and those cannot catch a wrong
assumption that the code and its own docstrings share. The regression tests
therefore still compare against its cached cell values, with its options
selected explicitly.

To run the two side by side:

```bash
python asclc.py --params params/MAPbBr3_S2.toml   --plot figures/current/
python asclc.py --params params/prototype_S2.toml --plot figures/prototype/
```

---

## 1. Where the prototype and the published equations disagree

Each is a switchable option. The default now follows the published equation; the
prototype's reading is still available.

### 1.1 Trap density of states — the largest difference

The prototype's trap column, `g(E)!M13`, reads

    ... / (1 + EXP(u)^2)

Excel binds `^` more tightly than `+`, so this evaluates `1 + e^{2u}`, not
`(1 + e^u)^2`. With its `N_t/(2 k_B T_t)` prefactor (`g(E)!M7`) the result is a
`sech(u)/2` profile.

| | Expression | ∫ g_t dE |
|---|---|---|
| Eq (S5) as printed — `SI_BIEXPONENTIAL`, **default** | `(N_t/k_BT_t)·e^u/(1+e^u)^2` = `sech²(u/2)/4` | `N_t` |
| Prototype — `PROTOTYPE_SECH` | `(N_t/2k_BT_t)·e^u/(1+e^{2u})` = `sech(u)/4k_BT_t` | `(π/4)·N_t` |

The π/4 is exact: `∫sech = π` against `∫sech² = 2`. The prefactor normalizes
neither form — under Eq (S5) it would give `N_t/2`.

**Why it matters.** Under the prototype profile `N_t` is a scale factor, not a
concentration: the trapped population saturates at 0.785·`N_t`, so `N_t` read off
a saturation plot (Supplementary Note A6) is 27 % low. Under the default it
saturates at exactly `N_t`, which is what makes step A6 self-consistent.

The two profiles share a peak height and an `e^{-|u|}` tail, so they are
indistinguishable on a log-scale DOS plot. The difference is entirely in the
area.

### 1.2 γ in the model branch

Supplementary Note M4 gives, as printed:

    γ = T/(T_t + T)   for T_t ≥ T,      γ = 0.5   for T_t < T

`T/(T_t+T)` is `1/(1 + T_t/T)`, which is `1/m` for the Mark–Helfrich
trap-filled-limit exponent `m = 1 + T_t/T` of an exponential trap distribution.
That reconciles Note M4 with the article's own definition `γ = 1/m`.

The prototype does not use it. Cell `j(U)!C6 = MODEL!B30/MODEL!B9` computes
`T_t/T`, which is the excess exponent `m − 1`, not `1/m`.

| `GammaModel` | Expression | At T_t = 30 K, T = 299 K |
|---|---|---|
| `SI_NOTE_M4` — **default** | `T/(T_t+T)` if `T_t ≥ T` else `0.5` | 0.5 |
| `TT_OVER_T` — prototype | `T_t/T` | 0.1003 |
| `TT_OVER_T_PLUS_TT` | `T_t/(T+T_t)` — in no source | 0.0912 |

Every configuration in the article and the prototype has `T_t < T`, so the
default selects the constant branch: **γ = 0.5 throughout**.

**Why it matters less than it looks.** γ is not identifiable from J-V data. It
enters only through `1/[(1−γ)(2−γ)]` in V and `(2−γ)` in J, neither depending on
`E_F`, so changing it rescales both axes by constants — on log-log a rigid
translation that leaves the curve's shape untouched. It is degenerate with `μ_0`
and the other scale parameters. Do not fit it.

### 1.3 Θ

Eq (S12), printed identically in the SI and the article, is

    Θ = p_f/p_s = p_f/(p_f + p_t)

The prototype computes `n(E)!N = ABS(L/J)`, which is the absolute free
population over the **injected** total — a mixed reference, and unbounded.

| `ThetaModel` | Reading | Bounded by 1? |
|---|---|---|
| `ABSOLUTE_OVER_TOTAL` — **default**, Eq (S12) | `p_f/(p_f + p_t)` | yes |
| `ABSOLUTE_OVER_INJECTED` — prototype | `p_f/p_inj` | **no** |
| `INJECTED_OVER_INJECTED` | `Δp_f/p_inj` | yes |

They agree at high injection and differ by up to a factor of six through the
trap-filling region, which is where the physics of interest sits. Under the
prototype reading Θ exceeds 1 near zero bias — so `μ_eff = μ_0Θ` exceeds the
microscopic mobility — and the `valid` mask flags it. Under the default Θ is a
fraction by construction and `valid` never fires on it.

---

## 2. Numerical defects the port corrects

These are not modelling choices and are not switchable.

### 2.1 Band quadrature

The prototype integrates with the rectangle rule on a uniform 3 meV grid. The
√ band edge has an infinite derivative, where that converges only as O(h^{3/2}),
giving about **1.5–3 % error** on the band integrals. The port substitutes
`E = E_v − t²`, absorbing the singularity into the Jacobian and leaving a smooth
integrand. The port's values are the more accurate ones; the regression
tolerances are set by the prototype's error, not the port's.

### 2.2 The slope window is degenerate

`Data-calculations!I1 = 0` is the binning width for every `AVERAGE(INDIRECT(...))`
and `LINEST(INDIRECT(...))` on the sheet. At zero, each `LINEST` regresses one
point against one point and **returns γ = 0 for all 210 rows**.

Since `(1−γ)(2−γ)` is maximal at γ = 0, the experimental `p_t` column collapses
to `2ε₀ε_r V/(eL²)` — the upper envelope of Eq (6), carrying no current data at
all. `p_f` and `μ_eff` still respond to the measurement but are biased low by
`(2−γ)/2` and `(1−γ)(2−γ)²/4`.

The port requires a window of at least 3 and defaults to 7. **The prototype's
entire analysis branch is therefore not reproducible, and is not a regression
target.** Only its model-branch columns are.

### 2.3 Direction of the slope fit

The prototype computes `Data-calculations!K = LINEST(ln V, ln j)`, regressing
`ln V` on `ln j` to get γ directly. The port regresses `ln|J|` on `ln|V|` and
takes `γ = 1/m`.

The article defines the two as equal — `γ = 1/m = d ln V/d ln J` — and they are
for a clean power law, but ordinary least squares is not symmetric: with scatter
`1/slope(y|x) ≠ slope(x|y)`, differing by more than 5 % at an 11-point window on
the bundled measurement. The J direction is kept because it stays conditioned
through the trap-filled-limit region, where `ln V` is nearly constant and the
prototype's direction regresses against a near-degenerate abscissa.

### 2.4 Temperature

Eq (7) uses the measured temperature of each data point. The prototype takes
`k_BT` per row from its recorded column (`Data-calculations!N`) but leaves `N_v`
at the nominal value. The port evaluates both at the same per-point temperature.

Worth about 0.002 eV in `E_F` over the bundled measurement's 282–314 K range —
small beside the 0.075 eV the `k_BT` term itself contributes, but there is no
reason to reproduce the inconsistency.

### 2.5 Physical constants

The prototype hard-codes rounded values (`e = 1.602e-19`, `k_B = 1.38e-23`,
`h = 6.62607004e-34`, `ε₀ = 8.854e-12`, `m₀ = 9.109e-31`). The port uses
CODATA 2018. The largest relative difference is `k_B` at 4.7e-4, which
propagates to about the same relative shift in every Boltzmann factor.
`Constants.prototype()`, `--prototype-constants` and `constants = "prototype"`
reproduce the rounded set.

---

## 3. What the port adds

- **A `valid` mask** on `AnalysisResult`, `ModelCurve` and `CarrierDensities`,
  meaning the same thing in each: this point's outputs are physically
  meaningful. The prototype has no such concept, and roughly half the points of
  a real measurement fail it while looking like ordinary numbers.
- **Exact charge conservation.** `n_inj = −p_inj` bit for bit. Computing the
  electron total as `(n_f + n_t)` minus its equilibrium value omits the
  valence-band term, which is 136 times larger than what remains.
- **Exactly zero injected charge at `E_F0`.** The equilibrium reference is
  evaluated inside the same array as the requested Fermi levels, so
  `p_s_injected[0]` is 0.0 rather than a small residual of either sign.
- **Strict parameter files.** Unknown keys are an error; a misspelled `E_t`
  would otherwise leave the model on its default while the file appeared to say
  otherwise.

---

## 4. Prototype quirks that are not differences

Recorded so anyone reading the spreadsheet is not misled by them.

- **`n(E)!O = ps_injected − pf0`.** The voltage column reads from `O`, which
  subtracts the equilibrium free-hole population from the injected total. At
  `pf0 ≈ 1.4e12 m⁻³` against injected totals of 1e19–1e27, this is fifteen orders
  down and numerically irrelevant. The port does not reproduce it.
- **Column labels.** `n(E)`'s `ps` and `ns` are injected totals over the full
  density of states, not `p_f + p_t`. `p_t` in Eq (6) is the total space charge,
  not the trapped part. `ptm` in Figs. 4c,d and 6 is the injected total —
  `ModelCurve.p_space_charge`, not `ModelCurve.p_t`.
- **Dead names.** All 34 defined names in the workbook are `OFFSET(#REF!,...)`,
  each duplicated four times.
- **Supplementary Table S4's `N_v`.** The tabulated values do not follow from
  Eq (S4), which gives ≈ 4.2e18 cm⁻³ for MAPbBr₃ against the table's
  4.52e16 cm⁻³; the MAPbI₃ entry of 5.17e10 cm⁻³ is six orders adrift. The
  formula agrees with the prototype and with the standard
  `2.5e19 cm⁻³ (m*/m₀)^{3/2}`. The port follows the formula.

---

## 5. Net effect on the numbers

Model branch, at the article's published parameter sets, prototype options
versus the current defaults:

| Configuration | γ | J max (A m⁻²) | p_t max |
|---|---|---|---|
| MAPbBr₃ S2 | 0.1003 → 0.5 | 0.752 → 0.157 | 0.785·N_t → **N_t** |
| MAPbBr₃ dark | 0.0267 → 0.5 | 0.0624 → 0.0100 | 0.785·N_t → **N_t** |
| MAPbBr₃ light | 0.1333 → 0.5 | 0.183 → 0.0618 | 0.785·N_t → **N_t** |
| MAPbI₃ light | 0.3167 → 0.5 | 0.00228 → 0.00054 | — (V_max truncates) |

The J shift is the γ change and is a rigid translation on log-log, absorbed by
`μ_0` on refitting. The `p_t` change is the substantive one: trapped charge now
saturates at `N_t` rather than 0.785·`N_t`.
