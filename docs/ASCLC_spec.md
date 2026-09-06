# A-SCLC model specification

Specification for `asclc.py`, a Python implementation of the advanced space-charge-limited current
model of Gavranovic, Zmeskal, Weiter and Pospisil, *Communications Physics* **8**, 280 (2025),
ported from the reference Excel workbook `SCLCKopecky.xlsx`.

Validated against the workbook's own computed values and the article's published parameter sets;
see section 6. Three items in section 7 are unresolved in the sources rather than in the port: the
trap distribution (7.1), the model-branch γ, which is not identifiable from J-V data (7.2), and which
of three readings of Eq (S12) defines Θ (7.3).

---

## 1. Structure

The model has two independent branches. They share equations but run in opposite directions and meet
only when results are plotted together.

**Analysis branch** — `analyse_jv`, steps A1–A7 of Supplementary Note 2. Voltage is the independent
variable. From a measured J-V curve the local logarithmic slope gives `m = d ln J / d ln V` and
`γ = 1/m`; effective mobility, carrier concentrations and Fermi level follow algebraically from
Eqs (4)–(7).

**Model branch** — `model_curve`, steps M1–M5. **The Fermi level is the independent variable.** For
each `E_F` the occupation integrals give the free and injected concentrations; the voltage follows by
inverting Eq (6) and the current from Eq (S14):

    V(E_F) = e L² p_inj(E_F) / (ε₀ ε_r (1−γ)(2−γ))
    J(E_F) = e μ₀ (2−γ) p_f(E_F) V(E_F) / L
    Θ(E_F) = p_f(E_F) / p_inj(E_F)
    μ_eff  = μ₀ Θ

No self-consistent solve appears anywhere: the J-V curve is generated parametrically in `E_F`.

Both branches are pure functions of their inputs, so wrapping an optimiser around either requires no
change to the module.

---

## 2. Equations

Energies in eV, concentrations m⁻³, mobilities m² V⁻¹ s⁻¹, lengths m, current density A m⁻².
Constants default to CODATA 2018; `Constants.workbook()` reproduces the workbook's rounded values.

### Density of states

Square-root transport bands, Eqs (S1), (S2):

    g_e(E) = C_e √(E − E_c),  E ≥ E_c        C = 4π (2 m₀ m* e)^{3/2} / h³
    g_h(E) = C_h √(E_v − E),  E ≤ E_v

Effective densities of states, Eqs (S3), (S4):

    N_{c,v} = 2 (2π m* m₀ k_B T / h²)^{3/2}

These are mutually consistent: `∫ g_h(E)[1−f] dE → N_v exp(−ΔE_F/k_BT)` in the Boltzmann limit, which
the tests check as an independent validation of both.

Localized states, Eq (S5), biexponential:

    g_t(E) = (N_t / k_B T_t) · e^u / (1 + e^u)²,    u = (E − E_t)/(k_B T_t)

This integrates to exactly `N_t`. Eq (S6) offers a Gaussian alternative with `σ = 2 k_B T_t`. The
workbook implements neither; see 7.1.

### Occupation

Fermi–Dirac, Eq (S9), and its hole counterpart:

    f(E − E_F)     = 1 / (1 + exp((E − E_F)/k_B T))
    1 − f(E − E_F) = 1 / (1 + exp(−(E − E_F)/k_B T))

Free concentrations integrate their own band only, Eqs (S7), (S10):

    p_f(E_F) = ∫ g_h(E) [1 − f(E − E_F)] dE
    n_f(E_F) = ∫ g_e(E) f(E − E_F) dE

### Injected charge

Eqs (S8), (S11) give the total concentration as one integral referenced to `E_F0`. The reference is
an **equilibrium** reference, not a range of integration, and it runs over the **full** density of
states — valence band, conduction band and traps together:

    p_inj(E_F) = ∫ g_total(E) [f(E − E_F0) − f(E − E_F)] dE
    n_inj(E_F) = − p_inj(E_F)

The second line is charge conservation: an injected hole is an electron removed from the same states.
It holds exactly and the implementation enforces it.

`p_inj` vanishes at `E_F = E_F0`, which zero bias requires, and it is the space charge driving Eq (6).
The cross-band terms matter asymmetrically: the conduction band contributes exactly zero to `p_inj`,
being empty at both Fermi levels, while the valence band dominates `n_inj`, exceeding the conduction
and trap terms by two orders of magnitude at high injection.

### Derived quantities

    Θ     = p_f / p_inj                                          (S12)
    μ_eff = μ₀ Θ = L³ J / (ε₀ ε_r (1−γ)(2−γ)² V²)                (4)
    p_f   = L J / (e μ₀ (2−γ) V)                                 (5)
    p_t   = ε₀ ε_r (1−γ)(2−γ) V / (e L²)                         (6)
    p_f   = N_v exp(−ΔE_F/k_B T),  ΔE_F = E_F − E_v              (7)
    E_F   = E_v + k_B T ln(N_v / p_f)                            (7, inverted)

Three points that are easy to get wrong, each pinned by a test:

- **Eqs (2) and (4) carry `(1−γ)(2−γ)²`**, Eq (6) the first power. The Mott–Gurney limit fixes the
  grouping: at γ = ½, `(1−γ)(2−γ)² = 9/8`, recovering Eq (3).
- **Θ divides by the injected total**, not by `p_f + p_t`. In the analysis branch this makes
  `Θ = μ_eff/μ₀ = p_f/p_t` an exact identity from Eqs (4)–(6), where `p_t` from Eq (6) is the total
  space charge despite its label. Which reading Eq (S12) intends is unresolved; see 7.3.
- **Eq (7) inverts to `ln(N_v/p_f)`**, since `ΔE_F ≥ 0` with `E_F` above `E_v`. Eq (7) is the
  Boltzmann limit, so where `p_f > N_v` the material is degenerate and the extraction does not apply;
  the code warns.

### γ

    Analysis branch:  γ = 1/m from the local slope of the data
    Model branch:     γ = T_t/T                          (workbook, j(U)!C6)

The model branch value is a positioning convention rather than a derived quantity; see 7.2.

---

## 3. Temperature

`analyse_jv` accepts `temperature` as a scalar or one value per data point, used in Eq (7).

Pass the measured per-point temperature whenever it is recorded. In the reference dataset the
recorded temperature spans 282–314 K, moving the extracted `E_F` by 0.075 eV. The Fermi level shifts
the article reports are 0.046 eV and 0.006 eV, so the temperature scatter exceeds the effect being
measured; holding T at a nominal value reassigns that scatter to the physics.

`N_v` and `k_B T` are evaluated at the same temperature. The workbook takes `k_B T` per row from the
measured column while leaving `N_v` at the nominal value, worth about 0.002 eV over its own range.

---

## 4. Validity

`AnalysisResult`, `ModelCurve` and `CarrierDensities` each carry a `valid` mask with the same
meaning: **this point's outputs are physically meaningful. Read every other field through it.**

A point is invalid when:

- the current changes sign inside its fit window, so the log-log slope is not a slope of anything
  (analysis branch only);
- `γ` falls outside (0, 1), where Eqs (4)–(6) are out of domain and return negative concentrations
  (analysis branch only);
- `Θ > 1`, since `μ_eff = μ₀Θ` cannot exceed the microscopic mobility (both branches).

The first two occur throughout the ohmic region, where the measured current is noise straddling zero:
on the reference measurement 114 of 214 points fail one of the three, and ungated `Θ` ranges over
−0.13 to 1.14 with nothing marking which values mean anything. The third has a different cause in
each branch — in the analysis branch it means the supplied `μ₀` is too small, in the model branch that
the equilibrium free population still rivals the injected charge — but the same consequence, so it is
flagged the same way.

Invalid points are flagged rather than blanked: plotting what was rejected is often the quickest way
to judge whether a measurement is usable. `analyse_jv` warns when more than half the points fail,
when none survive, and when `Θ > 1` or degeneracy affects otherwise-fittable points.

---

## 5. Numerics

**Band quadrature.** The √ band edge has an infinite derivative, where the trapezoid rule converges
only as O(h³ᐟ²); on a uniform 1 meV mesh the error in `p_f` reaches 2e-3 relative, larger than the
Fermi–Dirac correction the integral exists to resolve. Substituting `E = E_v − t²` absorbs the
singularity into the Jacobian:

    ∫ C √(E_v − E) [1−f] dE  =  ∫ 2 C t² [1−f] dt

leaving a smooth integrand that converges at a couple of thousand points.

**Grid size.** Defaults are `band_points=2001`, `trap_points=1001`, converged to about 1e-5. Larger
grids are slower without being more accurate: past roughly 4000 points the intermediate matrices
exceed cache and evaluation slows by nearly an order of magnitude.

**Conditioning of the injected charge.** The occupancy difference is taken inside the integral.
Integrating each band's absolute occupancy and subtracting afterwards differences far-band numbers of
order 1e28 whose ulp is near 1e12 — comparable to the result — making the answer depend on chunk
size.

**Memory.** The occupation integrals are matrix–vector products, chunked so peak memory stays bounded
regardless of sweep length. Results are chunk-independent to 2e-15.

**Cost.** `model_curve` at 3001 points runs in about 0.5 s, fast enough to tune against
interactively.

---

## 6. Validation

**Published parameter sets.** All four configurations of Supplementary Table S4 — both materials,
dark and illuminated — produce curves spanning the measured 0–3 V range with finite current
throughout, and Θ ≤ 1 across the in-domain region. Saturated space charge lands within a factor of
0.64–1.67 of the `p_t` values in Table 1.

**Workbook.** Against its own computed values for the MAPbBr₃ S2 configuration:

| Quantity | Workbook cell | Agreement |
|---|---|---|
| `p_f` | `n(E)!L` | 1.5–2 % |
| `p_inj` | `n(E)!O` | 0.07–1.5 % |
| `n_f` | `n(E)!E` | 3 % |
| `n_inj` | `n(E)!H` | 0.02–1.5 % |
| `V` | `j(U)!C` | 0.1–1.5 % |
| `J` | `j(U)!E` | 2–3 % |
| `Θ`, `μ_eff` | `n(E)!N`, `j(U)!K` | 1–2 % |
| analysis branch `μ_eff`, `p_f`, `p_t`, `Θ` | `MODEL!I,K,M,Q` | 1e-4 |

The residual on model-branch quantities is the workbook's rectangle rule on a 3 meV grid across the
√ band edge; the values here are the more accurate ones. The analysis branch agrees to rounding
because both evaluate the same closed-form expressions.

A round trip — model branch generates (V, J), analysis branch reads it back — recovers `p_f` to 4 %
and `E_F` to 1 meV, the residual being the γ definition mismatch of 7.2.

---

## 7. Open items

### 7.1 Trap density of states

The workbook's trap column reads `.../(1 + EXP(u)^2)`. Excel binds `^` before `+`, so this evaluates
`1 + e^{2u}` rather than `(1 + e^u)²`, giving a `sech(u)/2` profile where Eq (S5) specifies
`sech²(u/2)/4`. With the workbook's `N_t/(2 k_B T_t)` prefactor:

| | ∫ g_t dE |
|---|---|
| Eq (S5) as printed | `N_t` |
| Workbook | `(π/4) N_t` = 0.7854 `N_t` |

The π/4 is exact — `∫sech = π` against `∫sech² = 2` — and reproduces in the workbook's own cached
values. The prefactor normalizes neither form: with Eq (S5) it would give `N_t/2`.

The two profiles share a peak height and an `e^{-|u|}` tail, so they are indistinguishable on a
log-scale DOS plot; the difference is entirely in the area.

**Recommendation for new work: `TrapProfile.SI_BIEXPONENTIAL`.** It is what the SI documents and the
only variant whose integral is `N_t`, which is what makes `N_t` a concentration rather than a scale
factor. The default is `WORKBOOK_SECH` so the port reproduces existing curves; that default is for
validation.

**To decide empirically:** refit `N_t` under both profiles across dark/light or a temperature series.
The discriminating data is in the shoulders of `p_t(E_F)`, not at its peak.

### 7.2 γ in the model branch

`T_t/T` (workbook) or `T_t/(T + T_t)` (the most plausible reading of Note M4, whose PDF text is not
legible enough to be certain). At `T_t = 30 K, T = 299 K`: 0.1003 versus 0.0912.

**This is not identifiable from J-V data.** γ enters only through `1/[(1−γ)(2−γ)]` in V and `(2−γ)`
in J, neither depending on `E_F`. Changing it rescales both axes by constants: on log-log a rigid
translation, leaving the curve's shape untouched to 1e-12. The 1.5 % shift is degenerate with μ₀ and
the other scale parameters. A temperature series does not separate them either, since μ₀(T) is a free
function that absorbs the difference; that needs an independent constraint on μ₀.

**The more useful observation:** the analysis branch *defines* γ = 1/m from the local slope, while the
model branch assumes a constant. The curve the model branch produces has a median slope of about
m = 17, implying γ ≈ 0.06 against the 0.100 it is given. The two branches use the symbol for different
things. To make them consistent, use `γ(E_F) = 1/m(E_F)` from the model curve's own slope — a small
fixed-point iteration, left as a modelling decision.

### 7.3 Which Θ does Eq (S12) mean?

Three readings are defensible, they agree at high injection, and they differ by up to a factor of six
through the trap-filling region — which is where the physics of interest sits.

| Reading | Bounded by 1? | Used by |
|---|---|---|
| absolute `p_f / (p_f + p_t)` | yes | Eq (S12) as printed |
| absolute `p_f / p_inj` | **no** | the workbook, and this port |
| injected `Δp_f / p_inj` | yes | neither |

The port follows the workbook, which also keeps the two branches consistent: the analysis branch's
`Θ = μ_eff/μ₀ = p_f/p_t` likewise puts an absolute `p_f` from Eq (5) over the Eq (6) space charge.

That form is unbounded. Where the equilibrium free population rivals the injected charge, Θ exceeds
1 and `μ_eff = μ₀Θ` exceeds the microscopic mobility, which cannot happen. The workbook stays under 1
for its own parameters but the article's MAPbBr₃ illuminated set reaches Θ = 1.41 immediately above
`E_F0`.

This is treated as a domain limit rather than hidden by redefinition: Θ > 1 clears the `valid` mask
in both branches (section 4). In the model branch it is false only within a few millivolts of zero
bias, so whether any grid point lands there depends on the sampling. `p_f`, `p_t` and `p_s_injected`
are all returned, so the other two readings are one division away.

### 7.4 Naming hazards in the sources

Several quantities are labelled as one thing and computed as another. These account for most of the
subtle discrepancies between the workbook, the article and a naive reading:

| Label | What it is |
|---|---|
| `p_t` in Eq (6) | total space charge, not the trapped part |
| `ps`, `ns` in `n(E)` | injected totals over the full DOS, not `p_f + p_t` |
| `ptm` in Figs 4c,d and 6 | the injected total, not `p_t` |
| `nt` in `n(E)!H` | injected electron total, equal to `−p_inj` |

In the port, `ModelCurve.p_space_charge` is the article's `ptm`; `ModelCurve.p_t` is the absolute
trapped population, which starts at 0.54·`N_t` rather than zero.

---

## 8. Notes on the reference workbook

Findings that affect how the workbook's own outputs should be read.

**`Data-calculations!I1 = 0`.** This is the binning width for every `AVERAGE(INDIRECT(...))` and
`LINEST(INDIRECT(...))` on the sheet. At zero, each `LINEST` regresses one point against one point
and returns γ = 0 for all 210 rows. Since `(1−γ)(2−γ)` is maximal at γ = 0, the experimental `p_t`
column collapses to `2ε₀ε_r V/(eL²)` — the upper envelope of Eq (6), carrying no current data. `p_f`
and `μ_eff` still respond to the measurement but are biased low by `(2−γ)/2` and `(1−γ)(2−γ)²/4`.
Set a window before using the analysis branch; `local_loglog_slope` requires at least 3 points.

**Configuration.** MAPbBr₃ sample "S2", dark, measured 2024-09-27: L = 0.6 mm, T = 299 K, T_t = 30 K,
E_F0 = −4.84 eV, E_t = −4.82 eV, N_t = 4.7e16 m⁻³, μ₀ = 0.0027 m² V⁻¹ s⁻¹. Neither sample reported in
the article, so Table 1 is not a validation target.

**Supplementary Table S4 `N_v`.** The tabulated values do not follow from Eq (S4), which gives
≈ 4.2e18 cm⁻³ for MAPbBr₃ against the table's 4.52e16 cm⁻³; the MAPbI₃ entry of 5.17e10 cm⁻³ is six
orders adrift. The formula agrees with the workbook and with the standard
`2.5e19 cm⁻³ (m*/m₀)^{3/2}`. The port follows the formula.

**Dead names.** All 34 defined names are `OFFSET(#REF!,...)`, each duplicated four times.

---

## 9. API

    Constants, Material, Device, ModelParams    inputs
    EnergyGrid                                  quadrature configuration
    TrapProfile, GammaModel, SpaceCharge        modelling choices
    trap_dos, valence_band_dos,                 density of states
      conduction_band_dos, effective_dos
    fermi_dirac, carrier_densities              occupation
    model_curve  -> ModelCurve                  model branch
    analyse_jv   -> AnalysisResult              analysis branch
    local_loglog_slope                          log-log slope
    load_jv      -> Measurement                 measurement files

`python asclc.py` runs the model branch on the workbook configuration and prints a summary.

    --data PATH      also load and analyse a measurement file
    --celsius        its temperature column is in degrees Celsius
    --bin N          average N consecutive points
    --window N       points per local slope fit
    --csv PATH       write the model curve
    --trap-profile, --gamma-model, --v-max, --workbook-constants

---

## 10. Not implemented

**Extraction model.** Supplementary Note 4 describes an extraction counterpart to the injection
model, needed for the MAPbI₃ dark curve where an injection barrier makes extraction dominant. The
model branch sweeps `E_F` downward from `E_F0` only.

**Parameter fitting.** The five parameters are set by hand. Both branches are pure functions, so a
residual wrapped around `model_curve` and handed to `scipy.optimize.least_squares` is the natural next
step. Note the degeneracy of 7.2 when choosing a parameterisation.
