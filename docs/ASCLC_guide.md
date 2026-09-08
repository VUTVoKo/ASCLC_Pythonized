# A-SCLC: a guide to building the script

How to implement the advanced space-charge-limited-current (A-SCLC) model of
Gavranovic, Zmeskala, Weiter & Pospisil, *Commun. Phys.* **8**, 280 (2025),
following the algorithm of Supplementary Note 2 (steps A1–A7 and M1–M5).

Everything below was read out of the three sources in `data/`:

| Source | What it is | What it settles |
| --- | --- | --- |
| `s42005025022021.pdf` | main article, 10 pp. | Eqs (1)–(7), the physical meaning of Θ and γ |
| `42005_2025_2202_MOESM2_ESM.pdf` | Supplementary Information, 23 pp. | Eqs (S1)–(S16), Tables S2/S4, **the A1–A7 / M1–M5 algorithm** (Note 2, Fig. S7) |
| `SCLCKopecky.xlsx` | the authors' working workbook | the *executable* reference: exact formulas, grids, constants |

**No source is authoritative on its own.** The PDFs are the considered account
but contain at least one equation that contradicts the others and one table
entry off by 10²; the workbook runs, but it is a working file — simplified in
places, and left in a half-finished state in others (its saved analysis branch
was computed with γ ≡ 0; see §3, A2). So each disagreement is adjudicated on
evidence in §6, and where the evidence is numerical it is quoted. The workbook
remains the numerical oracle for regression tests (§7) whatever you decide.

The workbook's live case is **MAPbBr₃, sample S2, dark, T = 299 K** — which is
what `examples/` was produced from.

---

## 0. Scope

Two independent branches share one parameter set. Keep them separate in code.

```
                 ┌────────────────────────────────────────────────┐
  VJT.dat  ───▶  │  A-branch  "analysis"   (A1–A7)                │ ──▶ one row per
 (V, I, T)       │  measured J-V → γ, μ_eff, p_f, p_t, Θ, E_F     │     measured point
                 └────────────────────────────────────────────────┘
                 ┌────────────────────────────────────────────────┐
 .params  ───▶   │  M-branch  "modeling"   (M1–M5)                │ ──▶ one row per
 (5 knobs)       │  E_F sweep → g(E) ⊛ f → p_f, p_t → V, J        │     E_F grid point
                 └────────────────────────────────────────────────┘
```

The A-branch inverts measured data through the SCLC equations. The M-branch
generates a curve from the density of states. They meet only when plotted
together. **Neither branch fits anything** — the five variables of Table S2 are
hand-selected and read from the parameter file. Automated fitting is out of scope.

---

## 1. Constants and unit conventions

Use the workbook's constants verbatim, not CODATA, or the last digits will not match:

```python
e    = 1.602e-19       # C      (workbook, not 1.602176634e-19)
k_B  = 1.38e-23        # J/K    (workbook, not 1.380649e-23)
eps0 = 8.854e-12       # F/m
h    = 6.62607004e-34  # J s
m_e  = 9.109e-31       # kg
T0   = 273.15          # K      (degC -> K offset)
```

Conventions throughout:

* **Energies in eV on the vacuum scale** (negative: `E_c = -3.36`, `E_v = -5.58`).
  Only `k_B T` appears, as `kT_eV = k_B*T/e`.
* **Density of states in m⁻³ eV⁻¹**, so `∫ g dE` comes out in m⁻³ directly when
  `dE` is in eV.
* Concentrations m⁻³, lengths m, mobility m² V⁻¹ s⁻¹, current density A m⁻².
* Hole ("p") quantities are what the MAPbBr₃ device needs; the electron ("n")
  quantities are the mirror image and are worth computing anyway — Fig. 4c,d of
  the article plots both.

---

## 2. The equations, resolved

The PDF text layer mangles sub/superscripts and drops Symbol-font glyphs. The
forms below were recovered glyph-by-glyph and cross-checked against the workbook
formulas; they are correct.

### Main text

| # | Equation | Note |
| --- | --- | --- |
| (1) | `J = e*mu0*p_f*V/L` | Ohm's law |
| (2) | `J = eps0*eps_r*mu0*Theta*(1-g)*(2-g)^2 * V^2/L^3` | SCLC with traps |
| (3) | `J = (9/8)*eps0*eps_r*mu0*Theta * V^2/L^3` | Mott–Gurney; is (2) at γ = ½ |
| (4) | `mu_eff = mu0*Theta = L^3*J / (eps0*eps_r*(1-g)*(2-g)^2*V^2)` | |
| (5) | `p_f = L*J / (e*mu0*(2-g)*V)` | note the `(2-γ)`, absent from (1) |
| (6) | `p_t = eps0*eps_r*(1-g)*(2-g)*V / (e*L^2)` | see §6.1 — this is really `p_s` |
| (7) | `p_f = N_v*exp(-dE_F/kT)`, `dE_F = E_F - E_v` | invert for `E_F` |
| (S12) | `Theta = p_f/p_s = p_f/(p_f + p_t)` | |
| (S14) | `J = e*mu0*p_f*(2-g)*V/L` | (5) rearranged; used by the M-branch |

`γ ≡ 1/m = dlnV/dlnJ`. Regions: γ = 1 ohmic, γ = ½ Mott–Gurney, γ → 0 trap-filled limit.

### Supplementary — density of states

```
(S1)  g_c(E) = (8*pi*me*/h^3) * sqrt(2*me* * (E - E_c))      for E >= E_c, else 0
(S2)  g_v(E) = (8*pi*mh*/h^3) * sqrt(2*mh* * (E_v - E))      for E <= E_v, else 0
(S3)  N_c = 2*(2*pi*me* *k_B*T/h^2)^(3/2)
(S4)  N_v = 2*(2*pi*mh* *k_B*T/h^2)^(3/2)
(S5)  g_t(E) = (N_t/kT_t) * e^u/(1 + e^u)^2,  u = (E-E_t)/kT_t      [biexponential]
(S6)  g_t(E) = (N_t/(s*sqrt(2pi))) * exp(-(E-E_t)^2/(2s^2)),  s = 2*k_B*T_t   [Gaussian]
```

(S5) is `(N_t/4kT_t)*sech²(u/2)` and integrates to exactly `N_t`. The workbook
implements something else — see §6.2.

With energies in eV the band DOS prefactor collapses to one constant:

```python
C_v = 4*pi*(2*m_e*m_h_rel*e)**1.5 / h**3    # m^-3 eV^-3/2
g_v = C_v * sqrt(E_v - E)                   # m^-3 eV^-1
```

### Supplementary — occupation

```
(S9)   f(E - E_F) = 1/(1 + exp((E - E_F)/kT))       electrons
       1 - f      = 1/(1 + exp((E_F - E)/kT))       holes
(S7)   n_f(E_F) = INT_{E_c}^{+inf}   g(E)*f dE
(S8)   n_t(E_F) = n_s - n_f = INT_{E_F0}^{E_c} g(E)*f dE
(S10)  p_f(E_F) = INT_{-inf}^{E_v}   g(E)*(1-f) dE
(S11)  p_t(E_F) = p_s - p_f = INT_{E_v}^{E_F0} g(E)*(1-f) dE
```

The `E_F0` integration limits in (S8)/(S11) are the article's way of saying
"count only what was *injected* relative to equilibrium". The workbook achieves
the same by subtracting the equilibrium value instead of moving the limit — see
§4.2. That subtraction is not optional (§4.2), but *what* the workbook subtracts
in one of the two places is wrong (§6.5).

### The γ of the model branch (M4)

The Supplementary states, unambiguously once the Symbol-font glyphs are decoded:

```
gamma = T/(T + T_t)    for T_t >= T
gamma = 0.5            for T_t <  T
```

The workbook instead hardcodes `gamma = T_t/T` (`'j(U)'!C6`). For the live case
(T_t = 30 K, T = 299 K) the two give **0.5 vs 0.1003**, and the resulting model
currents differ by ~1500× at 3 V. This is by far the most consequential
disagreement between the sources, and the evidence goes against the published
formula — see §6.3 before writing a line of the M-branch.

---

## 3. The A-branch (A1–A7): analysing measured data

One output row per measured point. Reference: `MODEL` sheet columns E–W, fed from
`Data-calculations` columns E, G, K, N.

**A1 — read the measurement.** `VJT.dat` is tab-separated, header `U(V) I(A) T(C)`.
Convert `V = U + V_min`, `J = I/S`, `T_K = T_C + 273.15`. Optionally average
`bin_size` consecutive rows first (the workbook's `Data-calculations!I1` bin
width; the shipped case uses 1, i.e. no binning — `examples/` keeps all 214 raw points).

**A2 — local slope.** `γ_i` = slope of `ln|V|` against `ln|J|` over a short window
centred on point *i* (`Data-calculations!K`, a `LINEST` over `I1+1` rows; a
symmetric `window`-point least-squares fit is the sane generalisation). `m = 1/γ`.

> **The workbook ships degenerate here.** `Data-calculations!I1 = 0`, so every
> `LINEST` window is a *single cell* and the saved `γ` column is **0 at every
> point** (so is the `Ea` column, for the same reason). Every downstream number
> in the workbook's analysis branch was therefore computed at `γ = 0`, i.e. with
> `(1-γ)(2-γ)² = 4` and `(2-γ) = 2`. Reproduce that only to hit the §7 targets;
> for real work the window must be ≥ 2 points.

This is the noisiest step: in the ohmic region the measured current sits on the
noise floor and γ is meaningless — which is why most low-V rows of
`examples/alt/MAPbBr3_S2_model.dat` come out flagged invalid.

**A3 — effective mobility**, Eq (4). Plotting `μ_eff(V)` should give the U-shape of
Fig. 3; its high-V plateau is where μ₀ is read off by hand (or from the ohmic
region, where `p_t ≈ p_f`).

**A4 — concentrations**, Eqs (5) and (6), then `p_s = p_f + p_t` and `Θ = p_f/p_s`.
The workbook computes Θ the other way, as `μ_eff/μ₀` — the two differ (§6.1);
pick one and name it in the output header.

**A5 — Fermi level**, Eq (7) inverted:

```python
E_F = E_v + abs(kT_eV * log(abs(p_f)/N_v))     # workbook MODEL!S
```

Mind which temperature goes where. `MODEL!S` uses the **per-point measured**
`T` (`Data-calculations!N` = `T(°C) + 273.15`) in `kT_eV`, while `N_v` was
computed once at the **scalar** `[device] temperature`. Reproducing
`E_F = -4.933849` at row 9 requires exactly that mix: `kT` at 305.82 K,
`N_v` at 299 K. Since `VJT.dat`'s T column is ambient noise (§6.6), the per-point
choice sprays scatter across `E_F`; make it an option
(`ef_temperature = per_point | scalar`) and default to `per_point` to match.

**A6 — cross-plots.** `p_t = f(p_f)`, and `p_t`, `p_f` against `E_F`. The saturation
of `p_t` is where `N_t` is read off by hand (Fig. S8b,c). No code decision beyond
emitting the columns.

**A7 — p- or n-type.** A human judgement from `E_F0` versus the contact work
functions. The script only needs to emit `E_F0` and the band edges.

**Validity mask.** Emit a `valid` column rather than silently dropping rows. A
point is physical when `J > 0`, `V > 0`, `0 < γ < 1`, and `Θ ≤ 1`. Downstream plots
and any statistic must read `Θ` and `μ_eff` through that mask.

---

## 4. The M-branch (M1–M5): building the model curve

One output row per point of the E_F sweep. Reference: sheets `g(E)`, `FD-funkce`,
`n(E)`, `j(U)`, and `MODEL` columns Y–AR.

### 4.1 Grids (M1)

The workbook uses two aligned uniform grids with the **same 3 meV step**:

| grid | range | n | sheet |
| --- | --- | --- | --- |
| energy `E` | `0` → `-9.0` eV, step `-0.003` | 3000 | `g(E)!A13:A3012` |
| Fermi level `E_F` | `-9.0` → `0` eV, step `+0.003` | 3001 | `n(E)!A9:A3009` |
| kernel `dE = E - E_F` | `+9` → `-9` eV, step `-0.003` | 6001 | `FD-funkce!A9:A6009` |

Build `g(E)` once on the energy grid: `g_tot = g_v + g_c + g_t`, each behind an
on/off flag (`g(E)!E3, H3, L3`) so a run can isolate one contribution.

**On the 3 meV step.** It looks alarming — `k_B T_t = 2.58 meV` at `T_t = 30 K`, so
the trap is sampled at barely one point per width — but the integrals are fine:
refining to 0.2 meV moves J by 0.2 % at both 1 V and 3 V. The real cost is
*sampling of the output curve*: at 3 meV only 9–70 model points (depending on γ)
land inside the measured 0–3 V window, which makes for a visibly chunky J-V.
Make `dE` a `[numerics]` knob, default finer than the workbook (0.5–1 meV is
ample), and state the value used in the output header. Reproducing the §7 targets
exactly still requires the workbook's 3 meV.

### 4.2 The convolution (M2)

`f(E - E_F)` depends only on the difference, so tabulate it *once* on the ΔE grid
and slide a window along it — that is exactly what shifting `'FD-funkce'!B9:B3008`
by one row per E_F step does, and it is why the two grids must share a step:

```python
# f_kernel[j] = f(dE[j]), dE running +9 down to -9
# for the k-th Fermi level (k = 0 at E_F = -9) the window is f_kernel[k : k+n_E]
p_f[k]   = dot(dE_w, g_v   * (1 - f_kernel)[k:k+n_E])            # n(E)!L
p_s[k]   = dot(dE_w, g_tot * (1 - f_kernel)[k:k+n_E]) - p_s0     # n(E)!J
p_t[k]   = p_s[k] - p_f[k]        # n(E)!O subtracts p_f0 instead -- see §6.5
Theta[k] = p_f[k] / p_s[k]        # n(E)!N takes abs()           -- see §6.1
```

`scipy.signal.fftconvolve` is equivalent, but the sliding dot product is only
3001×3000 and is far easier to check line-by-line against the sheet.

The `p_s0` subtraction carries the whole "injected charge" idea and is **not
optional**: `p_s` counts every state below `E_F` as hole-occupied, including the
entire empty conduction band, so the raw sum is `5.0654e27 m⁻³` — dominated by a
large constant that has nothing to do with injection. Subtracting the value at
`E_F0` (`n(E)!K3`, gated by the flag `n(E)!K4`) removes it and leaves
`p_s(E_F0) = 0`. Skip this and the model is meaningless, not merely offset.

`p_f` needs no such treatment: it sums the valence band only, and comes to
`p_f0 = 1.3867e12 m⁻³` at equilibrium.

This is why `examples/MAPbBr3_S2_model.dat` carries both a `p_t` column and a
separate `ptm` column that starts at exactly 0 at V = 0.

The workbook has a *second* subtraction here — `n(E)!O` takes `p_s - p_f0` rather
than `p_s - p_f` — which is a slip; see §6.5.

### 4.3 Voltage and current (M3, M4)

Invert Eq (6) for V, then use (S14) for J:

```python
V      = p_s * e * L**2 / (eps0*eps_r*(1-gamma)*(2-gamma))   # j(U)!C
J      = e * mu_0 * (2-gamma) * V * p_f / L                  # j(U)!E
mu_eff = Theta * mu_0                                        # j(U)!K
```

`V` is driven by the **injected total** space charge `p_s` — that is Poisson's
equation, and it is what `j(U)!C` effectively computes too, since the `p_f0` it
subtracts (1.39e12) is negligible against `p_s` (~1e16). `J` then takes the
**absolute** free density `p_f`, because Ohm's law wants the carriers actually
present, not the increment.

`γ` here is the single scalar of §2, not a per-point quantity — that is the whole
point of M4 ("this parameter does not fundamentally affect the model"). That
claim is false for this parameter set: γ moves J at 3 V by ~1500× (§6.3).

**M5** is the human step: adjust `μ₀` until the modelled `μ_eff` matches the
A-branch plateau. The script only re-runs with a new `μ₀`.

### 4.4 Sweep direction

Injection moves `E_F` **down** from `E_F0` toward `E_v` (holes accumulate). That is
the branch to emit for MAPbBr₃: `examples/MAPbBr3_S2_model.dat` starts at
`E_F = E_F0 = -4.84` with `V = 0` and walks down. The extraction model
(Supplementary Note 4, needed for MAPbI₃ in the dark) is the opposite branch — so
keep the sweep direction a parameter rather than a hardcoded sign.

---

## 5. Suggested layout

```
asclc/
  constants.py    # the workbook constants, nothing else
  params.py       # parse the .params file -> frozen dataclass; validate units
  measurement.py  # A1: read VJT.dat, bin -> V, J, T arrays   (writes *_CDVF.dat)
  analysis.py     # A2-A5: gamma, mu_eff, p_f, p_t, Theta, E_F, valid
  dos.py          # M1: energy grid, g_v, g_c, g_t (3 profiles), N_c, N_v
  model.py        # M2-M4: FD kernel, sliding convolution, V, J, mu_eff
  io.py           # the commented tab-separated .dat writer
  cli.py          # asclc <params> -> CDVF + model + analysis outputs
tests/
  test_against_workbook.py
```

Input format: keep the existing `examples/MAPbBr3_S2.params` INI-style layout
(`[material] [device] [params] [model] [numerics]`, `#` comments, units named in
the comments). It is readable, diffable, and documents itself. Output: keep the
existing commented TSV — a `#`-prefixed provenance block naming the params file
and the model options actually used, then a header row, then data.

---

## 6. Decision points — adjudicated

Nine places where the article, the SI, and the workbook disagree. Each is
labelled by what the evidence says:

* **Artifact** — a slip in the workbook. Do not inherit it; reproduce it only
  behind a `workbook_compat` flag so the §7 regression targets still pass.
* **Load-bearing** — the workbook departs from the SI *and is right to*. Keep it.
* **Open** — genuinely undecided; pick, record the choice in the output header,
  and move on.

Make each a named option in `[model]` and print the active set into the output
provenance block. The recommended default is given for each.

### 6.0 The evidence that settles most of it

Model J at two voltages, against the measured curve, for the shipped parameter
set (`μ₀ = 2.7e-3`, `N_t = 4.7e16`, `E_t = -4.82`, `T_t = 30`, `E_F0 = -4.84`):

| variant | J at 1 V | J at 3 V | vs. measured at 3 V |
| --- | --- | --- | --- |
| **measured** (`examples/..._CDVF.dat`) | ~5.1e-06 | ~1.18e-02 | — |
| γ = T_t/T = 0.1003 (workbook) | 4.2e-06 | 2.14e-02 | **1.8× high** |
| γ = 0.5 (SI M4) | 2.1e-06 | 1.43e-05 | **~830× low** |
| trap `sech` (workbook) vs `sech²` (S5) | ±9 % | ±1.7× | small |
| grid 3 meV vs 0.2 meV | ±0.2 % | ±0.2 % | negligible |

`examples/MAPbBr3_S2_model.dat` gives J(3 V) = 1.45e-05, i.e. it was generated
with γ = 0.5 and undershoots the data it is meant to model by ~800×. That is the
single most consequential choice in the whole port, and §6.3 is where it lives.

### 6.1 Eq (6) is labelled `p_t` but is `p_s` — *open, low stakes*

Dividing Eq (5) by Eq (4) gives `p_f/Θ = eps0*eps_r*(1-γ)(2-γ)V/(eL²)`, which is
*exactly* Eq (6). So the quantity Eq (6) computes is the **total** space charge
`p_s`, and the trapped part is `p_s - p_f`. But the SI (A4, S12) says
`Θ = p_f/(p_f + p_t)` with `p_t` from Eq (6), and the workbook agrees with the SI:
`MODEL!O = K + M`, i.e. `p_s = p_f + Eq(6)`.

Consequence: the workbook carries two mutually inconsistent Θ, and uses the second:

* `Theta_ratio = p_f/(p_f + p_t)` (S12) — computed nowhere in the sheet
* `Theta_mob   = mu_eff/mu_0` (`MODEL!Q`) — what actually feeds `E_F` and the plots

Reading Eq (6) as `p_s` makes Eqs (4), (5), (6) and (S12) one self-consistent
system with a single Θ, and makes the §7 round-trip test exact. Tested on the
model branch, the three readings give J and V agreeing to better than 0.2 %;
they differ only in the floor of Θ (4.7e-05 / 1.5e-04 / 4.6e-04 at the bottom of
the sweep). So this changes the Θ and `p_t` figures (Figs. 4, S19) and essentially
nothing else.

**Default `theta_model = free_over_total`**, i.e. read Eq (6) as `p_s` and set
`p_t = p_s - p_f`, `Θ = p_f/p_s`. Keep `mobility_ratio` (workbook) and
`absolute_over_total` (`|p_f|/(|p_f|+|p_t|)`, which keeps noisy negative-J points
finite in the A-branch) as alternatives.

### 6.2 The workbook's trap DOS is not Eq (S5) — *artifact, low stakes*

`g(E)!M13` reads `N' * e^u/(1 + e^u^2)` — in Excel precedence `^` binds before `+`,
so it is `e^u/(1 + e^(2u))`, which is `½·sech(u)`, not the `e^u/(1+e^u)²` of
Eq (S5). Two independent oddities point the same way: that missing parenthesis,
and a prefactor of `N_t/(2*kT_t)` where the normalised form wants `N_t/kT_t`.
Together they integrate to `(π/4)*N_t ≈ 0.785 N_t` instead of `N_t` — verified
numerically (0.785398 vs 1.000000).

It costs little either way: switching to Eq (S5) moves J by 9 % at 1 V and 1.7× at
3 V, both far inside the γ effect of §6.3.

**Default `trap_profile = si_biexponential`** (Eq S5, normalised, published).
Keep `prototype_sech` for exact workbook reproduction and `gaussian` (Eq S6,
`σ = 2 k_B T_t`) because the paper offers it.

### 6.3 γ in the model branch — *load-bearing; the workbook is right*

SI M4 says `γ = T/(T+T_t)`, clamped to `0.5` when `T_t < T`. The workbook
hardcodes `γ = T_t/T` (`'j(U)'!C6`). Live case: **0.5 vs 0.1003**.

The theory argument favours the SI: for an exponential trap distribution the
characteristic exponent is `l = T_t/T`, the slope is `m = l+1`, and
`γ = 1/m = T/(T+T_t)` exactly. On that reading the workbook mistook the exponent
`l` for `γ` itself.

**The data says otherwise.** With the shipped parameters, γ = 0.1003 lands within
1.8× of the measured current at 3 V while γ = 0.5 is ~830× low (§6.0). The
elegant reading is the one that fails.

The resolution is that **γ and the hand-selected parameters are not separable**.
`N_t`, `E_t`, `T_t` and `μ₀` were tuned by hand *against* the workbook's γ, so the
parameter file and the γ convention travel together. Swapping γ without re-tuning
the other four produces the ~800× miss that `examples/MAPbBr3_S2_model.dat`
exhibits — which is very likely why the previous port, whose commit message reads
"Follow the published equations", was scrapped.

**Default `gamma_model = workbook` (`T_t/T`)**, and refuse to run `SI-M4` against
a parameter file tagged as the workbook set without at least a loud warning.
Note also that `T_t/T` exceeds 1 whenever `T_t > T`, which flips the sign of
`(1-γ)` and yields negative voltages — so it is unusable in exactly the regime
SI M4 was written for. Validate the resulting γ is in `(0,1)` and fail loudly if not.

A third option is worth recording but not yet recommending: take γ per-point from
the model's own local slope, `γ_i = 1/(dlnJ/dlnV)`, making the model internally
consistent with the A-branch definition. A quick fixed-point attempt reached
J(3 V) = 1.08e-02 against a measured 1.18e-02 — the best of the three — but it
**did not converge** (oscillating at ~0.5 amplitude after 60 damped iterations),
so it needs a proper solver before anyone trusts it.

### 6.4 The saved analysis branch ran at γ ≡ 0 — *artifact*

`Data-calculations!I1 = 0` collapses every `LINEST` window to one cell, so the
workbook's saved `γ` and `Ea` columns are 0 at every point, and its whole
analysis branch was evaluated with `(1-γ)(2-γ)² = 4` and `(2-γ) = 2` (§3, A2).
Reproduce only to hit the §7 targets. Real work needs a window of ≥ 2 points.

### 6.5 `p_t = p_s - p_f0` in the model branch — *artifact*

`n(E)!O9 = J9 - $N$3` subtracts the **equilibrium** free density `p_f0` from the
injected total, where (S11) subtracts the free density *at the same `E_F`*. With
`p_f0 = 1.39e12` against a `p_f` reaching ~7e14 the numerical effect is small, but
it is plainly a slip. **Use `p_t(E_F) = p_s(E_F) - p_f(E_F)`.**

### 6.6 Which temperature enters `E_F` — *artifact, but it is what reproduces §7*

`MODEL!S` uses the per-point measured `T` in `kT` while `N_v` was built once at the
scalar `[device] temperature` — two different temperatures in one equation. Since
`VJT.dat`'s `T(C)` column is ambient noise (§6.9), the per-point choice sprays
scatter across `E_F` for no physical reason.

**Default `ef_temperature = scalar`.** Reproducing the §7 A-branch target
(`E_F = -4.933849`) requires `per_point`, so keep that path alive for the test.

### 6.7 Published parameters ≠ workbook parameters — *not a bug*

The workbook's live case is simply a different working state from Table 1 / Table S4:

| | Table 1 (MAPbBr₃, dark) | workbook `MODEL!B` |
| --- | --- | --- |
| μ₀ | 50 cm²/Vs = 5.0e-3 m²/Vs | 2.7e-3 |
| N_t | 3.03e10 cm⁻³ = 3.03e16 m⁻³ | 4.7e16 |
| E_t | -4.768 eV | -4.82 |
| T_t | 8 K | 30 |
| E_F0 | -4.828 eV | -4.84 |
| T | 300 K | 299 |

Ship both as parameter files, each tagged with the γ convention it was tuned
against (§6.3). Reproducing the article's figures needs Table 1; reproducing
`examples/` and the workbook needs the workbook column.

### 6.8 `N_c,v` in Table S4 is off by 10² — *typo, settled*

Table S4 lists `4.52e16 cm⁻³`. Eq (S3) with `m_e* = 0.32`, `T = 299 K` gives
`4.5164e24 m⁻³ = 4.5164e18 cm⁻³` — the workbook's value, mantissa matching
exactly. Read the table as `4.52e18 cm⁻³`. **Always compute `N_c`, `N_v` from
Eqs (S3)/(S4); never read them from the table.**

### 6.9 Minor

* `E_g` from the band edges is 2.22 eV; the article says 2.20 eV. Derive it.
* `MODEL!S` is labelled `EF - Ec` but computes `E_F` on the vacuum scale (it adds
  `B14 = E_v`). Trust the formula, not the label.
* `Data-calculations!O` ("Ea") regresses `ln j` against `e/k_BT` — meaningful only
  for a temperature-modulated sweep, and 0 as shipped (§6.4). Skip it until a
  real T-sweep is loaded.
* The `T(C)` column of `VJT.dat` is ambient noise, not a controlled variable:
  consecutive points read 32.7, 21.3, 11.0, 13.6 °C. Keep it plumbed through,
  keep it optional, and do not read physics into its scatter.

## 7. Validation

Test against the workbook's cached values — they are stored in the file and
readable with `openpyxl(data_only=True)`, no Excel needed. Targets for the live
case (MAPbBr₃ S2 dark; `μ₀ = 2.7e-3`, `N_t = 4.7e16`, `E_t = -4.82`, `T_t = 30`,
`E_F0 = -4.84`, `T = 299`):

```
g(E)!F7   N_v    4.2025506442440113e+24 m^-3
g(E)!I7   N_c    4.5163559637947533e+24 m^-3
n(E)!N3   p_f0   1386704721143.2336 m^-3
n(E)!K3   p_s0   5.065430627403648e+27 m^-3
j(U)!C6   gamma  0.10033444816053512        (workbook convention)

E_F = -4.899 (n(E) row 1376):  p_s 1.001435e+16  p_f 1.370292e+13  Theta 1.368329e-03
                    j(U):      V   1.496545e+00  J   2.808376e-05  mu_eff 3.694487e-06
E_F = -5.001 (n(E) row 1342):  p_s 1.234312e+16  p_f 7.189317e+14  Theta 5.824555e-02
                    j(U):      V   1.844604e+00  J   1.816113e-03  mu_eff 1.572630e-04

A-branch, first measured point (MODEL row 9, V = 7.0311e-3 V, gamma = 0, T = 305.8211 K):
  J 9.492278e-07  mu_eff 4.592328e-03  p_f 9.363531e+13  p_t 5.505167e+13
  p_s 1.486870e+14  Theta 1.700862e+00 (>1 -> invalid)  E_F -4.933849e+00
```

The A-branch row above only reproduces with `γ = 0` (§3, A2) and with `kT` taken
at the point's own 305.8211 K against an `N_v` built at 299 K (§3, A5). All six
numbers were re-derived from Eqs (4)–(7) under exactly those two conventions and
match the workbook to every digit shown; if yours drift, one of the two is wrong.

Row *r* of `n(E)`/`j(U)` corresponds to `E_F = -9 + (r-9)*0.003`; note the grid
lands on `-4.839`, not exactly `-4.84`, so compare at grid points rather than at
nominal energies.

Also assert the invariants — they catch sign and unit errors fast:

* `Θ ≤ 1` wherever a point is flagged valid, and `Θ → 1` in the ohmic limit
* `μ_eff → μ₀` as Mott–Gurney is approached
* `∫ g_t dE = N_t` for `si_biexponential` and `gaussian`; `(π/4)·N_t` for `prototype_sech`
* Eq (2) at `γ = ½` reproduces Eq (3) to machine precision
* round-trip: feeding the M-branch's `(V, J)` back into the A-branch returns the
  `p_f`, `p_t`, `Θ` it started from. This is the strongest single test, and it
  holds only if §6.1 is resolved consistently in both branches.

And one test the workbook cannot provide, because it is the thing the workbook
was built to achieve — **the model must track the data it models**:

```
J_model(1 V) within ~2x of 5.1e-06 A/m2
J_model(3 V) within ~2x of 1.18e-02 A/m2      (measured, examples/..._CDVF.dat)
```

Assert it as a regression guard on the shipped parameter set. It is the only
check that would have caught the ~800× miss in `examples/MAPbBr3_S2_model.dat`,
and none of the workbook-fixture tests above go anywhere near it.

---

## 8. Build order

0. Settle §6.3 (γ) first, out loud, and write the choice into the parameter file
   next to the five variables it is coupled to. Everything downstream inherits it.
1. `constants.py`, `params.py`, and the `.params` round-trip. No physics yet.
2. `measurement.py` → regenerate `examples/MAPbBr3_S2_CDVF.dat` byte-identically.
3. `dos.py` → assert `N_v`, `N_c` against §7 to all digits. Cheapest real check.
4. `model.py` → assert the four `n(E)`/`j(U)` rows of §7 in `workbook_compat`
   mode, then assert the data-agreement guard in default mode. This is the core.
5. `analysis.py` → assert `MODEL` row 9, then the round-trip test.
6. `cli.py`, plotting, and the second parameter set (Table 1).

Get step 3 exact before writing step 4: a wrong `N_v` silently shifts every `E_F`
by a constant, and everything downstream still *looks* plausible.

Two modes are worth carrying permanently, and they are cheap if built in from the
start: `workbook_compat` (every §6 option set to the workbook's behaviour, which
makes the §7 fixtures pass) and the default (the adjudicated settings). Any run
should name its mode in the output header, so a `.dat` file on disk always says
which conventions produced it.
