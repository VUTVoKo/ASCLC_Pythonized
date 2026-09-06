# Working on this repository

A Python port of the A-SCLC model from Gavranovic et al., *Commun. Phys.* **8**, 280
(2025). **The article and its Supplementary Information are the reference.**
`docs/ASCLC_spec.md` is the specification and is authoritative; keep it in step with the
code.

The model was first written as an Excel workbook, `data/SCLCKopecky.xlsx`. That was a
prototype, not a specification: it departs from the published equations in three places
and has several numerical defects. `docs/prototype_differences.md` lists every one. It is
kept because it is the project's only independent numerical oracle — see Verification.

## Three things the sources genuinely disagree on

Do not silently pick one. They are switchable options; **the default follows the
published equation**, and the prototype's reading stays available for reproducing
pre-port results. Section 7 of the spec has the evidence for each.

1. **Trap density of states.** Eq (S5) specifies `sech²(u/2)/4`, which integrates to
   `N_t`; the prototype computes `sech(u)/2`, which integrates to `(π/4)·N_t`, so `N_t`
   read off a saturation plot per step A6 comes out 27 % low.
   `TrapProfile.SI_BIEXPONENTIAL` is the default; `PROTOTYPE_SECH` is retained.
2. **γ in the model branch.** `GammaModel.SI_NOTE_M4` is the default: Supplementary
   Note M4 as printed, `T/(T_t+T)` for `T_t ≥ T` and `0.5` otherwise — and every
   published configuration has `T_t < T`, so **γ = 0.5 throughout**. `TT_OVER_T` =
   `T_t/T` is what the prototype computed. A third option, `TT_OVER_T_PLUS_TT` =
   `T_t/(T+T_t)`, is in **no source**; it was once mislabelled as the SI's and is kept
   only for backward comparison. γ is *not identifiable from J-V data*: it rescales V and
   J by constants, so it translates the curve on log-log without changing its shape. Do
   not propose fitting it.
3. **Θ.** Eq (S12) supports three readings that agree at high injection and differ by up
   to 6× through the trap-filling region; pick with `ThetaModel`. The default is Eq (S12)
   as printed, `ABSOLUTE_OVER_TOTAL` = `p_f/(p_f + p_t)`, which is bounded by 1. The
   prototype's `ABSOLUTE_OVER_INJECTED` = `p_f/p_inj` is unbounded and `valid` flags
   where it exceeds 1. The analysis branch takes no such switch — there Eq (4) fixes
   Θ = `μ_eff/μ₀`.

## Conventions that are easy to get wrong

The sources label several quantities as one thing and compute another. This accounts for
most of the subtle errors made while porting:

| Label | What it is |
|---|---|
| `p_t` in Eq (6) | total space charge, not the trapped part |
| `ps`, `ns` in the prototype's `n(E)` | injected totals over the full DOS |
| `ptm` in Figs 4c,d and 6 | the injected total (`p_space_charge`), not `p_t` |

Also:

- **Eqs (2) and (4) carry `(1−γ)(2−γ)²`**; Eq (6) the first power. The Mott–Gurney limit
  fixes it: at γ = ½ the factor is 9/8.
- **Eq (7) inverts to `ln(N_v/p_f)`**, not `ln(p_f/N_v)`.
- **The injected charge integrates the full DOS.** The conduction band contributes
  exactly zero on the hole side but the valence band dominates the electron side.
  `n_inj = −p_inj` exactly, by charge conservation.
- **Temperature is per-point.** Eq (7) uses the measured temperature of each data point.
  In the bundled measurement the 282–314 K spread moves `E_F` by 0.075 eV, larger than
  the Fermi level shifts the article reports.

## Validity

`AnalysisResult`, `ModelCurve` and `CarrierDensities` all carry a `valid` mask meaning
the same thing. **Every example, plot and derived number must read through it.** Roughly
half the points of a real measurement fail it, and the failures look like ordinary
numbers.

## Verification

Prefer an independent reference over self-consistency. The regression tests compare
against the prototype's own cached cell values and against the article's published
parameter sets; those catch things that checking the code against its own docstrings
cannot. The prototype is not authoritative — the defaults no longer follow it — so each
of those tests selects its options explicitly via `PROTOTYPE_DENSITIES`. Its analysis
branch is degenerate (`Data-calculations!I1 = 0` returns γ = 0 for every row) and is not
a regression target; only its model-branch columns are.

When adding a claim to the spec, add the test that verifies it.

```bash
pytest -q && ruff check .
```

Before changing anything numerical, run the model branch end to end on
`data/MAPbBr3_S2_dark.txt` and confirm the regression tests still hold.

## Not implemented

- **Extraction model** (Supplementary Note 4), needed for the MAPbI₃ dark curve where an
  injection barrier makes extraction dominant. The model branch sweeps `E_F` downward
  from `E_F0` only.
- **Parameter fitting.** Both branches are pure functions, so a residual around
  `model_curve` handed to `scipy.optimize.least_squares` is the natural next step. Note
  the γ degeneracy above when choosing a parameterisation.
