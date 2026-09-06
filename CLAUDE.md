# Working on this repository

A Python port of the A-SCLC model from Gavranovic et al., *Commun. Phys.* **8**, 280
(2025), reimplemented from a reference Excel workbook. `docs/ASCLC_spec.md` is the
specification and is authoritative; keep it in step with the code.

## Three things are unresolved in the sources, not in the code

Do not silently pick one. They are switchable options with documented defaults, and the
default is whichever reproduces the reference workbook so existing results stay
comparable. Section 7 of the spec has the evidence for each.

1. **Trap density of states.** The workbook computes `sech(u)/2` where Eq (S5) specifies
   `sech²(u/2)/4`, and its prefactor normalizes neither form — the workbook's trap
   integrates to `(π/4)·N_t` rather than `N_t`. `TrapProfile.SI_BIEXPONENTIAL` is the
   better choice for new work; `WORKBOOK_SECH` is the default for validation.
2. **γ in the model branch.** `T_t/T` (workbook) or `T_t/(T+T_t)` (the SI). This is *not
   identifiable from J-V data*: γ rescales V and J by constants, so it translates the
   curve on log-log without changing its shape. Do not propose fitting it.
3. **Θ.** Eq (S12) supports three readings that agree at high injection and differ by up
   to 6× through the trap-filling region. The port follows the workbook's, which is
   unbounded; `valid` flags where it exceeds 1.

## Conventions that are easy to get wrong

The sources label several quantities as one thing and compute another. This accounts for
most of the subtle errors made while porting:

| Label | What it is |
|---|---|
| `p_t` in Eq (6) | total space charge, not the trapped part |
| `ps`, `ns` in the workbook's `n(E)` | injected totals over the full DOS |
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
against the reference workbook's own computed values and against the article's published
parameter sets; those catch things that checking the code against its own docstrings
cannot. When adding a claim to the spec, add the test that verifies it.

```bash
pytest -q && ruff check .
```

Before changing anything numerical, run the model branch end to end on
`data/MAPbBr3_S2_dark.txt` and confirm the workbook regressions still hold.

## Not implemented

- **Extraction model** (Supplementary Note 4), needed for the MAPbI₃ dark curve where an
  injection barrier makes extraction dominant. The model branch sweeps `E_F` downward
  from `E_F0` only.
- **Parameter fitting.** Both branches are pure functions, so a residual around
  `model_curve` handed to `scipy.optimize.least_squares` is the natural next step. Note
  the γ degeneracy above when choosing a parameterisation.
