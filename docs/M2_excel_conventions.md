# M2: Excel conventions and known inconsistencies

## Implementation decision

Follow `data/SCLCKopecky.xlsx` for M2, including its equilibrium subtraction,
free-carrier definitions, labelled trapped-carrier columns, and absolute-value
theta formulas. This decision was confirmed on 2026-09-14. The equations below
specify the workbook behavior to reproduce; they do not establish that its
labels have the usual physical meanings. Do not silently replace them with
absolute trap occupations, excess free populations, or the printed SI formulas.
Keep all existing parameters and numerical settings fixed.

M2 is implemented by `asclc_backend.model_occupations`, connected through
`asclc_workflow.run_model_occupations`, and plotted in the notebook with
`asclc_plotting.plot_model_occupations`. The returned `ModelOccupations`
contains both carrier types and their equilibrium reference values. The
free and trapped populations `p_t`, `n_t`, `p_f`, and `n_f` share a single
chart. The
notebook reuses the existing `E_F0 = -4.84` eV from `MODEL!B22` and the example
parameter file, and uses the workbook's integration and Fermi-sweep grids.

The workflow retains the backend's existing physical constants and M1 DOS;
it does not replace them with Excel's rounded constants. Workbook regression
tests feed cached DOS arrays and the workbook's thermal energy directly into
the occupation integrator, isolating M2 from these pre-existing M1 differences.
Thus the workflow follows Excel's formulas but is not bit-for-bit Excel output.
Zero denominators return IEEE `inf` or `nan` in place of Excel division errors.
The plots leave gaps for nonpositive/nonfinite values; the returned arrays keep
them. Neither the populations nor theta are clipped.

## Workbook calculation

Source: [SCLCKopecky.xlsx](../data/SCLCKopecky.xlsx), sheets `g(E)`,
`FD-funkce`, and `n(E)`. Reference equations: S7–S11 in the
[Supplementary Information](../data/42005_2025_2202_MOESM2_ESM.pdf), pp. 3–4.

Write `G = g_v + g_c + g_t`, with the workbook's component switches included.
Let `f(E, EF) = 1 / (1 + exp((E - EF) / kBT))`. Energies and `kBT` are in eV;
DOS is in m⁻³ eV⁻¹. Define the discrete occupation sums

```
A_n(EF) = sum_i dE_i G(E_i) f(E_i, EF)
A_p(EF) = sum_i dE_i G(E_i) [1 - f(E_i, EF)]
n_f(EF) = sum_i dE_i g_c(E_i) f(E_i, EF)
p_f(EF) = sum_i dE_i g_v(E_i) [1 - f(E_i, EF)]
```

All four quantities have units m⁻³. `SUMPRODUCT` uses `g(E)!B13:B3012`
for positive energy widths and the corresponding DOS rows. There are 3000
rectangles of width 0.003 eV on the descending energy grid from 0 to −9 eV;
the final grid point supplies the last interval boundary. This is not a
trapezoidal integral. Moving the `FD-funkce` range by one row advances `EF`
by 0.003 eV while leaving the DOS range fixed.

The equilibrium values use `EF0 = MODEL!B22` through `g(E)!B9`, with the
occupation functions in `g(E)!R` and `T`:

| Quantity | `n(E)` cell |
| --- | --- |
| `A_n(EF0)` | `C3` (also calculated in `C2`) |
| `A_p(EF0)` | `K3` (also calculated in `K2`) |
| `n_f0 = n_f(EF0)` | `F3` |
| `p_f0 = p_f(EF0)` | `N3` |

With the existing subtraction switches `C4 = K4 = 1`, each sweep row is:

| Workbook label | Actual expression | First row |
| --- | --- | --- |
| `n_s` | `A_n(EF) - A_n(EF0)` | `C9` |
| `n_f` | Absolute `n_f(EF)` | `E9` |
| `Qn` | `ABS(n_f / n_s)` | `G9` |
| `n_t` | `n_s - n_f0` | `H9` |
| `p_s` | `A_p(EF) - A_p(EF0)` | `J9` |
| `p_f` | Absolute `p_f(EF)` | `L9` |
| `Qp` | `ABS(p_f / p_s)` | `N9` |
| `p_t` | `p_s - p_f0` | `O9` |

In particular, `H9 = C9-$F$3` and `O9 = J9-$N$3`. Preserve these fixed
equilibrium references when reproducing Excel.

## Inconsistencies and their consequences

1. **Total and free populations use different reference states.** The `s`
   columns are excess total occupations relative to equilibrium, while the
   `f` columns are absolute band occupations. They are not a consistent
   decomposition into absolute free and trapped populations.

2. **The trapped columns subtract equilibrium free carriers.** Normally
   `p_s = p_f + p_t` would require subtracting the free population at the same
   `EF`. Excel instead gives `p_f + p_t = p_s + p_f - p_f0`, and likewise for
   electrons. Its trapped columns include the change in free occupation and
   cannot be interpreted as direct trap occupation integrals.

3. **Theta is not a bounded occupation fraction.** `ABS(p_f/p_s)` and
   `ABS(n_f/n_s)` combine the reference states above. They can exceed one.
   In exact arithmetic at `EF = EF0`, both excess totals vanish while the
   absolute free populations remain finite, so theta is singular. The
   labelled trap populations there are `-p_f0` and `-n_f0`. A sampled sweep
   need not contain `EF0` exactly. `ABS` removes a ratio's sign but does not
   repair these definitions; do not clip the results into `[0, 1]`.

4. **The calculation differs from the printed SI trap integrals.** S8
   integrates from `EF0` to `Ec`; S11 integrates from `Ev` to `EF0`. Excel
   integrates the total DOS across its full finite grid and subtracts an
   equilibrium occupation. A fixed energy cutoff and an equilibrium
   occupation subtraction are different operations at finite temperature.

5. **Large equilibrium backgrounds create cancellation risk.** Cached
   `A_n(EF0)` and `A_p(EF0)` are approximately `4.83e27` and `5.07e27` m⁻³.
   Their differences can be much smaller, so subtracting separately computed
   sums can lose precision. This is a numerical limitation, distinct from
   the density-definition issues. Algebraically equivalent evaluation and
   agreement with cached Excel values must not be confused with proof of
   physical correctness.

## Trap normalization carried into M2

The existing M1 profile reproduces the workbook's symmetric sech formula.
As derived in the [project guide](ASCLC_guide.md), its full continuous area
is `(pi/4) N_t`, not `N_t`. `g(E)!M7` is labelled m⁻³ but its prefactor has
units m⁻³ eV⁻¹. Retain the formula and parameter values; do not renormalize
the profile to compensate. Its area limits direct trap occupations, but it
does not bound the workbook's labelled `n_t` and `p_t` columns, which contain
other contributions as described above.

## Validation meaning

`tests/data/m2_reference.npz` stores the unchanged cached numerical inputs and
outputs extracted from `SCLCKopecky.xlsx`: the 3001-point DOS grid, total and
band DOS arrays, equilibrium Fermi level, thermal energy, and numeric `n(E)`
cells. Tests load this fixture without depending on spreadsheet libraries or
opening the original workbook. Source attribution stays in this document;
code, notebook text, and plot labels use mathematical terminology.

Check the occupation sums and cell references against the workbook to verify
reproduction. Independently check units, electron/hole occupation
complementarity, and equilibrium subtraction. Report the failed physical
identities above as known properties of the Excel convention rather than
changing formulas to make those identities pass.
