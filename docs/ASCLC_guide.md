# A-SCLC: loading and displaying measured current density

## Scope

This repository currently does one thing: read a standard measurement CSV and
plot the measured current density against voltage. The DOS/occupation model,
the A3 effective-mobility and A4 carrier-density analysis, and the spatial
transport prediction have been removed. Only measured-data loading and the
measured J–V plot remain.

## Parameters are user-controlled

Never change the sample, device, or numerical parameters supplied in the
notebook or in `examples/`. Only an explicit user instruction changing a
specific parameter overrides this. See [AGENTS.md](../AGENTS.md).

## Sources

| Source | Role |
| --- | --- |
| `data/s42005025022021.pdf` | Gavranovic, Zmeskal, Weiter & Pospisil, *Communications Physics* **8**, 280 (2025) |
| `data/42005_2025_2202_MOESM2_ESM.pdf` | Supplementary Information |
| `examples/` | Earlier exported data and calculations, kept for reference |

The XLSX workbook is the user's independent validation tool. Code here does not
read it; this guide refers to it only to record where a parameter came from
(see *Local-slope smoothing window* below).

## Frontend and backend

`notebooks/asclc.ipynb` is the **frontend**: an editable input path and
`device` dictionary, backend calls, and the plot. It contains no file parsing
or unit conversion.

The **backend** modules take paths, arrays, and parameters explicitly and
return data for the frontend to display:

- `asclc_backend.py`: `load_measurements` (standard CSV), `current_density`,
  `trailing_mean` (workbook columns E/G/N), `centered_mean` (an aligned variant
  for plotting), and `local_gamma` (column K).
- `asclc_workflow.py`: `run_calculation` (load + `V`, `J`) and
  `export_results` (figures as PNG/PDF/SVG plus `measured.csv`).
- `asclc_plotting.py`: `plot_measured_jv` and `plot_smoothed_jv(result,
  window=...)` (a centred window-point mean of `J` over the raw markers, at the
  same voltages), each returning `(figure, axes)` in plain publication styling
  with the notation of the published equations.

## CSV format

Standard CSV with the exact header `U(V),I(A)` and two numeric columns:
volts, amperes. A UTF-8 BOM is accepted. Incorrect headers or column
counts are rejected. Signs, row order, and nonfinite measurements are retained.

`run_calculation(path, device=device)` returns `V = U + device['voltage_offset']`
and `J = I / device['area']` (A/m²). The `device` dict is not mutated.

## Constants and units

Voltage in V; current density in A/m². No physical constants are used by the
current code path.

## Notebook sequence

1. Imports and the editable input path and `device` dictionary.
2. `run_calculation(...)` to load the measurements.
3. `plot_measured_jv(result)` — measured J–V on logarithmic axes; nonpositive
   readings stay in the data but cannot appear on a logarithmic axis.
4. `plot_smoothed_jv(result, window=numerics["window"])` — the same data with a
   centred window-point mean of `J` overlaid at the same voltages.
5. `export_results(result, figures, directory)` to write the figures and
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
  column K: `gamma = d ln|U| / d ln|j|` from a least-squares fit over `window`
  points. `centered=False` is the sheet's trailing anchor; `centered=True`
  aligns the result with `voltage`. Not called by the notebook.

Both return an array the length of the input, with `nan` at the ends without a
full window and wherever a window covers a zero or nonfinite reading (the sheet
instead shrinks the window there and, for `local_gamma`, returns `0`).

For the MAPbBr3 S2 dark scan `local_gamma` on the raw current is ≈ 0 at every
voltage — `ln|j|` is dominated by noise, so `d ln|U| / d ln|j|` collapses. That
is the faithful column-K result for this dataset, not a bug: the dark J–V has no
resolved ohmic (`gamma = 1`) region.

`asclc_backend.centered_mean(values, window=...)` is a non-workbook helper: the
same average centred on `i` instead of trailing, `nan` for the points at each
end where the window runs off the data. `plot_smoothed_jv` uses it so the
smoothed line stays aligned with the raw markers on the voltage axis.

## Validation

`test_asclc_csv.py` checks standard CSV loading: unit/column enforcement and the
preservation of signs, order, and nonfinite rows. `test_asclc_smoothing.py`
checks `trailing_mean` and `centered_mean` against explicit windows, alignment,
`window = 1` identity, and `nan` handling. `test_asclc_gamma.py` checks `local_gamma` against analytic
power laws (`j proportional to V**k` gives `gamma = 1/k`), current scale/sign
invariance, the centred variant, and the `nan` handling. Checks against the
reference spreadsheet are the user's, run outside this repository.
