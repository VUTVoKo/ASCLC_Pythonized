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

The XLSX workbook is the user's independent validation tool. No module,
notebook cell, test, or document here reads or refers to it.

## Frontend and backend

`notebooks/asclc.ipynb` is the **frontend**: an editable input path and
`device` dictionary, backend calls, and the plot. It contains no file parsing
or unit conversion.

The **backend** modules take paths, arrays, and parameters explicitly and
return data for the frontend to display:

- `asclc_backend.py`: `load_measurements` (standard CSV) and `current_density`.
- `asclc_workflow.py`: `run_calculation` (load + `V`, `J`, `T_K`) and
  `export_results` (figures as PNG/PDF/SVG plus `measured.csv`).
- `asclc_plotting.py`: `plot_measured_jv`, returning `(figure, axes)` in plain
  publication styling with the notation of the published equations.

## CSV format

Standard CSV with the exact header `U(V),I(A),T(C)` and three numeric columns:
volts, amperes, Celsius. A UTF-8 BOM is accepted. Incorrect headers or column
counts are rejected. Signs, row order, and nonfinite measurements are retained.

`run_calculation(path, device=device)` returns `V = U + device['voltage_offset']`,
`J = I / device['area']` (A/m²), and `T_K = T_C + 273.15`. The `device` dict is
not mutated.

## Constants and units

Voltage in V; current density in A/m²; temperature returned in K. No physical
constants are used by the current code path.

## Notebook sequence

1. Imports and the editable input path and `device` dictionary.
2. `run_calculation(...)` to load the measurements.
3. `plot_measured_jv(result)` — measured J–V on logarithmic axes; nonpositive
   readings stay in the data but cannot appear on a logarithmic axis.
4. `export_results(result, figures, directory)` to write the figure and
   `measured.csv`.

Changing the dataset or the `device` values is a configuration edit, not a code
change. Any measurement row count is accepted.

## Validation

`test_asclc_csv.py` checks standard CSV loading: unit/column enforcement and the
preservation of signs, order, and nonfinite rows. Checks against the reference
spreadsheet are the user's, run outside this repository.
