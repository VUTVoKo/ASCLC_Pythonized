# asclc

Python implementation of the advanced space-charge-limited current (A-SCLC) model of

> S. Gavranovic, O. Zmeskal, M. Weiter, J. Pospisil,
> *Advanced space-charge-limited current model for analyzing Fermi level shift in the
> bandgap of halide perovskites*, **Communications Physics 8, 280 (2025)**,
> [doi:10.1038/s42005-025-02202-1](https://doi.org/10.1038/s42005-025-02202-1)

ported from the reference Excel prototype. It extracts mobility, carrier concentrations
and Fermi level position from measured J-V curves, and generates model curves from a
set of material parameters.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,plot]"
```

## Use

Run the model branch on the bundled reference configuration:

```bash
python asclc.py
```

Run a different configuration — `params/` ships the prototype's sample and the
article's four published sets:

```bash
python asclc.py --params params/MAPbBr3_light.toml
```

Analyse a measurement alongside it:

```bash
python asclc.py --data data/MAPbBr3_S2_dark.txt --celsius --window 7
```

Write the figures — model curves over the extracted points, the comparison the
method turns on:

```bash
python asclc.py --data data/MAPbBr3_S2_dark.txt --celsius --plot figures/
```

That produces `jv`, `mobility`, `concentrations`, `bandgap_map`, `fermi_level`,
`theta` and `pt_vs_pf`, matching the prototype's embedded charts and the
article's Figs. 2-6. Points failing the `valid` mask are drawn hollow rather
than dropped.

From Python:

```python
import asclc

material, device, params = asclc.MAPBBR3_S2

# model branch: sweep the Fermi level, get a J-V curve
curve = asclc.model_curve(params, material, device, V_max=3.0)

# analysis branch: read a measurement and extract parameters
m = asclc.load_jv("data/MAPbBr3_S2_dark.txt", device, temperature_in_celsius=True)
res = asclc.analyse_jv(m.V, m.J, material, device,
                       mu_0=params.mu_0, temperature=m.temperature)

# read results through the validity mask
print(res.E_F[res.valid])
```

**Always read results through `.valid`.** Points in the ohmic region, where the measured
current is noise straddling zero, produce slopes that are not slopes and quantities
outside their physical range. On the bundled measurement 114 of 214 points fail.

## Tests

```bash
pytest -q      # 97 tests
ruff check .
```

The suite includes regression tests against the prototype's own computed
values, and runs all four parameter sets published in the article's Table S4.

## Documentation

[`docs/ASCLC_spec.md`](docs/ASCLC_spec.md) is the specification: equations, numerics,
validation results, and the open items. Read section 7 before trusting extracted
parameters — three quantities are ambiguous in the published sources, not in this code.

[`docs/prototype_differences.md`](docs/prototype_differences.md) records where this port
departs from `SCLCKopecky.xlsx`, the Excel prototype the model was first written in. The
article and its Supplementary Information are the reference; the prototype is kept as an
independent numerical oracle for the regression tests. Run
`params/prototype_S2.toml` to reproduce its results.

## Data

`data/MAPbBr3_S2_dark.txt` is the raw measurement from the prototype spreadsheet
(MAPbBr3 single crystal, sample S2, dark, 2024-09-27), as voltage / current /
temperature columns.
