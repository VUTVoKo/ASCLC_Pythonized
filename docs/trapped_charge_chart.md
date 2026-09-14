# Trapped charge versus free charge

Implemented 2026-09-14 at the user's explicit request to reproduce the
existing chart logic. Source: `data/SCLCKopecky.xlsx`, MODEL sheet,
`xl/charts/chart5.xml`. Runtime reads no workbook or cached chart data.

## Series mapping

| Source X / Y columns | Calculated X / Y | Axis pair | Legend |
|---|---|---|---|
| MODEL K / M | carriers.p_f / carriers.p_t | primary | pt (m-3) |
| MODEL K / K | carriers.p_f / carriers.p_f | primary | pf (m-3) |
| MODEL AH / AJ | model.p_f / model.p_t | primary | ptm |
| MODEL AH / AH | model.p_f / model.p_f | primary | pfm |
| MODEL AH / AI | model.p_f / model.n_t | primary | ntm |
| MODEL AG / AH | model.n_f / model.p_f | secondary | pfm |

The primary axes are logarithmic and retain the labels Free charge, nf
(m-3) and Trapped charge, nt (m-3). The independent secondary horizontal
axis is linear and reversed, labelled Energy, EF (eV); its independent
logarithmic vertical axis is hidden. All six legend entries are retained,
including the two pfm entries. Limits are derived from measured pt and pf;
model curves cannot widen them. Secondary Y uses the same measured limits;
secondary X uses model rows within the measured free-charge window.
All full model arrays remain present and are clipped only by the axes.
There are no copied range
endpoints, fixed bounds, parameter changes, or source styling.

## Meaning and retained inconsistencies

Concentration comparisons pair values from the same measurement interval
or model energy row. Identity series obey y=x. Every actual coordinate
here has dimensions m-3. There is no additional physical inversion or
conversion in this plotting operation.

AG9 references n(E)!E9 (free-electron occupation integral); AH9 references
n(E)!L9 (free-hole occupation integral). AI/AJ reference the labelled
trapped populations in n(E)!H/O. Therefore the secondary series is pf(nf),
not pf(EF), despite the retained energy label. It is not a physically valid
energy calibration of the bottom axis. Autoscaling exposes concentration-
sized top-axis values, unlike the source's fixed negative energy bounds.
The primary n-labelled axes also retain their source labels for p series.
These inconsistencies are preserved deliberately under the user's latest
instruction, rather than repaired through an invented coordinate mapping.

The existing carrier calculations and their documented density conventions
are reused unchanged. This reproduces the chart's series wiring, not a claim
that every upstream physical convention equals the source formulas.
Nonfinite values and nonpositive logarithmic coordinates leave gaps; signed
input arrays are preserved. The full existing grids are used, without
reproducing unequal or overextended source range lengths.

## Integration and validation

`plot_trapped_charge` returns the figure and two independent axes. The
notebook calls it after model occupations and registers `trapped_charge`
with the existing figure exporter. A rendered preview is saved as
`outputs/trapped_charge.png`.

Executed all notebook calculation/plotting cells with the existing inputs,
skipping the export call. Verified all six plotted X/Y arrays against the
mapping above, six legend entries, log/log and linear/log scales, reversed
top axis and hidden independent secondary Y. After the measured-limit update,
verified primary bounds against a separate measured-only plot and equal
primary/secondary Y limits. Rendered
and visually inspected the preview. No physics implementation changed.
