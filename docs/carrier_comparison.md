# Calculated carrier comparison against saved source results

Current notebook parameters retained. Comparison only; no source data supplied to calculations.

| Series | Finite pairs | Median relative error | Max relative error (row) | Sign disagreements |
|---|---:|---:|---:|---:|
| p_f | 3001 | 0.022900% | 4.459081% (2180) | 0 |
| p_t | 3001 | 0.022877% | 0.716093% (1396) | 0 |
| n_f | 3001 | 0.022900% | 4.460891% (858) | 0 |
| n_t | 3001 | 0.022877% | 0.649317% (1294) | 0 |

### Model errors within the displayed free-charge range

p_f within source bottom-axis free-charge window: 159 rows; median relative error 0.962151%, max 1.297019%.
p_t within source bottom-axis free-charge window: 159 rows; median relative error 0.060313%, max 0.716093%.
n_f within source bottom-axis free-charge window: 159 rows; median relative error 2.325459%, max 2.664849%.
n_t within source bottom-axis free-charge window: 159 rows; median relative error 0.033019%, max 0.649317%.

## Measured series

p_f: 204 valid corresponding rows; median relative error 57.670775%, max 8430.234642%.
p_t: 204 valid corresponding rows; median relative error 10.843821%, max 1248.115537%.

## Visible model identity series

159 source pfm points fall in the primary 1e10..1e18 X/Y window; first/last rows 1280/1438. Both coordinates reference AH, so this is exactly y=x.

## Why measured-derived values differ

The row comparison above is by corresponding original measurement index,
not identical processed voltage. Inspected raw voltage/current-density rows
agree to CSV rounding, but processing differs. Source Data-calculations!I1
is 0, while the current notebook uses numerics["window"] = 11. Source row
K13 has saved gamma=0. The notebook uses windowed means and fitted gamma.
These existing settings were not changed for this comparison.

For MODEL row 109, source U=1.4132591682 V and j=2.7692527e-5 A/m2;
the notebook's corresponding processed interval is U=1.4835705696 V and
j=1.3337049e-5 A/m2. Source pf=1.3590510e13, Python pf=6.2071893e12 m-3.
Thus the measured points are not numerical reproductions of the saved
source-derived values. The current trapped-density convention is also
retained, not replaced during this plotting work.

## Why pfm appears differently

The primary pfm series uses MODEL!AH for both X and Y; 159 saved points lie
inside the source's primary bounds. It is a red dashed identity line. The
measured pf series uses MODEL!K for both coordinates and therefore lies on
the same diagonal irrespective of the numerical differences. Our solid
line makes that relationship more prominent. Without inspecting an actual
Excel rendering, coincidence/style is an explanation, not proof of why
the user cannot distinguish that line visually.

The other pfm series uses AG/AH on independent axes. Its saved X values
are all positive, whereas the source secondary X bounds are -5.4..-4.6.
It is entirely outside the source view. Our automatic secondary limits
made this extra curve visible and allowed apparent crossings with primary
series. Those crossings compare different horizontal coordinates, not equal
carrier states. This is a confirmed display mismatch, independent of the
small model value differences. No plotting code was changed in this audit.

Full 3001-row model comparison: carrier_comparison.csv. Relative errors in
very small populations can be larger than those in the chart window;
near-zero excess populations also magnify relative errors. The current
model uses its existing DOS and constants, so the first-row 0.0226% figure
must not be generalized to every row. All four model series agree in sign
across the saved numerical rows.
