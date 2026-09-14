# Chart reproduction scope

User clarification, 2026-09-14: reproduce the charts in
`data/SCLCKopecky.xlsx` with respect to labels, plotted content, and each
axis's linear or logarithmic scale. Content includes the actual series,
axis quantities, units, and energy references used by each chart.

Reproduce that content using this project's calculations and measurement
inputs. Do not copy or hard-code plotted data from the workbook, including
cached chart values. The source charts define which quantities and series
to show; they do not supply the implementation's plotted numerical data.

Use automatic axis limits; do not copy fixed bounds. Styling need not match
the source charts.

Inspect the actual chart definitions and series references before choosing
an energy offset; do not infer a chart's reference from a related formula
or column alone. Keep existing physical and numerical parameters unchanged.
Source provenance belongs in docs; code and displayed labels describe the
mathematical quantities directly.

Identify charts in discussion and task lists by their visible axis or legend
labels. Keep chart XML numbers only as secondary references. For unlabelled
charts, identify the sheet and source range explicitly.

For Trapped charge versus Free charge, determine display limits from the
measured pt and pf series only. Model curves must not widen these limits.
Align the hidden secondary Y limits with the primary Y limits and fit the
secondary X range using model rows within the measured free-charge window.
