# Current density versus voltage, linear/log

Implemented at the user's request on 2026-09-15 as
`plot_hole_current_semilog`, displayed in the notebook and registered as
`hole_current_semilog` for the existing export step.

The source chart is `xl/charts/chart3.xml` in `data/SCLCKopecky.xlsx`.
Its measured series uses MODEL E/G, labelled j (A/m²). Its model series
uses MODEL Z/AC, labelled jm (pf); Z9 references j(U)!C9 and AC9
references j(U)!E9, the hole voltage and current. The source's unequal
model range lengths are not propagated: calculated U_p and J_p retain
their paired energy rows. No cached source values supply plotted data.

The calculation and its physical assumptions, units, and unresolved density
convention are unchanged; see [M4 derivation](M4_model_current.md).
This is an additional view of the existing outputs, with measured markers
and a model hole curve, linear voltage, logarithmic current, and automatic
limits. Finite signed voltages are retained. Nonpositive/nonfinite currents
cannot be drawn on the logarithmic axis; invalid model rows leave gaps.
All physical and numerical parameters remain unchanged.

Validation: executed all notebook calculation and plotting cells, excluding
the export cell. Checked both plotted series against their calculated input
arrays, exactly two series, linear/log scales, and automatic limits on both
axes. The full model voltage span dominates the automatic horizontal range,
so measurements occupy a narrow region near zero on this linear-axis view.

User-specified display limits, 2026-09-15: the notebook now sets voltage
to 0–4 V and current density to 10^-8–10^0 A/m², overriding the initial
automatic limits. The plotting function and calculated arrays are unchanged.

Accepted by the user as complete on 2026-09-15 with these display limits.
