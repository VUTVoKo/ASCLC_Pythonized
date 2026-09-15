# pt (EF), nt (EF), pf (EF), nf (EF)

Display aligned on 2026-09-15; awaiting user review.

Inspected `data/SCLCKopecky.xlsx`, `xl/charts/chart11.xml`. All four
series use `n(E)!A9:A3009`, the absolute Fermi energy in eV. Ordinates
are respectively columns O, H, L and E, mapped to the existing calculated
arrays `p_t`, `n_t`, `p_f` and `n_f` in m⁻³. No energy offset is applied.
The legend uses the four source labels with mathematical typography.
The source has no axis titles or equilibrium marker, so neither is drawn.
X is linear and Y logarithmic.

The plotting function uses automatic limits. Existing notebook display
settings are preserved: reversed X from −2 to −7 eV and Y minimum
10e10 = 10¹¹ m⁻³. No source fixed bounds were copied. No parameters,
grids, occupation definitions or calculated arrays were changed.

## Mathematics and limitations

Each state contributes `g(E) dE f(E, EF)` electrons and
`g(E) dE [1-f(E, EF)]` holes, with Fermi–Dirac occupation `f`.
DOS in m⁻³ eV⁻¹ times energy width in eV gives number density in m⁻³.
Band integrals give absolute free populations. The implemented labelled
trapped populations subtract equilibrium total and equilibrium free
populations from the total occupation integral. They therefore are not
direct trap occupation integrals; see the definitions and limitations in
[M2 conventions](M2_excel_conventions.md). This display change does not
alter or resolve those definitions. Signed values remain in the model;
nonpositive and nonfinite values leave gaps on the logarithmic axis.

## Validation

Executed the notebook's imports, unchanged parameter cell and occupation
cell. Verified four series against the calculated arrays, their energy
coordinates, linear/log scales, absent axis titles and preserved notebook
bounds. Rendered and inspected `/tmp/occupation_chart.png`.
