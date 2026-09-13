# M3: Fermi level shift and source mapping

## Implementation decision

On 2026-09-14 the user explicitly requested implementation of the workbook
way. M3 therefore carries the M2 energy sweep directly into
`asclc_workflow.run_model_fermi_level`; it does not invert M2 densities using
the experimental MODEL!S formula. `ModelFermiLevel` exposes absolute `E_F`,
`E_F0`, and signed `delta_E_v`, `delta_E_c`, and `delta_E_F0` in eV.
The notebook runs this after M2 and plots `delta_E_v` against M2 `p_f`.
The figure participates in the existing export step. M2 already supplies
the temperature-dependent relationship between energy and population, so M3
does not take a second temperature or introduce a new energy grid.

No model, material, fitting, or numerical settings are changed. In particular,
the experimental ABS inversion in MODEL!S is separate from this model step.

Validation: executed the notebook calculation and plotting cells with the
project environment (leaving the existing export files alone). All 3001
M3 energy rows and their offsets match cached MODEL!AO and the corresponding
band/reference cells within 2e-12 eV. Checked the band-gap identity
`delta_E_v - delta_E_c = E_c - E_v`, unchanged input coordinates, and
independent output storage. All 69 existing unit tests pass. This verifies
the coordinate mapping; it does not resolve the M2 occupation-definition
issues documented in M2_excel_conventions.md.

## Source evidence

The user's request on 2026-09-14 refers to the article's model step M3,
not a standalone calculation using the example parameter file. Interpret this
step from the supplied article and identify its counterpart in the user's
workbook. No parameter change is implied.

Sources checked directly: Supplementary Information
(`data/42005_2025_2202_MOESM2_ESM.pdf`), p. 7, M3; main article
(`data/s42005025022021.pdf`), p. 5, equation (7); and
`data/SCLCKopecky.xlsx`. Both PDF pages were inspected visually.

Equation (7) defines the hole energy separation as dEF = EF - Ev:

    pf = Nv exp(-dEF / kBT)
    dEF = kBT ln(Nv / pf)

This follows from integrating the valence-band DOS with nondegenerate hole
occupation. The ratio is dimensionless; kBT must be expressed in eV when
the band energies are in eV. It is distinct from EF - EF0.

The model already uses EF as its independent occupation-sweep coordinate.
`n(E)!A` contains EF. `n(E)!B9 = A9-'g(E)'!$B$8` subtracts Ec;
`MODEL!AO9 = ('n(E)'!B9)+$B$13` adds Ec back, returning absolute EF.
Thus the model energy separation for holes is MODEL!AO - MODEL!B14.
This mapping is an interpretation of M3 from the workbook's data flow;
M3 itself gives no additional equation. Inverting equation (7) on the M2
free-hole integrals is a nondegenerate consistency check, not an independent
new model energy coordinate over the entire sweep.

The experimental inversion is in MODEL!S. For example:

    S9 = ABS($S$7*$F$2*'Data-calculations'!N13/$F$1*LN(ABS(K9/$B$18)))+$B$14

Here K9 is the inferred free-hole density, B18 is Nv from g(E)!F7, B14 is
Ev, F2 is kB in J/K, F1 is elementary charge, and S7 is 1. Despite its
header, S contains absolute EF. Its ABS operations agree with the signed
Boltzmann inversion only for positive pf <= Nv. The model AO column also
has a relative-energy header despite containing absolute EF.

Nv uses the sheet temperature Data-calculations!B6 through g(E), whereas
S uses per-row temperatures from Data-calculations!N. See ASCLC_guide.md
for the resulting temperature-convention difference. Preserve existing
temperatures and other parameters.
