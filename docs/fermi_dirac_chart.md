# fFD(E), 1-fFD(E)

Complete and accepted by the user on 2026-09-15, including axis labels
and X display limits of −1 to 1 eV. The notebook registers
`fermi_dirac` for export after the model occupation calculation.

## Source mapping

Inspected `data/SCLCKopecky.xlsx`, `xl/charts/chart14.xml` and
`FD-funkce` formulas. The two series are B and C against A. B9 is
`1/(EXP(A9/$B$4)+1)` and C9 is `1/(EXP(-A9/$B$4)+1)`;
B4 is thermal energy in eV at the measurement temperature. Thus A is
energy relative to the Fermi level, not absolute energy. Both axes are
linear; the source has no axis titles. Legend labels are retained with
mathematical typography.

Source X references extend to row 9009 whereas Y references end at 6009.
The implementation does not reproduce this mismatch or copy cached data
or source fixed limits. It uses the existing `model_energy` grid minus
the existing `E_F0`, without resampling, and the existing temperature and
physical constants. Both display limits are automatic. All pre-existing
parameters and arrays remain unchanged.

Subsequent user override, 2026-09-15: added relative-energy and occupation-
probability axis labels and notebook X display bounds of −1 to 1 eV.
This supersedes the absent source titles and automatic X limits above;
Y limits remain automatic. The energy grid and probabilities are unchanged.

## Derivation

Assume independent fermionic states in thermal equilibrium, each either
empty or occupied. Their grand-canonical weights are 1 and
`exp(-(E-E_F)/(k_B T))`, respectively. Normalizing their sum gives

    f = 1 / (1 + exp(x))
    1-f = 1 / (1 + exp(-x))
    x = (E-E_F) / (k_B T).

The numerator energy and thermal energy are both in eV, with
`k_B T / e` converting the project's SI constants to eV. Probabilities
are dimensionless; no DOS or density normalization enters this chart.
No spatial boundary conditions are needed for a single-state occupation.
The backend evaluates each branch as `exp(-logaddexp(0, ±x))`, avoiding
exponential overflow and subtraction cancellation in small tails.

## Validation

Executed the notebook imports, unchanged parameter and occupation cells,
and the new chart cell. Verified complementarity, probability bounds,
monotonicity, saturated tails, half occupancy at zero relative energy,
and the independently derived odds identity `log((1-f)/f) = x` where
both probabilities are resolved. Verified two lines, linear axes and
automatic limits; inspected the rendered preview.
