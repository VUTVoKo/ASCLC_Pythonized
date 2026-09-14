# Free-hole model population: physics check

Checked with the current notebook inputs on 2026-09-14. No production code,
physical parameters, processing settings or numerical grids changed.
Numerical results are in `pfm_physics_check.json`.

User decision: the current numerical accuracy is acceptable for now.
Defer improved numerical integration to later work; retain the existing
method and grids. This acceptance concerns the integration accuracy, not
the separately identified chart-mapping issues.

## Derivation independent of the source chart

For an isotropic, spin-degenerate parabolic valence band,
E = Ev - hbar^2 k^2/(2 m_h). Counting two spin states in k-space gives
N(k)/volume = k^3/(3 pi^2). Differentiating with respect to energy gives

    g_v(E) = 4 pi (2 m_h e)^(3/2) / h^3 * sqrt(Ev-E), E < Ev.

Here numerical energies are in eV, e converts eV to joules, and m_h is
the effective mass in kg. The result has units m^-3 eV^-1. A hole is an
unoccupied electron state, so

    pf(EF) = integral[-infinity, Ev] g_v(E) [1-f(E,EF)] dE
    1-f(E,EF) = 1 / (1+exp((EF-E)/kT)).

This is the free-hole expression S10 in the supplied supplementary text,
`42005_2025_2202_MOESM2_ESM.txt`, page 4. Its extracted limits and occupation
factor are unambiguous for this check. It uses valence-band DOS only, not
trap DOS, and is an absolute population rather than equilibrium-subtracted.
The implementation's stable logistic evaluation has the correct sign.

With t=(Ev-E)/kT and eta=(Ev-EF)/kT, independent continuum integration gives

    pf = Nv * (2/sqrt(pi)) * integral[0,infinity] sqrt(t)/(1+exp(t-eta)) dt
    Nv = 2 (2 pi m_h k_B T/h^2)^(3/2).

For eta much less than zero, the integral reduces to the Boltzmann result
pf = Nv exp((Ev-EF)/kT). At EF=Ev it gives pf=0.765147... Nv.
Increasing EF reduces hole occupation at every energy, so pf must decrease.
The model's full 3001-row pf array is nonnegative and monotone decreasing.

## Numerical checks

An independent adaptive integral of the dimensionless expression above was
used as a diagnostic reference; it does not resample or replace the model
grid. Existing values: T=299 K, m_h=0.305 electron masses, Ev=-5.58 eV,
EF0=-4.84 eV; kT=0.02576582645 eV and Nv=4.205781137e24 m^-3.

| EF (eV) | Current pf (m^-3) | Continuum reference (m^-3) | Relative difference |
|---|---:|---:|---:|
| -9.000 | 4.798822504e27 | 4.838535410e27 | -0.82076% |
| -5.580 (band edge) | 3.198582754e24 | 3.218040923e24 | -0.60466% |
| -5.187 | 9.900327869e17 | 9.992095498e17 | -0.91840% |
| -4.840 (equilibrium) | 1.402200272e12 | 1.415197479e12 | -0.91840% |
| -4.713 | 1.014304790e10 | 1.023706536e10 | -0.91840% |

Across all existing sweep rows with pf between 1e10 and 1e18 m^-3,
the discrepancy is approximately -0.918402%. The continuum and Boltzmann
references agree at equilibrium to negligible relative error. The bias is
consistent with the retained 0.003 eV rectangular DOS integration, including
sampling the square-root band edge. It is not an occupation-sign error or
a missing spin factor. No grid changes were made to investigate it.

At EF=-9 eV the finite integration boundary matters as well. Integrating
the continuum only over the existing -9..Ev domain gives 4.800468540e27,
versus 4.838535410e27 for the ideal infinite parabolic band. The current
rectangle sum gives 4.798822504e27. This separates the finite-domain effect
from the quadrature approximation without changing the model domain.

## Interpretation

The governing free-hole formula and implemented units/signs are correct
within a single parabolic-band effective-mass model. Numerical evaluation
has a measurable roughly 0.92% bias in the chart window. This is distinct
from the earlier comparison against rounded source inputs.

The huge pf near EF=-9 eV follows from placing EF 3.42 eV below Ev, deep
inside the extrapolated valence band. It is a mathematical result of that
model; accuracy of a parabolic effective-mass approximation over such a
large energy span is not established by this check.

These checks do not validate pt or its population convention. Nor do they
justify the extra secondary-axis pfm curve: pf plotted against itself is
an identity reference, and pf plotted against nf with an energy label is
a separate chart-mapping issue.
