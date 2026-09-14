# g(E), dn/dEF versus Energy, E, delta EF

Status: accepted by the user as complete. Retain the agreed X range
-5.8 to -2.8 eV, both Y ranges 1e11 to 1e29, and displayed DOS band-edge
extensions below 1e12. Further numerical improvements are deferred.

Implemented at the user's request. Source inspection: MODEL chart6.xml,
MODEL!P/O/AL/AQ/AR and n(E)!J/K. No physical parameters, smoothing settings,
energy grids or occupation conventions changed. Runtime calculations and
plots read no source workbook or cached chart values.

## Derivatives and their meaning

Measured MODEL!O=K+M is pf+pt. P9=ABS((O10-O9)/(S9-S10)) yields

    D_i = abs((p_sum[i+1]-p_sum[i])/(EF[i+1]-EF[i])).

This is the magnitude of the secant slope from the definition of the energy
derivative, with units m^-3/eV. It is paired with EF[i], matching the source,
not moved to the interval midpoint. The final interval has no successor
and returns NaN. Repeated energies or nonfinite endpoints also give NaN;
there is no sorting or bridging of gaps.

Model AL references n(E)!K, where K9=-LINEST(J9:J11,A9:A11). J contains
the model's existing excess total holes. For each three-point window,

    D_i = -sum((E-mean(E))*(p_s-mean(p_s))) / sum((E-mean(E))^2).

The fitted slope follows by minimizing squared residuals with an intercept.
The result is assigned to the first row, as in the source. For equally
spaced energy it is the negative centered secant about the middle row;
the first-row assignment is retained deliberately. The final two rows are
NaN, without extending the grid or inventing boundary values.

For an absolute hole population A_p=integral G(E)*(1-f) dE, differentiating
the Fermi occupation gives -dA_p/dEF=integral G(E)*f*(1-f)/kT dE. Subtracting
a constant equilibrium population does not alter this derivative. Thus
the model's negative slope has the DOS units and is a thermally broadened
DOS response, not generally identical to G(E). The measured absolute slope
retains noise and sign reversals in the underlying data; it does not prove
that every row is a physically valid DOS estimate.

## Series and axis mapping

| Legend | X / Y source | Calculated content | Y axis |
|---|---|---|---|
| g(E) | AQ / AR | existing total DOS versus its energy grid | left |
| dp/dEF | S / P | measured derivative versus measured EF | left |
| ntm | AO / AI | model trapped electrons versus model EF | right |
| dpm/dEF | AO / AL | negative fitted model total-hole slope | left |
| ptm | AO / AJ | model trapped holes | left |
| pfm | AO / AH | model free holes | left |
| nfm | AO / AG | model free electrons | left |
| pt (m-3) | S / M | measured trapped holes | right |
| pf (m-3) | S / K | measured free holes | right |

AQ restores the absolute DOS energy; AR is the total DOS, not just the
trap component. Its broken source name reference is replaced by g(E).
All nine entries retain their legend order. The source places several
concentration series on the derivative-labelled left axis; that grouping
is preserved and is not a conversion of concentration into DOS units.
Likewise the n-labelled concentration axis includes p data.

Both Y axes are logarithmic. Both source X quantities are absolute energy;
the implementation aligns them on one common linear X scale rather than
allowing the hidden secondary X autoscale to misalign the energy coordinates.
The retained X label reads Energy, E, delta EF, though the actual values
are absolute E/EF. Energy bounds follow finite measured EF. Left Y bounds
follow the measured derivative; right Y bounds follow measured pt and pf.
The full calculated series remain present and are clipped by display bounds.
No fixed source bounds or unequal/overextended source ranges are copied.

## Integration and validation

Backend: density_energy_derivative. Workflow: run_density_derivatives.
Notebook: density-energy-plot, registered as density_energy for export.
Preview: outputs/density_energy.png.

Subsequent explicit user display setting: X limits -5.8 to -2.8 eV,
applied in the notebook to the shared energy axis. This overrides the
measured-energy display window only; calculation grids remain unchanged.

Subsequent user display setting: both logarithmic Y axes use 1e11 to 1e29.

At the user's request, extend the displayed DOS band edges below 1e12:
include the exact Ev/Ec coordinates in a separate plotting-only energy
array and retain positive DOS-floor endpoints below the visible axis.
The parabolic band contribution vanishes at its edge; the existing floor
only makes that endpoint drawable on a log scale. No finite DOS below the
display floor is inferred from this rendering. The occupation integration
continues to use the original model_energy grid unchanged.

Executed the notebook with unchanged settings, omitting the export call:
203 finite measured derivatives and 2999 finite model derivatives, nine
plotted series, linear X/log Y scales and rendered output. All 72 tests pass,
including independent constant, affine and quadratic derivative cases,
sign checks, missing intervals and repeated-energy handling.

Separately evaluated the routines against saved source input columns as
a read-only comparison (not notebook inputs). All 209 comparable nonzero
measured slopes agree exactly; 2832 comparable nonzero model slopes agree
within 7.8e-16 relative difference. This verifies formula reproduction;
the notebook's actual values can differ because its upstream settings and
calculations remain those of this project.
