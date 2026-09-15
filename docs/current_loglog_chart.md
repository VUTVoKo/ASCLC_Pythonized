# Current density versus voltage, log/log

Implemented 2026-09-15 as the next chart after accepted Drift mobility.
Awaiting user review. All existing parameters and numerical grids are fixed.

## Source mapping

Inspected `data/SCLCKopecky.xlsx`, `xl/charts/chart2.xml` and MODEL formulas:

| Series | Source X / Y | Calculated inputs |
|---|---|---|
| j (A/m²) | MODEL E / G | measured V / J |
| jm (pf) | MODEL Z / AC | U_p / J_p |
| −jm (nf) | MODEL Z / AB | U_p / −J_n |
| m = 1 | MODEL I2:I3 / K2:K3 | Ohmic reference |
| m = 2 | MODEL I2:I3 / N2:N3 | Mott–Gurney reference |

Both model series explicitly use Z, the hole voltage, rather than the
electron voltage Y. AB negates AA, which references the electron current
calculated at its own voltage. Thus U_p versus −J_n is a pairing by shared
Fermi energy, not a self-consistent electron J(U_p) solution or a bipolar
current sum. The display negation leaves the signed calculated arrays intact.
Nonpositive or nonfinite log coordinates leave gaps; all model rows remain.

Both axes are logarithmic. Labels are Voltage, U (V), and Current density,
j (A/m²). Limits are captured from positive finite measured points before
adding either model curve or either reference. The expanded mobility bounds
do not apply to this chart.

## Reference equations and source differences

In the low-injection Ohmic limit, take the equilibrium free-hole density
pf0 and a uniform field F=U/L. Drift gives j=e μ0 pf0 U/L, hence slope
d ln(j)/d ln(U)=1. This assumes a maintained equilibrium reservoir and
negligible injected space charge, diffusion and recombination.

For the trap-free single-carrier space-charge limit, steady current and
Poisson's equation give j=e μ0 p F and dF/dx=e p/ε. Eliminating p yields
j=ε μ0 F dF/dx. With ideal injection F(0)=0, integration gives
F²=2jx/(ε μ0), and U=∫F dx implies j=9 ε μ0 U²/(8L³). Its log slope
is 2. These are distinct limiting regimes, not replacements for the
existing model gamma or occupation calculations. Both expressions have
units A/m². The derivation also appears in [M4](M4_model_current.md).

References use the existing μ0, L, εr and model equilibrium pf0. Only their
two plotted endpoints are chosen from the measured voltage display bounds;
a straight line in log/log coordinates represents the exact power law.

Source K multiplies the Ohmic expression by MODEL L3=2; source N divides
the Mott–Gurney expression by O3=0.4. These manual amplitude factors have
no derivation in the governing equations and are not project parameters.
They are retained at the user's instruction on 2026-09-15, because they are
what keeps the M4 curve inside the band the two references mark out. Writing
the M4 hole branch against the bare limits gives

    j/j_Ohmic     = (2-g) * pf/pf0
    j/j_Mott–Gurney = (8/9)(1-g)(2-g)^2 * Theta,   Theta = pf/pt

At the supplied g=T_t/T=0.10033 those prefactors are 1.8997 and 2.8859, so
the bare m=2 reference is crossed wherever Theta>0.347; over the plotted
sweep the curve reaches 1.62 times that reference. The L3 and O3 factors
shift the references by 2 and 2.5, and the curve then stays between them
across the whole displayed range (1.03 to 0.65 of the respective lines).
Both factors are therefore reproduced as constants carrying their source
cells, not as fitted quantities, and neither g nor any model parameter is
changed to accommodate them. Only g=1/2 would make the second prefactor
unity and bound the curve by the unshifted Mott–Gurney line; that value is
the printed M4 prescription for T_t<T but not the retained source
assignment, and adopting it would move the trap-filling rise from 1.7 V to
4.0 V, outside the measured range. See [M4](M4_model_current.md).

Source J2 uses a Boltzmann equilibrium density; the implementation uses the
already calculated equilibrium free-hole density, which agrees with J2 to
0.19 per cent. Reference endpoints evaluated at MODEL I2=1e-3 V and I3=100 V
reproduce the workbook's own K and N values to better than 2e-4 relative,
the residual being its rounded e and eps_0. No cached chart values or
alternate parameter configurations supply plotted values.

## Validation

Execute all notebook calculations and plots with unchanged settings,
excluding the export cell. Check five series, both model mappings against
all original energy rows, log/log scales and equality of display bounds
to the measured-only chart. Verify slopes from endpoint logarithms and
Ohmic conductivity jL/U=e μ0 pf0; independently invert the Mott–Gurney
reference using equation (4) at its analytical gamma=1/2 limit to recover
the supplied μ0. This analytical limit is a validation identity, not a
change to the model gamma. Existing tests also cover that inversion.

`test_asclc_references.py` additionally checks both reference slopes, both
endpoint amplitudes against the workbook cell values, the two factors as the
only departure from the bare equations, and that the M4 hole branch on the
notebook's own parameters and sweep stays between the references everywhere
it is displayed while exceeding the unshifted Mott–Gurney line.
