# M5: Model effective mobility

The user requested M5 on 2026-09-14. Implement the effective mobility
calculation with the existing microscopic mobility held fixed. The parenthetical
optimization in Supplementary Note 2 does not override the project's explicit
parameter restriction.

Source: [Supplementary Information](../data/42005_2025_2202_MOESM2_ESM.pdf),
p. 7, M5, and [main article](../data/s42005025022021.pdf), equation (4).
The main article's p. 5 identity `mu_eff = mu_0 Theta` was checked visually.

## Derivation and limits

For steady, single-carrier drift transport with immobile trapped charge,
neglect diffusion and recombination. In magnitude, the local current is
`J = e mu_0 p_f F`. Define an effective mobility per total carrier density
by `J = e mu_eff p_s F`. Equating the two gives
`mu_eff = mu_0 p_f/p_s = mu_0 Theta`. The corresponding electron expression
uses `n_f/n_s`; these are separate single-carrier descriptions, not a
combined bipolar mobility. Each ratio is dimensionless, leaving mobility
in m² V⁻¹ s⁻¹. This local definition requires no contact boundary condition.
The Mott–Gurney check additionally assumes ideal injecting contact `F(0)=0`,
uniform permittivity, constant mobility, and negligible equilibrium charge.
Combining `dF/dx=e p_s/epsilon` with drift and integrating across thickness
L yields `J=(9/8) epsilon mu_eff V²/L³`. For physically consistent
populations, no trapping gives `mu_eff=mu_0`, equal free and trapped
populations give `mu_eff=mu_0/2`, and a vanishing free fraction gives zero
mobility.

The retained M2 convention instead supplies `Theta=abs(p_f/p_s)` with
absolute free and excess total populations. Thus the calculation is the
algebraic continuation of the mobility identity using those supplied
ratios; a physical fraction interpretation is not justified everywhere.
Values above one and equilibrium singularities remain unchanged.
See [M2 conventions](M2_excel_conventions.md) for the unresolved density
definitions. No clipping, revised density definitions, or fitting is introduced.

## Data flow and validation

`run_model_mobility` consumes the M2 sweep and its two theta arrays and the
existing `params['mu_0']`. It returns independent energy and mobility arrays.
The same supplied microscopic mobility is applied separately to both
carrier types, without introducing another material parameter.
The notebook plots the two mobilities against the existing quasi-Fermi
energy coordinate and registers the figure with the existing export step.
M4 voltage/current construction is now present in the notebook; see
[M4 derivation and limitations](M4_model_current.md). M5 keeps its energy-domain
plot and its existing population ratio. It is not the voltage-domain
comparison shown in article Fig. 3.

Validation: all 69 existing tests (plus four subtests) pass, including the
independent Mott–Gurney check of equation (4). Executing the notebook's
calculation and plotting cells succeeds for all 3001 model rows with
`mu_0=0.0027 m² V⁻¹ s⁻¹`. Existing exports were left untouched. This
checks integration, not physical validity of the retained M2 ratios.
