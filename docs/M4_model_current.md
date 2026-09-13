# M4: Model voltage and current density

Implemented at the user's request on 2026-09-14. All existing material,
sample, model, and numerical parameters remain fixed. The notebook carries
the M2 energy rows into separate electron and hole curves through
`run_model_current`, plots them alongside measurements, and registers the
figure with the existing export step.

## Derivation

Assume steady one-dimensional single-carrier drift, negligible diffusion,
recombination and ionic current, constant microscopic mobility and uniform
permittivity epsilon. Use field and current magnitudes for injection, with
an ideal injecting contact F(0)=0 and thickness L. Poisson and drift give

    dF/dx = e q(x)/epsilon,    J = e mu_0 f(x) F(x),    dJ/dx = 0.

Here q is total space-charge number density, including mobile and immobile
charge, and f is mobile density, both in m^-3. Adopt the power profile
F(x)=F_L (x/L)^(1-gamma), 0<gamma<1. This is an additional closure, not
an exact solution for arbitrary energy-dependent trap occupations.
Integration gives U=F_L L/(2-gamma). At the collecting boundary,
q_L=epsilon (1-gamma)F_L/(e L), hence

    U = e L^2 q_L / [epsilon (1-gamma)(2-gamma)]
    J = e mu_0 f_L (2-gamma) U/L.

These are the backend equations. Their units are V and A/m^2. Eliminating
q_L using theta=f_L/q_L gives article equation (2). For constant theta,
Poisson and drift instead integrate directly to
F^2=2Jx/(epsilon mu_0 theta), yielding gamma=1/2 and
J=(9/8)epsilon mu_0 theta U^2/L^3 (Mott–Gurney).
For a general supplied gamma, constant current requires f(x)=J/(e mu_0 F(x));
the spatial density ratio need not be constant. The occupation sweep does
not solve that spatial problem. Its logarithmic J–U slope need not equal
the supplied gamma.

## Source mapping and unresolved differences

Checked the original PDFs visually: main article
[equations (2)–(4)](../data/s42005025022021.pdf), pp. 2–3, and
[Supplementary Information](../data/42005_2025_2202_MOESM2_ESM.pdf),
p. 7 M4 and p. 13 equation (S14).

The retained source assignment is `j(U)!C6 = MODEL!B30 / MODEL!B9`,
T_t/T (30/299 for the notebook). The printed M4 instead specifies
T/(T_t+T) for T_t >= T and 0.5 for T_t < T. This discrepancy remains
explicit; the implementation does not substitute the printed prescription
or tune gamma. The backend takes gamma explicitly without a default.

The workbook voltage columns `j(U)!B:C` use `n(E)!H:O`, which are M2
n_t=n_s-n_f0 and p_t=p_s-p_f0. Current columns D:E multiply those voltages
by the corresponding absolute free populations E:L and the drift prefactor.
The workflow retains this mapping and the existing backend physical constants.
The workbook's rounded constants are not substituted.

The required Poisson density is total charge; these labelled trapped columns
are not generally that density. Thus the workflow is an algebraic continuation
of the established source convention, not a validated self-consistent spatial
transport solution. Negative outputs are preserved; they do not establish
a physical extraction solution. Zero charge gives zero voltage and current.
Nonfinite density pairs give NaN, and singular/out-of-domain gamma is rejected.
No absolute values, clipping, sorting, resampling or branch summation is added.

In particular, applying equation (4) to a positive M4 branch gives mu_0 f/q,
whereas M5 retains mu_0 abs(f/s). Their discrepancy follows from q=s-f0;
M4 does not change M2 or M5 to force agreement. See
[M2 density definitions](M2_excel_conventions.md).

## Validation

Tests check the integrated power-field voltage, collecting-boundary drift,
charge conservation, Mott–Gurney limit, invalid inputs, signed/zero outputs,
and row-preserving workflow and plotting. Analytical profiles are test
fixtures, not alternate sample parameter configurations. The notebook is
executed through all calculation and plotting cells without running its
export cell, preserving existing output files.
