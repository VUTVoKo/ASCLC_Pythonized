# Theta versus voltage

Implemented measured/model comparison using existing arrays and parameters.
Source chart9.xml pairs MODEL!Z/AN for the model and E/Q for measurements.
Z references j(U)!C (hole voltage); AN references n(E)!N, abs(pf/ps).
Q is measured effective mobility divided by microscopic mobility, equal
to the current default measured theta convention. Existing upstream
processing remains unchanged; no source data are copied into the plot.

The model curve uses model_current.U_p and model.theta_p on the same EF
sweep. It retains the previously documented excess-total convention, not
a newly imposed pf/(pf+pt) formula. The identity mu_eff=mu0*Theta connects
the measured ratio to the current carrier calculation. Values above one
are retained, not numerically capped.

The notebook adds the combined chart after model voltage is available;
it replaces the carrier_fraction export entry. The earlier measured-only
diagnostic remains available. The combined chart has two series and no
extra diagnostic reference lines. Labels distinguish Theta and Theta_m;
the source model series was unnamed. X is linear, Y logarithmic, and
display bounds fit the measured points. Full model arrays remain present
outside that viewport. Nonfinite/nonpositive logarithmic values leave gaps.

Validation: full notebook calculation/plotting run without export; verified
both series, hole-branch X coordinates, abs(pf/ps) Y coordinates and scales.
All 72 existing tests pass. Preview: outputs/carrier_fraction.png.
