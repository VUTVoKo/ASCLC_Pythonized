# Fermi energy versus voltage

User-authorized implementation of the MODEL!S formula, with both ABS
operations retained. No physical, processing or grid parameters changed.

    EF = Ev + abs(kB*T/e * ln(abs(pf/Nv)))

The source multiplier S7 is 1; no tuning factor is introduced. The backend
entry point is `fermi_level_absolute_shift`. The existing signed inversion
remains available unchanged. Positive nondegenerate populations satisfy
pf=Nv*exp(-(EF-Ev)/kT), which derives the inversion. Negative populations
produce the requested absolute-value result but are flagged invalid for
physical interpretation; above-Nv populations are folded above Ev and
flagged outside the nondegenerate regime. Zero/nonfinite input and invalid
temperatures return NaN rather than a finite energy or a spreadsheet error.

`run_measured_fermi_level` uses the existing carrier intervals and configured
experiment temperature (299 K), with Nv calculated from the existing hole
mass. It does not import source temperatures or cached Nv. This preserves
current settings; it is not a promise of numerical identity to MODEL!S,
which references its own measured-temperature column and Nv cell.

The notebook calculates energies in the former M3 section, then plots the
comparison after model voltage becomes available. Measured EF is paired
with carriers.U. Model absolute EF is paired with model_current.U_p, tracing
MODEL!Z to j(U)!C and MODEL!AO to the restored absolute energy. EF0 remains
the configured reference. Both axes are linear; limits are derived from
measured energies/voltages and the EF0 reference, not the extreme model
sweep. Labels use EF, correcting the source's EF-Ec series-header mismatch.
The figure is registered for existing export as `fermi_energy`.

Validation: executed notebook calculation and plotting cells without the
export call; checked inverse recovery of measured pf, both ABS branches,
invalid-row gaps, model hole-voltage/energy pairing and linear scales.
All 69 existing unit tests pass. There are 204 valid nondegenerate measured
rows, spanning -5.0468423872 to -4.8144178687 eV. Preview:
`outputs/fermi_energy.png`.
