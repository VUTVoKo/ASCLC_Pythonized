# Alternative measured series

Two further measured series, for checking that the analysis chain is not tied
to `data/MAPbBr3_S2_dark.csv`.

| File | Sample | Condition | Source figure |
|---|---|---|---|
| `data/MAPbBr3_S2_light.csv` | C/MAPbBr3/C | illuminated | Fig. S17(a) |
| `data/MAPbI3_light.csv` | C/MAPbI3/C | illuminated | Fig. S18(a) |

Both hold 214 rows of `U(V),I(A)` on a uniform 14.04 mV sweep ending at 3.0 V,
matching the step and length of the existing dark series.

## Why these two

Supplementary Note 4 states that the **injection** model — the one implemented
here — "is sufficient to fully describe the behavior of MAPbBr3 single crystals
(dark, after illumination)", while MAPbI3 **under dark conditions** shows a hole
injection barrier so that "the extraction mechanism dominates" and "the
extraction diode model must be implemented".

MAPbI3 dark is therefore outside what this code models. Run through the current
chain it yields gamma ~ 1, theta ~ 166 and negative `p_t`. It was traced and
discarded; the illuminated MAPbI3 series is used instead.

## Provenance and method

The figure panels are raster images (~600x490 px, JPEG), so the points were
recovered by tracing, not read from vector paths:

1. The plot frame is located from the full-width/full-height dark rows and
   columns.
2. Axis calibration is confirmed independently from tick spacing rather than
   assumed. Major y-tick pitch gives 40.75 px/decade over 408 px (10 decades,
   1e-4..1e-14) for Fig. S17(a) and 51 px over 409 px (8 decades, 1e-4..1e-12)
   for Fig. S18(a). Interior x-ticks fall within 0.5 px of the positions implied
   by 0.001-10 V across four decades.
3. The coloured *averaged* markers are masked by hue, separating them from the
   grey raw cloud, and traced one pixel column at a time as the midpoint of the
   marker outline. The legend box is excluded by region.
4. The trace is resampled onto the uniform voltage grid by interpolation in
   log U - log I.

### Validation

The same procedure applied to Fig. S17(b), whose series is already held as
`data/MAPbBr3_S2_dark.csv`, reproduces it. Against a log-V binned mean of those
raw readings the traced curve has a median ratio of 0.9-1.0x across the sweep;
restricted to the signal-dominated region above 2 V, where the raw data is not
dominated by the noise floor, the median ratio is 0.93x with a maximum deviation
of 2.2x. Treat the two series above as accurate to roughly 10-20% in current
where the signal is well above the noise floor, and worse below it.

### Caveats

- These are the **averaged** curves, not raw readings. They arrive already
  smooth, so the trailing-mean step is close to a no-op on them; the existing
  dark series is raw and behaves differently under the same window.
- The legend inside Fig. S17(a) reads "MAPbBr3 dark"; the caption, the current
  magnitude and Fig. S19(a) all identify the panel as the illuminated series.
  The legend text is an error in the source.
- Coverage starts where the averaged markers become continuous: 0.0081 V for
  MAPbBr3 light and 0.0091 V for MAPbI3 light. Nothing is extrapolated below.

## Published parameters

From Table S4. Units are converted to those the notebook uses; the paper's
concentrations are cm^-3 and its mobilities cm^2 V^-1 s^-1.

### MAPbBr3, illuminated

Same sample as the existing dark series, so `material` changes only if the
thickness is taken from the paper (see below).

    params = dict(
        mu_0=1.6e-2,    # m^2 V^-1 s^-1  (160 cm^2 V^-1 s^-1)
        N_t=3.54e16,    # m^-3           (3.54e10 cm^-3)
        E_t=-4.768,     # eV
        T_t=40.0,       # K
        E_F0=...,       # not published per series
    )

### MAPbI3, illuminated

    material = dict(
        L=1.05e-3,      # m
        S=20.10e-6,     # m^2
        eps_r=32.00,
        E_c=-3.93,      # eV
        E_v=-5.43,      # eV
        E_g=1.50,       # eV
        m_eff_p=0.350,
        m_eff_e=0.350,
    )

    params = dict(
        mu_0=2.3e-2,    # m^2 V^-1 s^-1  (230 cm^2 V^-1 s^-1)
        N_t=3.45e16,    # m^-3           (3.45e10 cm^-3)
        E_t=-4.715,     # eV
        T_t=95.0,       # K
        E_F0=...,       # not published per series
    )

Table S4 also gives N_c,v = 5.17e10 cm^-3 for MAPbI3 and T = 300 K, V = 0-3 V at
1 mV s^-1 for every series.

These are the published values, listed so they do not have to be looked up
again. They are not a recommendation to replace any hand-selected parameter.
`E_F0` is given in Table S2 only as a fixed "intrinsic Fermi level", not per
series, so it has no published per-series value to quote.

### Divergence from the existing dark parameters

The existing MAPbBr3 setup does not use the Table S4 values throughout:
`L` is 6.00e-4 m against the paper's 0.80 mm, `E_t` is -4.82 eV against -4.768,
and `N_t` is 4.7e16 m^-3 against 3.03e16. These are hand-selected and stay as
they are; the point here is only that the MAPbBr3-light block above is the
paper's, so mixing it with the existing thickness is a choice to make
deliberately.

## Reference check

Measured-side chain with the published parameters above, at the top of each
sweep, against the `p_t` of Table S4:

| Series | p_t computed | p_t published | ratio |
|---|---|---|---|
| MAPbBr3 dark (existing, hand-selected params) | 1.89e16 | 7.58e15 | 2.49x |
| MAPbBr3 light | 1.98e15 | 8.84e15 | 0.22x |
| MAPbI3 light | 1.10e16 | 8.63e15 | 1.27x |

The new series sit either side of the published values by about the same factor
the existing series already differs by, so they are no worse conditioned than
the baseline. Both produce some rows with theta above 1 and a few negative `p_t`
at low bias, which the existing dark series does not.
