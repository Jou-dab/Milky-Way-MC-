# Monte Carlo Simulation of Milky Way Stars

Synthetic Stellar Population of the Milky Way:

This model uses a Kroupa IMF, constant star formation, stellar remnants and chemical evolution. The metallicity distribution function of G dwarf in this model is compared with observations of the solar neighbourhood to get conclusions about the G-Dwarf problem. All the code is compress in a notebook.

## How to run

The necessary packages to run the notebook are: numpy, scipy, matplotlib, astroquery.

## **Notebook Structure**

| Section | Content |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| Parameters | Model parameters and random seed (`seed = 42`). |
| I.a IMF | Kroupa (2001) IMF, inverse transform sampling, and KS test. |
| I.b SFR | Constant star formation rate, with `t_birth ~ U(0, 10 Gyr)`. Main-sequence lifetime and stellar classification. |
| II Stellar Remnants | WD masses from Miller et al. (2026) and NS/BH masses from Raithel et al. (2018). Initial and final mass distributions. |
| III Chemical Evolution | Closed-box model, gas fraction, metallicity evolution, and age–metallicity relation. |
| IV Observations | Haywood (2001) and GCS III (Holmberg et al. 2009) data from VizieR. |
| V MDF | G-dwarf MDF compared with observations. Fraction of stars with [Fe/H] < −0.5 and the G-dwarf problem. |
| Infall | Gas accretion model with primordial gas (`Z_in = 0`). Calibration of `nu` and `y_eff`, and MDFs for `tau = 2` and `7 Gyr`. |
| VI Convergence | Stellar population fractions for different sample sizes (`N = 10^2 ... 10^6`), using 10 random seeds. |

## **Main Parameters**

| Parameter | Value | Description |
| ------------------------------ | --------------------------------------- | -------------------------------------------- |
| `seed` | 42 | Random seed |
| `N` | 10^6 | Number of stars |
| `mass_ranges`, `alphas` | [0.01, 0.08, 0.5, 100], [0.3, 1.3, 2.3] | Kroupa IMF parameters |
| `T_MW` | 10 Gyr | Galaxy age |
| `T_MS_SUN`, `MS_EXP` | 10 Gyr, −2.5 | Main-sequence lifetime parameters |
| `M_H_BURN` | 0.08 Msun | Minimum mass for hydrogen burning |
| `M_WD_MAX` | 9 Msun | Maximum initial mass for WD formation |
| `F_EJ` | 0.9 | Mass ejection fraction during BH formation (Raithel+18) |
| `M_BH_BRANCH_II` | 45 Msun | Initial mass at which BH Branch II begins |
| `R_NS_KM` | 12 km | NS radius, used to calculate gravitational mass |
| `NS_EDGES`, `NS_PROB` | Table 2 of Raithel+18 | NS formation probabilities |
| `Z_sun`, `MU_F` | 0.014, 0.1 | Solar metallicity and final gas fraction |
| `y_eff` | `Z_sun / ln(1/MU_F)` = 0.00608 | Effective yield |
| `MH_FLOOR` | −3 | Minimum [M/H] value for Z = 0 |
| `M_G_MIN`, `M_G_MAX` | 0.8, 1.05 Msun | G-dwarf mass range |
| `TAU_INFALL`, `Z_IN`, `M_GAS0` | [2, 7] Gyr, 0, 0 | Infall model parameters |
| `N_CONV`, `N_SEEDS` | 10^2 ... 10^6, 10 | Sample sizes and number of seeds for convergence tests |

## **Functions**

| Function | Description |
| ------------------------------------------------------- | --------------------------------------------------- |
| `kroupa_imf(m)` | Kroupa IMF |
| `kroupa_cdf(m)` | Normalised IMF CDF |
| `sample_kroupa(n, rng)` | Samples stellar masses from the IMF |
| `sample_birth_times(n, rng)` | Generates birth times for a constant SFR |
| `ms_lifetime(M)` | Calculates main-sequence lifetime |
| `wd_mass_miller(Mi)` | Calculates WD mass using Miller+26 (Eqs. 6–9) |
| `ns_probability(Mi)` | Calculates NS formation probability (Raithel+18) |
| `ns_baryonic_mass(Mi, rng)` | Samples the baryonic mass of NSs |
| `ns_gravitational_mass(Mb)` | Converts baryonic mass into gravitational mass |
| `bh_mass_raithel(Mi)` | Calculates BH mass using Raithel+18 |
| `remnant_from_mass(Mi, rng)` | Determines remnant type and mass |
| `classify_population(M, t_age, rng)` | Classifies stars as BD, MS, WD, NS or BH |
| `gas_fraction(t)`, `closed_box(t)` | Calculates gas fraction and metallicity over time |
| `Z_to_MH(Z)` | Converts metallicity Z into [M/H] |
| `closed_box_frac_below`, `closed_box_mdf_analytic` | Analytical closed-box predictions |
| `load_feh_catalogs()` | Downloads or loads the Haywood and GCS catalogues |
| `mdf`, `frac_below`, `binomial_error` | Calculates MDFs, fractions and statistical errors |
| `integrate_infall`, `build_infall_model`, `check_model` | Builds and checks the infall model |
| `sample_birth_times_from_sfr(model, n, rng)` | Samples birth times using the model SFR |
| `simulate(N, seed)` | Runs the full simulation |

## **Assumptions**

- The IMF is the same at all times, and stellar mass is independent of birth time.
- All stars are treated as isolated systems, with no binary interactions.
- Remnant mass relations assume solar metallicity.
- The Miller relation is extrapolated between 8.39 and 9 Msun, with a maximum WD mass of 1.32 Msun. The Raithel Branch I relation is extrapolated between 40 and 45 Msun.
- NS or BH formation is randomly assigned using the probabilities from Raithel et al. (2018).
- Chemical evolution follows a one-zone model with instantaneous recycling and a well-mixed gas reservoir. Galactic outflows are not included.
- Model [M/H] is compared directly with observed [Fe/H], assuming [M/H] ≈ [Fe/H].
## References

- Haywood, M. 2001, MNRAS, 325, 1365
- Holmberg, J., Nordström, B. & Andersen, J. 2009, A&A, 501, 941
- Kroupa, P. 2001, MNRAS, 322, 231
- Miller, D. R., Caiazzo, I., Heyl, J. et al. 2026, arXiv:2510.24877
- Raithel, C. A., Sukhbold, T. & Özel, F. 2018, ApJ, 856, 35