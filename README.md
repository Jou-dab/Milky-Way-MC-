# Monte Carlo Simulation of Milky Way Stars

A synthetic stellar population of the Milky Way built with Monte Carlo sampling. The model combines a Kroupa IMF, a constant star formation rate, stellar remnants (white dwarfs, neutron stars and black holes) and chemical evolution (closed box and infall models). The metallicity distribution function (MDF) of the model G dwarfs is compared with observations of the solar neighbourhood to study the G-dwarf problem.

All the code is contained in a single Jupyter notebook
## How to run

```bash
pip install numpy scipy matplotlib astroquery jupyter
```

- An internet connection is required: the observational catalogues are downloaded from VizieR.
- All figures are saved in the `figs/` folder (it must exist before running).
- All random numbers come from `np.random.default_rng(seed)` with `seed = 42`, so the results are reproducible. The convergence study uses seeds `seed + 1 … seed + 10`.
- To run with a different number of stars, change `N` in the parameter cell. The function `simulate(N, seed)` runs the full population for any `N` and seed.

## Notebook structure

| Section | Content |
|---|---|
| Model Parameters | All physical and numerical parameters and the random seed |
| Part 1: Kroupa IMF | Kroupa (2001) IMF, inverse transform sampling and validation against the theoretical probability of each mass segment |
| Part 2: Star Formation and Main Sequence | Constant SFR, `t_birth ~ U(0, 10 Gyr)`, main-sequence lifetime `t_MS = 10 Gyr (M/Msun)^-2.5` and classification (BD / MS / evolved) |
| Part 3: Stellar Remnants | WD masses from Miller et al. (2026) and NS/BH masses from Raithel et al. (2018); initial–final mass relation and distributions of initial mass, final mass and age by type |
| Part 4: Closed Box Model | Gas fraction `mu(t)`, metallicity `Z = y ln(1/mu)`, validation against the analytic prediction and age–metallicity relation |
| Part 5: Observations | Haywood (2001) and GCS III (Holmberg et al. 2009) [Fe/H] from VizieR |
| Part 6: G-Dwarf Problem | G-dwarf MDF (MS stars with 0.8 < M < 1.05 Msun) compared with observations; metal-poor fraction [Fe/H] < −0.5 and MDF peak |
| Part 7: Infall Model | Exponential infall of pristine gas (`Z_in = 0`), gas and metal equations, calibration of `nu` and `y_eff`, and MDFs for several values of `tau` |
| Part 8: Convergence | MS/WD/NS/BH fractions for `N = 10^2 … 10^6` with 10 seeds, compared with the binomial error |

## Main parameters

| Parameter | Value | Description |
|---|---|---|
| `seed` | 42 | Random seed |
| `N` | 10^6 | Number of objects in the main simulation |
| `mass_ranges`, `alphas` | [0.01, 0.08, 0.5, 100], [0.3, 1.3, 2.3] | Kroupa IMF break masses and slopes |
| `T_MW` | 10 Gyr | Age of the Galaxy |
| `T_MS_SUN`, `MS_EXP` | 10 Gyr, −2.5 | Main-sequence lifetime parameters |
| `M_H_BURN` | 0.08 Msun | Minimum mass for hydrogen burning |
| `M_WD_MAX` | 9 Msun | Maximum initial mass for WD formation |
| `F_EJ` | 0.9 | Fraction of the envelope ejected when a BH forms (Raithel+18) |
| `M_BH_BRANCH_II` | 45 Msun | Initial mass at which BH Branch II begins |
| `R_NS_KM`, `GM_SUN_C2_KM` | 12 km, 1.4766 km | NS radius and G Msun / c², used for the NS gravitational mass |
| `NS_EDGES`, `NS_PROB` | Table 3 of Raithel+18 | NS formation probability in each initial-mass interval |
| `Z_SUN`, `MU_F` | 0.014, 0.1 | Solar metallicity and final gas fraction |
| `Y_EFF` | `Z_SUN / ln(1/MU_F)` = 0.00608 | Effective yield of the closed box |
| `MH_FLOOR` | −3 | [M/H] assigned to stars with Z = 0 |
| `MH_CUT` | −0.5 | Metal-poor threshold |
| `G_MIN`, `G_MAX` | 0.8, 1.05 Msun | G-dwarf mass range |
| `DBIN` | 0.1 dex | MDF bin width |
| `TAU_INFALL` | [1, 3, 5, 7, 9] Gyr | Infall timescales |
| `Z_IN`, `M_GAS0` | 0, 0 | Metallicity of the infalling gas and initial gas mass |
| `STEP` | 5000 | Time steps of the Euler integration |
| `N_CONV`, `N_SEEDS` | 10^2 … 10^6, 10 | Sample sizes and number of seeds for the convergence study |

## Functions

| Function | Description |
|---|---|
| `kroupa_imf(m)` | Unnormalized Kroupa IMF |
| `sample_kroupa(n, rng)` | Samples initial masses from the IMF (inverse transform) |
| `sample_birth_times(n, rng)` | Birth times for a constant SFR |
| `ms_lifetime(M)` | Main-sequence lifetime |
| `wd_mass_miller(Mi)` | WD mass, Miller+26 Full Sample Fit 1 (Eqs. 6–9) |
| `ns_probability(Mi)` | NS formation probability (Raithel+18, Table 3) |
| `ns_baryonic_mass(Mi, rng)` | NS baryonic mass for each branch (Raithel+18, Eqs. 11–16) |
| `ns_gravitational_mass(Mb)` | Converts baryonic into gravitational mass (Raithel+18, Eqs. 20–21) |
| `bh_mass_raithel(Mi)` | BH mass (Raithel+18, Eqs. 1–4) |
| `remnant_from_mass(Mi, rng)` | Remnant type and mass of evolved stars |
| `classify_population(M, t_age, rng)` | Classifies all objects as BD, MS, WD, NS or BH |
| `gas_fraction(t)`, `closed_box(t)` | Closed-box gas fraction and metallicity |
| `Z_to_MH(Z)` | Converts Z into [M/H] |
| `mdf`, `frac_below`, `binomial_error` | Normalized MDF, metal-poor fraction and binomial error |
| `integrate_infall`, `build_infall_model`, `check_model` | Integrates, calibrates and checks the infall model |
| `sample_birth_times_from_sfr(model, n, rng)` | Samples birth times from the infall SFR `psi(t)` |
| `simulate(N, seed)` | Runs the full simulation for the convergence study |

## Assumptions

- The IMF is the same at all times, and stellar mass is independent of birth time.
- Constant star formation rate in the closed box model.
- Brown dwarfs (M < 0.08 Msun) are excluded from the stellar fractions.
- All stars are treated as single stars (no binary interactions). The remnant relations assume solar metallicity.
- The Miller relation is extrapolated between 8.39 and 9 Msun (maximum WD mass 1.32 Msun). Raithel Branch I is extrapolated between 40 and 45 Msun.
- NS or BH formation is drawn at random with the probabilities of Raithel et al. (2018). In the intervals 21.7–25.2 and 27.5–60 Msun, where Raithel define no NS branch, all stars form BHs. Fallback BHs are not treated separately.
- NS masses are gravitational masses (R = 12 km). BH masses use a single ejection fraction `f_ej = 0.9`.
- Chemical evolution follows a one-zone model with instantaneous recycling and a well-mixed gas reservoir. Galactic outflows are not included. In the infall model the infalling gas is pristine (`Z_in = 0`).
- Model [M/H] is compared directly with observed [Fe/H] ([M/H] ≈ [Fe/H]).

## Main results

| Sample | f([M/H] < −0.5) |
|---|---|
| Closed box | 0.573 ± 0.003 |
| Infall, tau = 1 Gyr | 0.561 ± 0.003 |
| Infall, tau = 3 Gyr | 0.326 ± 0.003 |
| Infall, tau = 5 Gyr | 0.160 ± 0.002 |
| Infall, tau = 7 Gyr | 0.096 ± 0.002 |
| Infall, tau = 9 Gyr | 0.071 ± 0.002 |
| Haywood (2001) | 0.069 ± 0.017 |
| GCS (Holmberg et al. 2009) | 0.060 ± 0.002 |

The closed box predicts about 8 times more metal-poor G dwarfs than observed (G-dwarf problem). A slow infall of pristine gas (tau ≳ 7–9 Gyr) brings the model close to the observations.

## References

- Haywood, M. 2001, MNRAS, 325, 1365
- Holmberg, J., Nordström, B. & Andersen, J. 2009, A&A, 501, 941
- Kroupa, P. 2001, MNRAS, 322, 231
- Miller, D. R., Caiazzo, I., Heyl, J. et al. 2026, ApJ, 996, 69
- Raithel, C. A., Sukhbold, T. & Özel, F. 2018, ApJ, 856, 35