# BREAD 0.0.0.9000 (development)

Initial MVP — milestone 1.

## New features

* `fit_bread()` — end-to-end pipeline: validate → map → summarize → fit
  Bayesian region-level model → classify. Default backend is a conjugate
  Normal-Inverse-Gamma posterior computed analytically (no MCMC).
* `backend = "brms"` — optional full MCMC backend via Stan; compiles the
  model once, then reuses it across regions with `update()`. Shares the
  same downstream posterior/classification path as the conjugate backend.
* `bread_prior()` — constructor for weakly-informative priors.
* `validate_bread_input()`, `map_probes_to_features()`, `summarize_features()`,
  `posterior_summary()`, `classify_regions()` — exported pipeline helpers.
* `plot_region_posterior()`, `plot_region_data()`, `plot_feature_set()` —
  publication-oriented `ggplot2` helpers for BreadFit output.
* `bread_kycg()` — KnowYourCG enrichment on the probes in hyper- or
  hypo-classified regions via `knowYourCG::testEnrichment()`.
* `bread_colors()` — MetBrewer `Cross` palette embedded as hex (no runtime
  dependency on the MetBrewer package).
* `BreadFit` / `BreadResults` S4 classes with accessors `results()`,
  `classifications()`, `posterior_draws()`, and `show()`.

## Known scope

* `backend = "cmdstanr"` and `mode = "hierarchical"` are scaffolded with
  "planned for a later release" errors.
* Both backends fit each region independently (summary mode). Partial
  pooling across regions and class-level priors are deferred to M3.
