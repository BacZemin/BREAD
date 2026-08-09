# BREAD 0.99.0

First release candidate, prepared for Bioconductor submission.

## Breaking changes

0.99.0 is unreleased, so these land without a deprecation period.

* **`classification` now has four levels**, not three:
  `hypermethylated`, `hypomethylated`, **`unchanged`**, `inconclusive`.
  `unchanged` means the posterior is concentrated inside the region of
  practical equivalence — evidence that a region did *not* move, which
  `inconclusive` previously absorbed along with genuinely uninformative
  regions. Membership of `hypermethylated` and `hypomethylated` is
  unchanged for any `prob_cutoff > 0.5`, so directional filters are
  unaffected; only code that counts or filters `inconclusive` sees a
  difference, and it sees a more honest one.
* **`results()` gains six columns**: `prob_rope`, `ref_beta`, `mean_dbeta`,
  `dbeta_lo`, `dbeta_hi`, `delta_beta`.
* **`n_features_in` now counts distinct region IDs**, not `GRanges` ranges.
  When several ranges share a `region_id` — the only way to pin a region to
  an exact probe set — the old counter reported the range count, printing
  e.g. `n_regions: 788 (of 790 input)` for 30 regions built from 790 probes.
  `show()` had the same defect. Results were always correct; the counters
  were not.
* **`fit_bread()`'s first argument is renamed `se` to `x`**, since it now
  also accepts a matrix.
* The `inconclusive` swatch in `bread_colors("classification")` moved from
  olive to neutral grey; `unchanged` took the olive. Absence of information
  should not look like a finding.

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
  hypo-classified regions via `knowYourCG::testEnrichment()`. Default
  knowledgebase selection is now platform-aware: the previous pattern
  required a literal `.` after `TFBS`, so it could never match the real
  MM285 titles (`KYCG.MM285.TFBSconsensus.20220116`) and mouse users got a
  silently empty result. A no-match now warns and lists what *is* available,
  and `mtc_by_group` / `mtc_method` are passed through when the installed
  knowYourCG supports them.
* `refit_bread()` — re-fit or re-threshold an existing `BreadFit` without
  recomputing the probe-to-region mapping or the region summary. Makes
  label-permutation calibration a first-class workflow instead of a reason
  to call `BREAD:::fit_bread_summary()`.
* **Matrix input.** `fit_bread()` and the new `bread_se()` accept a
  probe-by-sample matrix with `colData` plus either `rowRanges` or
  `platform`, as well as sesame's `list(betas =, sampleInfo =)` shape —
  so `openSesame()` output no longer has to be hand-assembled into a
  `SummarizedExperiment` first. The platform is never inferred from probe
  IDs: `cg`-numbers are shared across arrays, and a wrong guess would give
  wrong coordinates silently.
* `bread_delta_beta()` / `bread_delta_m()` — convert effect sizes between
  the M and beta scales, and `results()` now reports a per-region
  `delta_beta`. `delta = 0.10` on the M scale is a beta change of about
  0.017 at mid-methylation and less toward the extremes; `fit_bread()` says
  so once when handed beta-scale input.
* `ci` and `rope_cutoff` are now `fit_bread()` arguments. `ci` was
  previously hardcoded at 0.95, distinct from `prob_cutoff` but not
  reachable.
* A rank-deficient design now warns. The conjugate prior absorbs the
  deficiency rather than erroring, so such fits previously looked normal
  while returning prior-driven estimates for the collinear coefficients.
* `bread_colors()` — MetBrewer `Cross` palette embedded as hex (no runtime
  dependency on the MetBrewer package).
* `BreadFit` / `BreadResults` S4 classes with accessors `results()`,
  `classifications()`, `posterior_draws()`, and `show()`.
* `BreadResults()` — constructor for the `BreadResults` class, with a
  `show()` method. The class previously had no way to build one.
* `posterior_summary()` now accepts a `BreadFit` directly, so callers no
  longer need to reach into the object's `model` slot.
* Every exported object now carries a runnable `@examples` block built on
  the packaged EPICv2 example data.
* Posterior probability columns in `results()` are named `prob_pos`,
  `prob_neg`, `prob_hyper` and `prob_hypo` — deliberately not `p_*`, which
  invites reading them as p-values. They are posterior probabilities of the
  parameter given the data (`prob_hyper` = P(effect > +delta),
  `prob_hypo` = P(effect < -delta)), not tail probabilities of a statistic
  under a null hypothesis.

## Known scope

* `backend = "cmdstanr"` and `mode = "hierarchical"` are scaffolded with
  "planned for a later release" errors.
* The unimplemented `report_feature_set()` stub has been removed;
  `plot_feature_set()` covers the same ground.
* Both backends fit each region independently (summary mode). Partial
  pooling across regions and class-level priors are deferred to M3.
