# Classify regions as hyper / hypo / inconclusive

Applies the BREAD decision rule to the output of
[`brms::posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html):

- `hypermethylated` if `p_gt_delta >= prob_cutoff`

- `hypomethylated` if `p_lt_neg_delta >= prob_cutoff`

- `inconclusive` otherwise

## Usage

``` r
classify_regions(post, delta = 0.1, prob_cutoff = 0.95)
```

## Arguments

- post:

  Output of
  [`brms::posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html).

- delta:

  Effect-size threshold used for the rule. Default `0.10`. Stored as an
  attribute; does not re-evaluate the posterior probabilities (those
  must have been computed at this same `delta` upstream).

- prob_cutoff:

  Posterior probability cutoff. Default `0.95`.

## Value

The input `data.frame` with an added `classification` factor column
(levels: `hypermethylated`, `hypomethylated`, `inconclusive`).
Attributes `delta` and `prob_cutoff` are updated.

## Details

In the rare case that both probabilities exceed the cutoff (only
possible for very low `prob_cutoff`), the region is assigned to
whichever side has the larger posterior probability.
