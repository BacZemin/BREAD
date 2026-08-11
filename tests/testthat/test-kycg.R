# bread_kycg() reaches out to KnowYourCG reference databases, so these tests
# deliberately cover only the validation that happens before any network or
# annotation-hub access.

.toy_fit <- function() {
  fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
}

test_that("a non-BreadFit is rejected before any database access", {
  expect_error(bread_kycg(data.frame(x = 1)), "BreadFit")
  expect_error(bread_kycg(NULL), "BreadFit")
})

test_that("`which` and `platform` are matched against their allowed values", {
  skip_if_not_installed("knowYourCG")
  fit <- .toy_fit()
  expect_error(bread_kycg(fit, platform = "NotAPlatform"))
  expect_error(bread_kycg(fit, which = "sideways"))
})

test_that("the fit exposes the mapping columns bread_kycg() reads", {
  # Not a network test: pins the columns bread_kycg() depends on, so a
  # refactor of map_probes_to_features() cannot silently break it.
  fit <- .toy_fit()
  expect_true(all(c("probe_id", "region_id") %in% colnames(fit@mapping)))
  expect_type(fit@mapping$probe_id, "character")
})

# The real group titles as registered by knowYourCG (verified against the
# installed release). Used as a fixture so default-database selection is
# testable without touching ExperimentHub.
.KYCG_TITLES_MM285 <- c(
  "KYCG.MM285.chromHMM.20210210",
  "KYCG.MM285.chromosome.mm10.20210630",
  "KYCG.MM285.designGroup.20210210",
  "KYCG.MM285.HMconsensus.20220116",
  "KYCG.MM285.Mask.20220123",
  "KYCG.MM285.metagene.20220126",
  "KYCG.MM285.probeType.20210630",
  "KYCG.MM285.seqContext.20210630",
  "KYCG.MM285.seqContextN.20210630",
  "KYCG.MM285.TFBSconsensus.20220116",
  "KYCG.MM285.tissueSignature.20211211"
)

test_that("mouse default databases actually match the real MM285 titles", {
  # The regression: the old pattern required a literal "." after "TFBS", so
  # it could never match "KYCG.MM285.TFBSconsensus.20220116" and mouse users
  # silently received an empty data.frame.
  dbs <- BREAD:::.kycg_default_dbs("MM285", .KYCG_TITLES_MM285)

  expect_true("KYCG.MM285.TFBSconsensus.20220116" %in% dbs)
  expect_true("KYCG.MM285.chromHMM.20210210"      %in% dbs)
  expect_true("KYCG.MM285.HMconsensus.20220116"   %in% dbs)
  expect_length(dbs, 6L)
})

test_that("technical annotation groups are deliberately excluded", {
  dbs <- BREAD:::.kycg_default_dbs("MM285", .KYCG_TITLES_MM285)
  for (junk in c("Mask", "chromosome", "probeType", "seqContext")) {
    expect_false(any(grepl(junk, dbs, fixed = TRUE)), info = junk)
  }
})

test_that("default selection returns nothing for an unlisted platform", {
  expect_length(BREAD:::.kycg_default_dbs("EPIC", .KYCG_TITLES_MM285), 0L)
  expect_length(BREAD:::.kycg_default_dbs("MM285", character(0)), 0L)
})

test_that("human platforms select the documented families", {
  titles <- c("KYCG.EPIC.TFBS.20210210", "KYCG.EPIC.chromHMM.20211020",
              "KYCG.EPIC.CGI.20210713",  "KYCG.EPIC.Mask.20220123")
  dbs <- BREAD:::.kycg_default_dbs("EPIC", titles)
  expect_length(dbs, 3L)
  expect_false(any(grepl("Mask", dbs, fixed = TRUE)))
})
