set.seed(1)
anl_local <- data.frame(
  USUBJID = c(paste("id", seq(1, 12), sep = "")),
  ARMCD = c(rep("ARM A", 6), rep("ARM B", 6)),
  SEX = c(rep("Female", 6), rep("Male", 6)),
  AVAL = as.numeric(sample(seq(1, 5), 12, replace = TRUE)),
  stringsAsFactors = TRUE
)

adsl_local <- data.frame(
  USUBJID = c(paste("id", seq(1, 12), sep = "")),
  ARMCD = c(rep("ARM A", 6), rep("ARM B", 6)),
  SEX = c(rep("Female", 6), rep("Male", 6)),
  stringsAsFactors = TRUE
)

testthat::test_that("s_count_patients_sum_exposure works as expected", {
  df <- anl_local
  adsl <- adsl_local
  result <- s_count_patients_sum_exposure(df = df, .N_col = nrow(adsl)) # nolintr

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})


testthat::test_that("summarize_patients_exposure_in_cols works well with default arguments", {
  df <- anl_local
  adsl <- adsl_local

  result <- basic_table() %>%
    split_cols_by("ARMCD", split_fun = add_overall_level("Total", first = FALSE)) %>%
    summarize_patients_exposure_in_cols(var = "AVAL", col_split = TRUE) %>%
    split_rows_by("SEX") %>%
    summarize_patients_exposure_in_cols(var = "AVAL", col_split = FALSE) %>%
    build_table(df = df, alt_counts_df = adsl)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("summarize_patients_exposure_in_cols works well with custom arguments", {
  df <- anl_local
  adsl <- adsl_local

  result <- basic_table() %>%
    split_cols_by("ARMCD", split_fun = add_overall_level("Total", first = FALSE)) %>%
    summarize_patients_exposure_in_cols(
      var = "AVAL",
      col_split = TRUE,
      custom_label = "xyz",
      .stats = "sum_exposure"
    ) %>%
    split_rows_by("SEX") %>%
    summarize_patients_exposure_in_cols(
      var = "AVAL",
      col_split = FALSE,
      .stats = "sum_exposure"
    ) %>%
    build_table(df = df, alt_counts_df = adsl)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that(
  "summarize_patients_exposure_in_cols returns correct column label when no variable split and only one statistic",
  code = {
    df <- anl_local
    adsl <- adsl_local

    table <- basic_table() %>%
      summarize_patients_exposure_in_cols(
        var = "AVAL",
        col_split = TRUE,
        custom_label = "xyz",
        .stats = "n_patients"
      ) %>%
      build_table(df = df, alt_counts_df = adsl)

    invisible(capture.output(result <- col_paths_summary(table)$label))

    res <- testthat::expect_silent(result)
    testthat::expect_snapshot(res)
  }
)
