testthat::test_that("format_fraction works with healthy inputs", {
  result <- format_fraction(c(num = 2, denom = 3))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction works with 0 numerator input", {
  result <- format_fraction(c(num = 0L, denom = 3L))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction_fixed_dp works with healthy inputs", {
  result <- format_fraction_fixed_dp(c(num = 2, denom = 3))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction_fixed_dp works with whole number percentages", {
  result <- format_fraction_fixed_dp(c(num = 2, denom = 8))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction_fixed_dp works with 0 numerator input", {
  result <- format_fraction_fixed_dp(c(num = 0L, denom = 3L))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_count_fraction works with healthy inputs", {
  result <- format_count_fraction(c(2, 0.6667))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_count_fraction works with count of 0", {
  result <- format_count_fraction(c(0, 0))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_count_fraction_fixed_dp works with healthy inputs", {
  result <- format_count_fraction_fixed_dp(c(2, 0.5))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_count_fraction_fixed_dp works with healthy inputs", {
  result <- format_count_fraction_fixed_dp(c(2, 0.6667))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_count_fraction_fixed_dp works with count of 0", {
  result <- format_count_fraction_fixed_dp(c(0, 0))

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction fails with bad inputs", {
  x <- list(
    c(num = c(1L, 2L, 3L), denom = 5L),
    c(num = NA_integer_, denom = 2L)
  )
  for (i in x) {
    testthat::expect_error(format_fraction(i))
  }
})

testthat::test_that("format_xx works with easy inputs", {
  test <- list(c(1.658, 0.5761), c(1e1, 785.6))
  z <- format_xx("xx (xx.x)")
  result <- sapply(test, z)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_fraction_threshold works with easy inputs", {
  test <- list(c(100, 0.1), c(10, 0.01), c(0, 0))
  format_fun <- format_fraction_threshold(0.02)
  result <- sapply(test, format_fun)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("h_get_format_threshold works with easy inputs", {
  # Test default.
  result <- h_get_format_threshold()

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)

  # Test non-default value.
  result <- h_get_format_threshold(1L)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("h_format_threshold works with easy inputs", {
  test <- c(0.782, 0.127, Inf, 0, 0.009, NA)
  result <- sapply(test, h_format_threshold)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_extreme_values works with easy inputs", {
  test <- c(0.127, Inf, 0, 0.009, NA)
  format_fun <- format_extreme_values(2L)
  result <- sapply(test, format_fun)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})

testthat::test_that("format_extreme_values_ci works with easy inputs", {
  test <- list(c(0.127, Inf), c(0, 0.009), c(NA, NA))
  format_fun <- format_extreme_values_ci(2L)
  result <- sapply(test, format_fun)

  res <- testthat::expect_silent(result)
  testthat::expect_snapshot(res)
})
