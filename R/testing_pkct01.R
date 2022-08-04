library(scda)
library(rtables)
library(dplyr)
library(tern)

############################
# Example of pkct01 code with only the mean
# This is only using summarize_row_groups and is not transposing the data
############################
adpc <- scda::synthetic_cdisc_data("latest")$adpc
adpc <- adpc %>%
  mutate(NRELTM1 = as.factor(NRELTM1),
         AVALC = as.factor(AVALC)) %>% filter(ACTARM %in% c("A: Drug X")) %>%
  mutate(ACTARM = factor(ACTARM, levels = c("A: Drug X")))

adpc_1 <- adpc %>%
  mutate(NRELTM1 = as.factor(NRELTM1),
         AVAL = AVALC) %>% filter(ACTARM %in% c("A: Drug X")) %>%
  mutate(ACTARM = factor(ACTARM, levels = c("A: Drug X")))

l_rows <- basic_table() %>%
  split_rows_by(
    var = "ACTARM",
    split_label = "Cohort/Treatment",
    label_pos = "topleft") %>%
  split_rows_by(
    var = "VISIT",
    split_label = "Visit",
    label_pos = "topleft") %>%
  split_rows_by(
    var = "NRELTM1",
    split_label = "Norminal Time from First Dose",
    label_pos = "topleft"
  )

l1 <- l_rows %>%
  summarize_vars_in_cols(
    var = "AVAL",
    col_split = TRUE,
    .stats = c("mean"),
    .formats = c(
      mean = sprintf_format("%.3e")),
    .labels = c(
      mean = "Mean"),
    .format_na_strs = c(
      mean = "Mean")
  )


result1 <- build_table(l1, df = adpc)
result1



############################
# Example 2
# This is only using summarize_row_groups and is not transposing the data
############################
s_summary <- function(x) {
  stopifnot(is.numeric(x))

  list(
    n = sum(!is.na(x)),
    mean = mean(x),
    min_max = range(x)
  )
}

a_summary <- make_afun(
  fun = s_summary,
  .formats = c(n = "xx", mean = "xx.xx", min_max = "xx.xx - xx.xx"),
  .labels = c(n = "n", mean = "Mean", min_max = "min - max")
)

a_summary3 <- make_afun(a_summary,
                        .formats = c(mean = "xx.xxx"),
                        .format_na_strs = c(mean = "Ridiculous"))

l <- basic_table() %>% split_cols_by("ARM") %>%
  split_rows_by("ACTARM",
                split_label = "Cohort/Treatment",
                label_pos = "topleft") %>%
  split_rows_by(var = "VISIT",
                split_label = "Visit",
                label_pos = "topleft") %>%
  split_rows_by(var = "NRELTM1",
                split_label = "Norminal Time from First Dose",
                label_pos = "topleft"
  ) %>%
  summarize_row_groups(label_fstr = "%s (n)") %>%
  analyze("AVAL", afun = a_summary3, format = "xx.xx")

tbl <- build_table(l, adpc)
tbl


