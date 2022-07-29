library(scda)
library(rtables)
library(dplyr)
library(tern)

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

result_l_rows <- build_table(l_rows, df = adpc)
result_l_rows

l1 <- l_rows %>%
    summarize_vars_in_cols(
        var = "AVAL",
        col_split = TRUE,
        .stats = c("n", "mean"),
        .formats = c(
            n = "xx.",
            mean = sprintf_format("%.3e")),
        .labels = c(
            n = "n",
            mean = "Mean"),
        .format_na_strs = c(
            n = "Ridiculous",
            mean = "Ridiculous_2"
        )
    )


result1 <- build_table(l1, df = adpc)
result1


