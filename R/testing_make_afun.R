
s_foo <- function(df, .N_col, a = 1, b = NA) {
  list(
    nrow_df = nrow(df),
    .N_col = .N_col,
    a = a,
    b = b
  )
}

s_foo(iris, 40)

a_foo <- make_afun(s_foo,
                   .formats = c(nrow_df = "xx.xx", ".N_col" = "xx.", a = "xx", b = "xx.x"),
                   .labels = c(nrow_df = "Nrow df", ".N_col" = "n in cols", a = "a value", b = "b value"),
                   .indent_mods = c(nrow_df = 2L, a = 1L),
                   .format_na_strs = c("Test", "Test", "Test", "Test")
)
a_foo(iris, .N_col = 40)

##
# Make a table where you have rtable object and search for format_na_strs when the .format_na_strs
# actually works
# Then use the summarize_vars_in_cols


