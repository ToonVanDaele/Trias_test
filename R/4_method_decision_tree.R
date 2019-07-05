# Method decision tree

# The decision wether a species is emerging or not
# is based on multiple simple decision rules.

spDT <- function(df){

  # If there is only a single value in the time series (> 0)
  # (If there is only one value it should always be > 0 (preprocessing))
  lyear <- max(df$year)

  ptitle <- paste0("/dt/", df[[1,1]], "_DT_",lyear)
  #g <- plot_ts(df, ptitle, printplot = FALSE)

  if (nrow(df) == 1) {
    out <- "3"
  }

  # More than one value in the time series
  # The value of the last year is > than the year before
  # (Can be fine tuned, example: more than 10 % increase)
  if (nrow(df) > 1) {

    temp <- df %>%
    arrange(year) %>%
    .$ncells

    if (nth(temp, -1) - nth(temp, -2) > 0) {
      out <- TRUE
    }else{out <- FALSE}
  }
  df$fit <- df$ucl <- df$lcl <- NA  # decision has no fit, ucl or lcl
  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear, em = out)
  outlist <- list(df = df, em = df_em)
  return(outlist)
}

