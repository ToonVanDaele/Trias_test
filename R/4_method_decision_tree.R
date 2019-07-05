# Method decision tree

# The decision wether a species is emerging or not
# is based on multiple simple decision rules.

spDT <- function(df){

  spec <- df[[1,1]]
  lyear <- max(df$year)
  ptitle <- paste0("/dt/", spec, "_DT_",lyear)

  # One value only  -> always emerging ("3")
  # If the 'time series' consists of only one value
  # it is always > 0 (due to preprocessing)

  if (nrow(df) == 1) {
    df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "DT", em = "3")
    return(list(df = df, em = df_em))
  }

  # Two values
  # Second value > 0  -> emerging ("3")
  # Second value = 0  -> possibly emerging ("2")
  if (nrow(df) == 2) {
    if (df[[2,"ncells"]] > 0)
    {df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "DT", em = "3")
      return(list(df = df, em = df_em))
    }else{
      df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "DT", em = "2")
      return(list(df = df, em = df_em))
    }
  }

  # Three or more values
  # first value > 0 and all other = 0 -> not emerging "-3

  # last value > 120% of the mean of the previous values -> "3"
  # last value between 80% and 120% of the mean of previous values -> "2"

  if (nrow(df) >= 3){
    if (sum(df$ncells[2:nrow(df)]) == 0) {
      df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "DT", em = "-3")
      return(list(df = df, em = df_em))
    }
  }

  if (nrow(df) >= 3) {
    pmean <- mean(filter(df, year < lyear)$ncells)
    if (df[df$year == lyear, "ncells"] > pmean * 1.2) {
      df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = "3")
      return(list(df = df, em = df_em))
    }else{
      if (df[df$year == lyear, "ncells"] > pmean * 0.8) {
        df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = "2")
        return(list(df = df, em = df_em))
      }else{
        df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = "0")
        return(list(df = df, em = df_em))
      }
    }
  }



  # More than one value in the time series
  # The value of the last year is > than the year before
  # (Can be fine tuned, example: more than 10 % increase)
  # if (nrow(df) > 1) {
  #
  #   temp <- df %>%
  #   arrange(year) %>%
  #   .$ncells
  #
  #   if (nth(temp, -1) - nth(temp, -2) > 0) {
  #     out <- TRUE
  #   }else{out <- FALSE}
  # }
  #
  # df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear,  method_em = "DT", em = out)
  # outlist <- list(df = df, em = df_em)
  # return(outlist)
}

