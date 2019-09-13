#' Method: decision tree

#' Decision of emerging status based on some decision rules.
#' The emergency status based on decision rules is needed for time series
#' that cannot be handled with statistical algorithms (too short or too few observations)

#' First the time series are characterised with some tests (TRUE/FALE)
#' Second the results of these tests are combined to decided the emergency status

#' Only the occupancy (ncells) information is used.

#' @param df data frame with time series
#' @return list with dataframe, emergency status and tests

spDT <- function(df){

  spec <- df[[1,1]]     # species name
  lyear <- max(df$year) # last year
  nb <- nrow(df)        # length of time series
  ptitle <- paste0("/dt/", spec, "_DT_",lyear)

  dt <- rep(FALSE, 6)  # Vector to store results of tests (default = FALSE)

  # Time series with only one value?   -> appearing species  (always > 0)
  if (nb == 1) dt[1] <- TRUE

  # 50% of time series > 0?
  if (sum(df$ncells > 0) / nb >= 0.5) dt[2] <- TRUE

  # Second value = 0  -> possibly emerging ("2")
  #if (df[[2,"ncells"]] > 0) dt[3] <- TRUE

  # 0 since 5 years -> not emerging
  if (sum(df$ncells[max(nb - 4, 0):nb]) == 0) dt[4] <- TRUE

  # > 0 and 5 consecutive years 0 before  -> (re)appearing
  if (nb > 1 & df$ncells[nb] > 0 & sum(df$ncells[max(nb - 5, 0):(nb - 1)]) == 0) dt[5] <- TRUE

  # Maximum of ncells == 1
  if (max(df$ncells) <= 1) dt[6] <- TRUE

  # Increase? Last value > before last value
  if (nb > 1 & df$ncells[nb] > df$ncells[nb - 1]) dt[7] <- TRUE

  #em status codes:
  # 0 = not emerging
  # 1 = appearing
  # 2 = re-appearing
  # 3 = possibly emerging
  # 4 = emerging

  em <- case_when(
    dt[4] == TRUE ~ 0,                 # zeros since > 5 years => not emerging
    dt[1] == TRUE ~ 1,                 # One value > 1 => appearing
    dt[5] == TRUE ~ 2,                 # > 0 after 5 consecutive 0 => re-appearing
    dt[6] == TRUE & dt[2] == TRUE ~ 3,  # possibly emerging
    dt[7] == TRUE ~ 4,                 # Last year increase => emerging
  )

  if (is.na(em)) em <- "0"

  ptitle <- paste0(spec, "_", lyear, "_", em)
  g <- plot_ts(df = df, ptitle = ptitle)

  df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = em)
  outlist <- list(df = df, dt = dt, plot = g, em = df_em)
  return(outlist)
}


## Rule set can be simplyfied by cutting any series with > 4 years 0 and only consider
# the series herafter.



  # One value
  # x         (x > 0)    ->  "appearing"

  # Two values
  # x 0                  ->  "possibly emerging"
  # x y  &  y > x        ->  "emerging"
  # x y  &  y <= x       ->  "possibly emering"

  # Only two values > 0, not consecutive.
  # x 0 y      (y > x)     -> "emerging" (3)
  # x 0 y      (y <= x)    -> "possibly emerging" (2)
  # x 0...0 y              ->
  # x 0...0 y 0...0        ->


  # 1 0.(min 5) 0                   -> not emerging
  # 1 0.(min 5)..0 1                -> (re)appearing species
  # 1 00(max4) 1                    -> not emerging
  # 1 11 1                          -> not emerging

  # laatste waarde > 1 en > voorlaatste waarde

  # Three or more values
  # first value > 0 and all other = 0 -> not emerging "-3


