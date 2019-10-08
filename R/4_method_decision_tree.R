#' Method: decision tree

#' Decision of emerging status based on some decision rules.
#' The emergency status based on decision rules is needed for time series
#' that cannot be handled with statistical algorithms (too short or too few observations)

#' First the time series are characterised with some tests (TRUE/FALE)
#' Second the results of these tests are combined to decided the emergency status

#' Only the occupancy (ncell) information is used.

#' @param df data frame with time series
#' @return list with dataframe, emergency status and tests

spDT <- function(df){

  spec <- df[[1,1]]     # species name
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  lyear <- max(df$year) # last year
  nb <- nrow(df)        # length of time series
  ptitle <- paste0("/dt/", spec, "_", spn, "_DT_",lyear)
  cat(ptitle, "\n")

  dt <- rep(FALSE, 8)  # Vector to store results of tests (default = FALSE)

  # DT_1: Time series with only one value?   -> appearing species  (always > 0)
  if (nb == 1) dt[1] <- TRUE

  # DT_2: 50% of time series > 0?
  if (sum(df$ncell > 0) / nb >= 0.5) dt[2] <- TRUE

  # DT_3: last value above median value -> possibly emerging ("2")
  if (df$ncell[nb] > median(df$ncell)) dt[3] <- TRUE

  # DT_4: 0 since 5 years -> not emerging
  if (sum(df$ncell[max(nb - 4, 0):nb]) == 0) dt[4] <- TRUE

  # DT_5: > 0 and 5 consecutive years 0 before  -> (re)appearing
  if (nb > 1 && df$ncell[nb] > 0 & sum(df$ncell[max(nb - 5, 0):(nb - 1)]) == 0) dt[5] <- TRUE

  # DT_6: Maximum of ncell == 1
  if (max(df$ncell) <= 1) dt[6] <- TRUE

  # DT_7: Increase? Last value > before last value
  if (nb > 1 && df$ncell[nb] > df$ncell[nb - 1]) dt[7] <- TRUE

  # DT_8: Maximum ever observed?
  if (df$ncell[nb] > max(df$ncell)) dt[8] <- TRUE

  #em status codes:
  # 0 = not emerging
  # 1 = unclear
  # 2 = potentially emerging
  # 3 = emerging
  # 4 = appearing / re-appearing

  em <- case_when(
    dt[4] == TRUE ~ 0,                  # zeros since > 5 years => not emerging
    dt[1] == TRUE ~ 1,                  # One value > 1 => appearing / re-appearing
    dt[5] == TRUE ~ 1,                  # > 0 after 5 consecutive 0 => re-appearing
    dt[6] == TRUE & dt[2] == TRUE ~ 1,  # appearing / re-appearing
    dt[7] == TRUE & dt[8] == FALSE ~ 2, # Last year increase but not maximum => potentially emerging
    dt[3] == TRUE ~ 2,                  # potentially emerging
    dt[8] == TRUE & dt[1] == FALSE ~ 3  # maximum ever observerd => emerging
  )

  if (is.na(em)) em <- 0

  ptitle <- paste0(spec, "_", lyear, "_", em)
  g <- plot_ts(df = df, ptitle = ptitle)

  df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = em)
  outlist <- list(dt = dt, plot = g, em = df_em)
  return(outlist)
}


## Rule set can be simplyfied by cutting any series with > 4 years with 0 observations
#and only consider the observations after.



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


