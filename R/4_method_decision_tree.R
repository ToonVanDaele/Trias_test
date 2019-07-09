# Method decision tree

# The decision wether a species is emerging or not
# is based on multiple simple decision rules.

spDT <- function(df){

  spec <- df[[1,1]]
  lyear <- max(df$year)
  nb <- nrow(df)
  ptitle <- paste0("/dt/", spec, "_DT_",lyear)
  dt <- rep(FALSE, 6)  # vector with results of decision rules tests

  # One value only?   -> appearing species
  if (nb == 1) dt[1] <- TRUE

  # Two values
  if (nb == 2) dt[2] <- TRUE

  # Second value > 0  -> emerging ("3")
  # Second value = 0  -> possibly emerging ("2")
  if (df[[2,"ncells"]] > 0) dt[3] <- TRUE

  # 0 till 5 years back in time -> not emerging
  if (sum(df$ncells[max(nb - 4, 0):nb]) == 0) dt[4] <- TRUE

  # 1 or more with 0 at least 5 years back in time -> (re)appearing
  if (nb > 1 & df$ncells[nb] > 0 & sum(df$ncells[max(nb - 5, 0):(nb - 1)]) == 0) dt[5] <- TRUE

  # Maximum ncells 1 -> sowieso niet emerging
  if (max(df$ncells) <= 1) dt[6] <- TRUE

  # Last value > before last value
  if (nb > 1 & df$ncells[nb] > df$ncells[nb - 1]) dt[7] <- TRUE

  # Combine rules to 0 (not emerging), 1 (possibly emering), 3 (emerging), 4 (re-appearing)
  if (dt[1] == TRUE){
    out <- 3
  }

  # Verder uit te werken...

  out <- 0  # emergence status


  df_em <- tibble(taxonKey = spec, eyear = lyear,  method_em = "DT", em = out)
  outlist <- list(df = df, dt = dt, em = df_em)
  return(outlist)
}



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




# categorie '(re)appearing' species  # nieuwe soort (niet gezien de laatste vijf jaar)
