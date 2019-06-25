# Increasing time series

# For each time series in the data frame
# duplicate the time series with increasing length
# This is for testing purpose.
# Verify how the algoritm behaves with increasing length of time series.

# Create time series for each evaluation year

incrts <- function(df){

  #eval_year <- 2017
  eval_year <- c(seq(from = 1900, to = 2017, by = 1))

  df_temp <- expand.grid(taxonKey = unique(df$taxonKey),
                         eval_year = eval_year)


  filter_year <- function(eyear, tK, df){
    df_out <- df %>%
      filter(taxonKey == tK & year <= eyear)
    if (nrow(df_out) > 0) {
      df_out <- cbind(df_out, eyear)
    }
    return(df_out)
  }


  t <- map2_dfr(.x = df_temp$eval_year,
                .y = df_temp$taxonKey,
                .f = filter_year,
                df)
}
