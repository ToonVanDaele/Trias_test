# Increasing time series
#
# backward argument indicates how many timesteps
#
# The time series in df is duplicated with increasing length

# This allows to verify how the algoritm behaves with increasing length of time series.

# Value
# eyear = evaluation year (the last year of each time series)

filter_year <- function(eyear, tK, df){
  df_out <- df %>%
    filter(taxonKey == tK & year <= eyear)
  if (nrow(df_out) > 0) df_out$eyear <- eyear
  return(df_out)
}


incrts <- function(df, backward = 4){

  lyear <- max(df$year)

  eval_year <- c(seq(from = lyear - backward, to = lyear, by = 1))

  df_temp <- expand.grid(taxonKey = unique(df$taxonKey),
                         eval_year = eval_year)

  t <- map2_dfr(.x = df_temp$eval_year,
                .y = df_temp$taxonKey,
                .f = filter_year,
                df)
}


## New set of functions

filter_year2 <- function(eyear, df){
  df_out <- df %>%
    filter(year <= eyear)
  if (nrow(df_out) > 0) df_out$eyear <- eyear
  return(df_out)
}


# # Create incresing time series
# dfincr <- function(df, eval_year, meth) {
#
#   t <- map_dfr(.x = eval_year,
#                .f = filter_year2,
#                df)
#
#   out <- t %>%
#     group_by(eyear) %>%
#     group_split() %>%
#     map(.f = get(meth)) %>%
#     set_names(unique(t$eyear))
#
#   return(out)
# }
