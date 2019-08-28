# Functions

# Calculate first and second derivative from smoother 'year'
#
# g gam model
# length length of data frame for prediction
#
# return list with derivatives and emerging evaluation

syear_deriv <- function(g, dfrows){

  deriv1 <- derivatives(g, term = "year", type = "central",
                        order = 1, level = 0.8, n = dfrows, eps = 1e-4)
  deriv2 <- derivatives(g, term = "year", type = "central",
                        order = 2, level = 0.8, n = dfrows, eps = 1e-4)

  em = em_level(deriv1, deriv2)

  return(list(deriv1 = deriv1, deriv2 = deriv2, em = em))

}

# Evaluate emerging based on first and second derivative
#
# df1 result from gratia::derivatives - order 1
# df2 result fomr gratia::derivatives - order 2
#
# Return emerging evaluation

em_level <- function(df1, df2){

  em1 <- df1 %>%
    mutate(em1 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em1)

  em2 <- df2 %>%
    mutate(em2 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em2)

  # Optie om hier continue score van te maken?
  # Combineren van 1ste en 2de afgeleide

  em <- bind_cols(em1, em2) %>%
    mutate(em = case_when(
      em1 == 1   & em2 == 1  ~  "4",
      em1 == 1   & em2 == 0  ~  "3",
      em1 == 1   & em2 == -1 ~  "2",
      em1 == 0   & em2 == 1  ~  "1",
      em1 == 0   & em2 == 0  ~  "0",
      em1 == 0   & em2 == -1  ~ "-1",
      em1 == -1  & em2 == 1  ~  "-2",
      em1 == -1  & em2 == 0  ~  "-3",
      em1 == -1  & em2 == -1  ~ "-4")) %>%
    dplyr::select(-year1)
  return(em)

}

