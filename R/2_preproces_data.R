#### Preprocessing and filtering of data

preproc <- function(df){

  # Total number of species before preprocessing
  length(unique(df$taxonKey))

  # Extract time series with number of cells per year & species
  df_pp <- df %>%
    group_by(taxonKey, year) %>%
    summarise(ncells = n())

  #### Period
  # Minimum and maximum year
  # - year > 1950
  # - year <= 2017
  df_pp <- filter(df_pp, year > 1950 & year <= 2017)

  #### Gaps
  # Gaps of less than 5 years (no obs during < 5 years)
  # Between first and last observation replace NA with 0
  df_pp <- df_pp %>%
    spread(key = taxonKey, value = ncells) %>%
    gather(key = taxonKey, value = ncells, -year) %>%
    mutate_all(list(~replace(., is.na(.), 0)))

  # Remove zeros before first observation
  df_pp <- df_pp %>%
    filter(ncells > 0) %>%
    group_by(taxonKey) %>%
    summarise(minyear = min(year)) %>%
    left_join(df_pp, by = "taxonKey") %>%
    mutate(drop = year < minyear) %>%
    filter(drop == FALSE) %>%
    select(-drop, -minyear)


  #### Minimum time series
  # At least 5 consecutive observations till 2017
  # Species that are not observed since a long time are not considered
  # just for testing purpose. Requires fine tuning

  # specfilter <- df_pp %>%
  #   summarise(maxyear = max(year)) %>%
  #   filter(maxyear >= 2017) %>%
  #   .$taxonKey
  #
  #
  # df <- df %>%
  #   filter(taxonKey %in% specfilter)
  #
  #
  #
  #
  #
  #
  #
  #
  # # When there is a gap of > 5 years (no observations during > 5 years),
  # # the data before the gap is not considered.
  #
  #
  #
  #
  # # Total number of species after preprocessing
  # length(unique(df$taxonKey))

  return(df_pp)
}
