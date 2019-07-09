#### Preprocessing and filtering of data

# df <- input data frame from occ-processing
#
# df_pp$ncells <- area of occupancy (number of cells) per species, year
# df_po$occ    <- sum of occurrences per species, year

preproc <- function(df){

  # Total number of species before preprocessing
  print(length(unique(df$taxonKey)))

  # Extract time series with number of cells per year & species
  df_pp <- df %>%
    group_by(taxonKey, year) %>%
    summarise(ncells = n(),
              occ = sum(n))

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

  return(df_pp)
}
