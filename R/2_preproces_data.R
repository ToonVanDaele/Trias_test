#### Preprocessing and filtering of data

# df <- input data frame from occ-processing
#
# df_pp$ncells <- area of occupancy (number of cells) per species, year
# df_po$occ    <- sum of occurrences per species, year

preproc <- function(df){

  # Total number of species before preprocessing
  print(length(unique(df$taxonKey)))

  # Extract time series with number of cells (ncells) and number
  # of occurrences (occ) per year & species
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
  # spread and gather on ncells to add the years with 0 observations
  df_temp <- df_pp %>%
    dplyr::select(-occ) %>%
    spread(key = taxonKey, value = ncells) %>%
    gather(key = taxonKey, value = ncells, -year)

  # Remove all zeros before first observation. Define the first year with
  # data (> 0) and drop all records before that year.
  # Join again the occurence data and replace NA's by 0
  df_temp <- df_temp %>%
    filter(ncells > 0) %>%
    group_by(taxonKey) %>%
    summarise(minyear = min(year)) %>%
    left_join(df_temp, by = "taxonKey") %>%
    mutate(drop = year < minyear) %>%
    filter(drop == FALSE) %>%
    select(-drop, -minyear) %>%
    left_join(select(df_pp, - ncells), by = c("taxonKey", "year")) %>%
    mutate_all(list(~replace(., is.na(.), 0)))

  return(df_temp)
}
