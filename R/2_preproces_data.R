#### Preprocessing and filtering of Beglium_cube data

# df <- input data frame
#
# df_pp$ncells <- area of occupancy (number of cells) per species, year
# df_pp$obs    <- sum of occurrences per species, year

preproc <- function(df_in, df_bl, spec_names){

  # Add xy
  df_bl$x <- as.integer(substr(df_bl$eea_cell_code, start = 5, stop = 8))
  df_bl$y <- as.integer(substr(df_bl$eea_cell_code, start = 10, stop = 13))

  # Extract time series with number of cells (ncells) and number
  # of occurrences (occ) per year & species
  df_s <- df_in %>%
    filter( year > 1950 & year <= 2017) %>%
    left_join(spec_names %>%
                select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl %>%
                select(year, eea_cell_code, classKey, cobs = n),
              by = c("year", "eea_cell_code", "classKey"))

  df_in2 <- df_s %>%
    group_by(taxonKey, year, classKey) %>%
    summarise(ncells = n(),
              obs = sum(n)
              ncobs = n(),
              cobs = sum(n))

  #### Gaps
  # spread and gather on ncells to add the years with 0 observations
  df_in2 <- df_in2 %>%
    dplyr::select(-obs) %>%
    spread(key = taxonKey, value = ncells) %>%
    gather(key = taxonKey, value = ncells, -year)

  # Remove all zeros before first observation. Define the first year with
  # data (> 0) and drop all records before that year.
  # Join again the occurence data and replace NA's by 0
  df_in2 <- df_in2 %>%
    filter(ncells > 0) %>%
    group_by(taxonKey) %>%
    summarise(minyear = min(year)) %>%
    left_join(df_in2, by = "taxonKey") %>%
    mutate(drop = year < minyear) %>%
    filter(drop == FALSE) %>%
    select(-drop, -minyear) %>%
    left_join(select(df_in2, - ncells), by = c("taxonKey", "year")) %>%
    mutate_all(list(~replace(., is.na(.), 0)))

  # Extract time series with number of occurrences (occ) per year & species
  df_bl2 <- df_bl %>%
    group_by(classKey, year) %>%
    summarise(cobs = sum(n),
              ncobs = n())

  # Join class observations with species observations by classKey
  df_pp <- df_in2 %>%
    left_join(spec_names %>%
                select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl2, by = c("classKey", "year"))

  return(list(df_s, df_pp))
}


## Preprocessing baseline data by year, cell_code, kingdomKey

preproc2 <- function(df_in, df_bl, spec_names) {

  # Filter on period 1950-2017
  df_in2 <- df_in %>%
    filter(year > 1950 & year <= 2017) %>%
    left_join(spec_names %>%
                select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl %>%
                select(year, eea_cell_code, classKey, cobs = n),
              by = c("year", "eea_cell_code", "classKey"))


}


## Preprocessing species data

ppsp <- function(df){

  df_temp <- df %>%
    filter(year > 1950 & year <= 2017)


}


#### Preprocessing
#
# Number of species per year, cellcode

# next step
