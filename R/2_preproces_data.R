#### Preprocessing and filtering of Beglium_cube data

# df <- input data frame
#
# df_pp$ncells <- area of occupancy (number of cells) per species, year
# df_po$occ    <- sum of occurrences per species, year

preproc <- function(df){

  # Extract time series with number of cells (ncells) and number
  # of occurrences (occ) per year & species
  df_pp <- df %>%
    group_by(taxonKey, year) %>%
    summarise(ncells = n(),
              occ = sum(n))

  #### Period 1950 < year <= 2017
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

#### Preprocessing and filtering of Beglium_cube_baseline data

# df <- input data frame from getdata

preprocbl <- function(df){

  # Extract time series with number of occurrences (occ) per year & species
  df_temp <- df %>%
    group_by(kingdomKey, year) %>%
    summarise(occ = sum(n))

  #### Period between 1950 and 2017
  df_temp <- filter(df_temp, year > 1950 & year <= 2017)

  #### Gaps
  # spread and gather on occ to add the years with 0 observations
  df_temp2 <- df_temp %>%
    spread(key = kingdomKey, value = occ) %>%
    gather(key = kingdomKey, value = occ, -year)

  # Remove all NA's before first observation. Define the first year with
  # data (> 0) and drop all records before that year.
  # Join again the occurence data and replace NA's by 0
  df_temp2 <- df_temp2 %>%
    filter(occ > 0) %>%
    group_by(kingdomKey) %>%
    summarise(minyear = min(year)) %>%
    left_join(df_temp, by = "kingdomKey") %>%
    mutate(drop = year < minyear) %>%
    filter(drop == FALSE) %>%
    select(-drop, -minyear) %>%
    mutate_all(list(~replace(., is.na(.), 0)))

  return(df_temp2)
}


## Preprocessing baseline data by year, cell_code, kingdomKey

ppbl <- function(df) {

  # Filter on period 1950-2017
  df_temp <- df %>%
    filter(year > 1950 & year <= 2017) %>%
    filter(n >= 1)
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
