#### TRIAS - preprocessing and filtering of Belgium_cube data


## Preprocessing data by ~ year + taxonKey + cellID

# df_in raw data (from getdata.R)
# df_bl baseline data (species class) (from getdata.R)
# spec_names list with species names (from getdata.R)
# firstyear first year of data to filter (1950)
# lastyear last year of data to filter (2017)
#
# return df_pp time series with number of cells (ncells) and number of
#          observations (obs) ~ year + species + cellID

preproc_s <- function(df_in, df_bl, df_xy, spec_names,
                      firstyear = 1980, lastyear = 2017) {

  # Select observations within the requested period (years) and
  # within Belgium (rectangle). Records outside belgium should not be here.
  df_in2 <- df_in %>%
    left_join(df_xy,
              by = "eea_cell_code") %>%
    filter(year > firstyear & year <= lastyear) %>%
    filter(x >= 3788 & x <= 4065 & y >= 2930 & y <= 3201) %>%
    dplyr::select(taxonKey, year, eea_cell_code, obs = n)

  # Sum of invasive species observations grouped by cell_code + year + class
  df_ispec <- df_in2 %>%
    left_join(spec_names %>%
                dplyr::select(taxonKey, classKey), by = c("taxonKey")) %>%
    group_by(classKey, eea_cell_code, year) %>%
    summarise(ispec = sum(obs))

  # Substract the invasive species observations from the baseline data
  df_bl2 <- df_bl %>%
    left_join(df_ispec, by = c("classKey", "eea_cell_code", "year")) %>%
    replace_na(list(ispec = 0)) %>%
    mutate(cobs = n - ispec) %>%
    mutate(new_cobs = ifelse(cobs < 0, 0, cobs)) %>%
    select(year, eea_cell_code, classKey, cobs = new_cobs)

  df_bl2$cobs <- as.integer(df_bl2$cobs)

  # The second last line corrects an error which results in negative cobs.
  # i.e. more observed invasive species than total observed species.
  # This is not possible. Still need to check were the error comes from.
  # It concerns only a small number of records (~200).
  # There's no impact on the results

  # Identify for each species all cells with at least one observation
  df_cc <- df_in2 %>%
    group_by(taxonKey) %>%
    distinct(eea_cell_code)

  # Identify for each species the first year with an observation
  df_by <- df_in2 %>%
    group_by(taxonKey) %>%
    summarize(begin_year = min(year))

  # Join begin_year to unique location + species combination
  df_cc2 <- df_cc %>%
    left_join(df_by, by = "taxonKey")

  # Create a time series from begin_year till lastyear for each cell and species.
  # Only empty records here. No data yet.

  # Remove unused data frames to release some memory
  remove(df_in, df_cc, df_bl)

  epg <- function(eea_cell_code, taxonKey, begin_year, lastyear) {

    expand.grid(eea_cell_code = eea_cell_code,
                taxonKey = taxonKey,
                year = seq(from = begin_year, to = lastyear),
                stringsAsFactors = FALSE)

  }

  df_s <- pmap_dfr(df_cc2, .f = epg, lastyear = lastyear)

  # Join observations of the class ~ cell_code + year  and
  # join the observations of invasive species ~ cell_code + year + taxonKey
  df_s <- df_s %>%
    left_join(spec_names %>%
                dplyr::select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl2 %>%
                dplyr::select(year, eea_cell_code, classKey, cobs),
              by = c("year", "eea_cell_code", "classKey")) %>%
    left_join(df_in2,
              by = c("taxonKey", "year", "eea_cell_code")) %>%
    replace_na(list(cobs = 0, obs = 0))

  #Set obs & cobs to integer
  df_s$obs <- as.integer(df_s$obs)
  df_s$cobs <- as.integer(df_s$cobs)

  # Add columns for presence/absence data
  df_s <- df_s %>%
    mutate(pa_cobs = ifelse(cobs > 0, 1, 0),
           pa_obs = ifelse(obs > 0, 1, 0))

  df_s$pa_obs <- as.integer(df_s$pa_obs)
  df_s$pa_cobs <- as.integer(df_s$pa_cobs)
  return(df_s)

}


#### Preprocessing data summarised over all cells

preproc_pp <- function(df_s){

  # Total observations and class observations
  df_pp <- df_s %>%
    group_by(taxonKey, year) %>%
    summarise(obs = sum(obs),
              cobs = sum(cobs),
              ncell = sum(pa_obs),
              ncobs = sum(pa_cobs)) %>%
    ungroup()

  return(df_pp)
}

