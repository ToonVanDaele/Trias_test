#### TIRAS - preprocessing and filtering of Belgium_cube data


## Preprocessing data by ~ year + taxonKey + cellID

# df_in raw data (from getdata.R)
# df_bl baseline data (species class) (from getdata.R)
# spec_names list with species names (from getdata.R)
# firstyear first year of data to filter (1950)
# lastyear last year of data to filter (2017)
#
# return df_pp time series with number of cells (ncells) and number of
#          observations (obs) ~ year + species + cellID

preproc_s <- function(df_in, df_bl, spec_names, df_xy,
                      firstyear = 1950, lastyear = 2017) {

  # Select observations within the requested period (years) and
  # within Belgium (rectangle). Records outside belgium should not be here.
  df_in2 <- df_in %>%
    left_join(df_xy,
              by = "eea_cell_code") %>%
    filter(year > firstyear & year <= lastyear) %>%
    filter(x >= 3768 & x <= 4079 & y >= 2926 & y <= 3236) %>%
    dplyr::select(taxonKey, year, eea_cell_code, obs = n)

  # Calculate the total number of invasive species observations
  # grouped by cell_code + year + class
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
    mutate(new_cobs = if_else(cobs < 0, 0, cobs)) %>%
    select(year, eea_cell_code, classKey, cobs = new_cobs)

  df_bl2$cobs <- as.integer(df_bl2$cobs)

  # The second last line corrects an error which results in negative cobs.
  # i.e. more observed invasive species than total observed species.
  # This is not possible. Still need to check were the error comes from.
  # However, it concerns only a small number of records (~200).
  # There's no impact on the results

  # Identify unique locations with at least one observation for each species
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

  # Joint observations of the class ~ cell_code + year
  # and the observations of invasive species ~ cell_code + year + taxonKey
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

}


# Preprocessing data summarised over all cells

preproc_pp <- function(df_s){

  df_1 <- df_s %>%
    group_by(taxonKey, year) %>%
    summarise(obs = sum(obs),
              cobs = sum(cobs))

  df_2 <- df_s %>%
    filter(obs > 0) %>%
    group_by(taxonKey, year) %>%
    summarise(ncells = n_distinct(eea_cell_code))

  df_3 <- df_s %>%
    filter(cobs > 0) %>%
    group_by(taxonKey, year) %>%
    summarise(ncobs = n_distinct(eea_cell_code))

  df_pp <- df_1 %>%
    left_join(df_2,
              by = c("taxonKey", "year")) %>%
    left_join(df_3,
               by = c("taxonKey", "year")) %>%
    replace_na(list(ncells = 0, ncobs = 0)) %>%
    ungroup()

  return(df_pp)
}


## Preprocessing for presence/absence (not used)

# preproc_pa <- function(df_in, df_bl, spec_names, firstyear, lastyear) {
#
#   spec <- "3172100"
#   specc <- spec_names[spec_names$taxonKey == spec, "classKey"]
#   spn <- spec_names[spec_names$taxonKey == spec, "spn"]
#
#   df_temp <- df_in %>%
#     filter(taxonKey == spec) %>%
#     filter(year > firstyear & year <= lastyear) %>%
#     rename(obs = n)
#
#   fyear <- min(df_temp$year)
#
#   df_temp <- df_temp %>%
#     full_join(df_bl %>%
#                 filter(classKey == specc & year >= fyear),
#               by = c("year", "eea_cell_code")) %>%
#     rename(cobs = n) %>%
#     filter(!is.na(obs) | (is.na(obs) & cobs > 20))
#
#   df_temp$taxonKey <- spec
#   df_temp$spn <- spn
#   df_temp$obs[is.na(df_temp$obs)] <- 0
#   df_temp$pa <- ifelse(df_temp$obs > 0, 1, 0)
#
#   # Add xy
#   df_temp$x <- as.integer(substr(df_temp$eea_cell_code, start = 5, stop = 8))
#   df_temp$y <- as.integer(substr(df_temp$eea_cell_code, start = 10, stop = 13))
#
#   df_temp %>%
#     filter(year == 2010) %>%
#     ggplot(aes(x = x, y = y, colour = obs)) + geom_point()
#
#   df_bl %>%
#     filter(classKey == 220) %>%
#     ggplot(aes(x = year, y = n)) + geom_point() + geom_line()
#
#
# }
#
#

# Preprocessing data by ~ year + taxonKey
#
# df_in raw data (from getdata.R)
# df_bl baseline data (species class) (from getdata.R)
# spec_names list with species names (from getdata.R)
# firstyear first year of data to filter
# lastyear last year of data to filter
#
# return df_pp time series with number of cells (ncells) and number of
#          observations (obs) ~ year + species.

# preproc <- function(df_in, df_bl, spec_names, firstyear, lastyear){
#
#   df_in2 <- df_in %>%
#     filter( year > firstyear & year <= lastyear) %>%
#     group_by(taxonKey, year) %>%
#     summarise(ncells = n(),
#               obs = sum(n))
#
#   # Years without observation of the species have no records
#   # create data frame with all taxonkeys and all possible years
#   df_f <- expand.grid(taxonKey = unique(df_in$taxonKey),
#                       year = seq(from = firstyear, to = lastyear),
#                       stringsAsFactors = FALSE)
#
#   df_in3 <- df_in2 %>%
#     full_join(df_f, by = c("taxonKey", "year")) %>%
#     replace_na(list(ncells = 0, obs = 0)) %>%
#     arrange(taxonKey, year)
#
#   # Remove all zeros before the first observation.
#   # The first year is the first year with obs > 0.
#   # All zeroes before that year are dropped.
#   df_in4 <- df_in3 %>%
#     filter(obs > 0) %>%
#     group_by(taxonKey) %>%
#     summarise(minyear = min(year)) %>%
#     left_join(df_in3, by = "taxonKey") %>%
#     mutate(drop = year < minyear) %>%
#     filter(drop == FALSE) %>%
#     dplyr::select(-drop, -minyear) %>%
#     arrange(taxonKey, year)
#
#   # Add class occurrences (cobs) & cell class occurrences (ncobs)
#   # Warnings about introduction of NA's is ok, they're filtered later.
#
#   df_bl2 <- df_bl %>%
#     mutate(x = as.integer(substr(eea_cell_code, start = 5, stop = 8)),
#            y = as.integer(substr(eea_cell_code, start = 10, stop = 13))) %>%
#     filter(x >= 3768 & x <= 4079 & y >= 2926 & y <= 3236) %>%
#     filter(year > firstyear & year <= lastyear) %>%
#     group_by(classKey, year) %>%
#     summarise(cobs = sum(n),
#               ncobs = n())
#
#   # Add classKey
#   df_in5 <- df_in4 %>%
#     left_join(spec_names %>%
#                 dplyr::select(taxonKey, classKey),
#               by = "taxonKey")
#
#   # The number of observations of invasive species in each class
#   df_bl3 <- df_in5 %>%
#     group_by(classKey, year) %>%
#     summarise(invobs = sum(obs)) %>%
#     left_join(df_bl2,
#               by = c("classKey", "year")) %>%
#     replace_na(list(invobs = 0))
#
#   # Join class observations with species observations (via classKey)
#   df_pp <- df_in5 %>%
#     left_join(df_bl3, by = c("classKey", "year")) %>%
#     mutate(cobs = cobs - invobs)
#
#   return(df_pp)
# }
