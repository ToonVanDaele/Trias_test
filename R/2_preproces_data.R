#### Preprocessing and filtering of Beglium_cube data

# df_in, df_bl, spec_names
#
# df_pp <- time series with number of cells (ncells) and number of
#          occurrences (occ) ~ year + species
# df_pp$ncells <- area of occupancy (number of cells) per species, year
# df_pp$obs    <- sum of occurrences per species, year

preproc <- function(df_in, df_bl, spec_names, firstyear, lastyear){

  df_in2 <- df_in %>%
    filter( year > firstyear & year <= lastyear) %>%
    group_by(taxonKey, year) %>%
    summarise(ncells = n(),
              obs = sum(n))

  # Years without observation of the species have no records
  # create data frame with all taxonkeys and all possible years
  df_f <- expand.grid(taxonKey = unique(df_in$taxonKey),
                      year = seq(from = firstyear, to = lastyear),
                      stringsAsFactors = FALSE)

  df_in3 <- df_in2 %>%
    full_join(df_f, by = c("taxonKey", "year")) %>%
    replace_na(list(ncells = 0, obs = 0)) %>%
    arrange(taxonKey, year)

  # Remove all zeros before the first observation.
  # The first year is the first year with obs > 0.
  # All zeroes before that year are dropped.
  df_in4 <- df_in3 %>%
    filter(obs > 0) %>%
    group_by(taxonKey) %>%
    summarise(minyear = min(year)) %>%
    left_join(df_in3, by = "taxonKey") %>%
    mutate(drop = year < minyear) %>%
    filter(drop == FALSE) %>%
    dplyr::select(-drop, -minyear) %>%
    arrange(taxonKey, year)

  # Add class occurrences (cobs) & cell class occurrences (ncobs)
  # Warnings about introduction of NA's is ok, they're filtered later.
  df_bl2 <- df_bl %>%
    mutate(x = as.integer(substr(eea_cell_code, start = 5, stop = 8)),
           y = as.integer(substr(eea_cell_code, start = 10, stop = 13))) %>%
    filter(x >= 3768 & x <= 4079 & y >= 2926 & y <= 3236) %>%
    filter(year > firstyear & year <= lastyear) %>%
    group_by(classKey, year) %>%
    summarise(cobs = sum(n),
              ncobs = n())

  # Join class observations with species observations (via classKey)
  df_pp <- df_in4 %>%
    left_join(spec_names %>%
                select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl2, by = c("classKey", "year"))

  return(df_pp)
}


## Preprocessing baseline data by year, cell_code, kingdomKey

preproc_s <- function(df_in, df_bl, spec_names, firstyear, lastyear) {

  df_in2 <- df_in %>%
    filter(year > firstyear & year <= lastyear) %>%
    mutate(x = as.integer(substr(eea_cell_code, start = 5, stop = 8)),
           y = as.integer(substr(eea_cell_code, start = 10, stop = 13))) %>%
    left_join(spec_names %>%
                select(taxonKey, classKey),
              by = "taxonKey") %>%
    left_join(df_bl %>%
                select(year, eea_cell_code, classKey, cobs = n),
              by = c("year", "eea_cell_code", "classKey"))
}


## Preprocessing for presence/absence

preproc_pa <- function(df_in, df_bl, spec_names, firstyear, lastyear) {

  spec <- "3172100"
  specc <- spec_names[spec_names$taxonKey == spec, "classKey"]
  spn <- spec_names[spec_names$taxonKey == spec, "spn"]

  df_temp <- df_in %>%
    filter(taxonKey == spec) %>%
    filter(year > firstyear & year <= lastyear) %>%
    rename(obs = n)

  fyear <- min(df_temp$year)

  df_temp <- df_temp %>%
    full_join(df_bl %>%
                filter(classKey == specc & year >= fyear),
              by = c("year", "eea_cell_code")) %>%
    rename(cobs = n) %>%
    filter(!is.na(obs) | (is.na(obs) & cobs > 20))

  df_temp$taxonKey <- spec
  df_temp$spn <- spn
  df_temp$obs[is.na(df_temp$obs)] <- 0
  df_temp$pa <- ifelse(df_temp$obs > 0, 1, 0)

  # Add xy
  df_temp$x <- as.integer(substr(df_temp$eea_cell_code, start = 5, stop = 8))
  df_temp$y <- as.integer(substr(df_temp$eea_cell_code, start = 10, stop = 13))

  df_temp %>%
    filter(year == 2010) %>%
    ggplot(aes(x = x, y = y, colour = obs)) + geom_point()

  df_bl %>%
    filter(classKey == 220) %>%
    ggplot(aes(x = year, y = n)) + geom_point() + geom_line()


}


