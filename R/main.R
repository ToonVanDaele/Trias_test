### Main2

# 5 approaches for the analysis

# obs = number of observations
# ncells = number of cells (obs > 1)
# cobs = number of observations of the species class
# ncobs = number of cells (cobs > 1)

# A. obs ~ year   (neg. binom - loglink )
# B. ncells ~ year  (neg. binom - loglink)  ncells <<<< tot. number of cells
# C. obs ~ year + cobs  (neg. binom - loglink)
# D. ncells ~ year + cobs  (neg. binom - loglink) ncells <<<< tot. number of cells
# E. ncells ~ year + ncobs  (neg. binom - loglink) ncells <<<< tot. number of cells

#### By year + cellID
# obs = number of observations
# pobs = presence / absence
# cobs = number of observations of species class
# ncobs = number of cells (cobs > 1)

# E. obs ~ year + cobs  (neg. binom. - log)
# F. pobs ~ year + cobs  (binomial - logit link)
# G. pobs ~ year + ncobs (binomial - logit link)

# Three dataframes are generated:
# df_s = data by year and cellID
# df_pp = data sumarized by year
# df_spa = data by year and cellID (presence/absence)

# Init
library(tidyverse)
source(file = "./R/2_preproces_data.R")
source(file = "./R/3_select_species.R")
source(file = "./R/3b_incr_times_series.R")
source(file = "./R/5c_method_GAM.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

## Load data
df_in <- readRDS(file = "./data/cube_belgium.RDS")
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

## Preprocessing
temp <- preproc(df_in, df_bl, spec_names)
df_s <- temp[[1]]
df_pp <- temp[[2]]
remove(temp)

# presence absence to be added


# baseline data
# add x, y
df_bl$x <- as.integer(substr(df_bl$eea_cell_code, start = 5, stop = 8))
df_bl$y <- as.integer(substr(df_bl$eea_cell_code, start = 10, stop = 13))

# filter within xy belgian limits
# year between 1950 and 2017
# minimum 10 observations per class/year/cell
df_bl <- df_bl %>%
  filter(year > 1950 & year < 2017) %>%
  filter(x >= 3768 & x <= 4079 & y >= 2926 & y <= 3236)

head(df_bl)
head(df_in)
head(spec_names)

# Join class and species name to df_in
df_in <- df_in %>%
  rename(obs = n) %>%
  left_join(select(spec_names,
                   -kingdomKey), by = "taxonKey")

head(df_in)

# For each species join the baseline with the proper classKey
# and create a data frame with presences and absences

spec <- "3171948"
specc <- spec_names[spec_names$taxonKey == spec, "classKey"]
spn <- spec_names[spec_names$taxonKey == spec, "spn"]

temp <- df_in %>%
  filter(taxonKey == spec) %>%
  full_join(df_bl %>%
              filter(classKey == specc),
            by = c("year", "eea_cell_code", "classKey")) %>%
  rename(cobs = n) %>%
  filter(!is.na(obs) | (is.na(obs) & cobs > 20))

temp$taxonKey <- spec
temp$spn <- spn
temp$obs[is.na(temp$obs)] <- 0

temp %>%
  ggplot(aes(x = cobs, y = obs)) + geom_point()





#listlenght = aanal soorten per jaar en per cel  smoother

#enkel cellen en jaren waar minstens 10 soorten gezien zijn

# gamxy

# smoother -> listlength
# x= aantal soorten gezien per cel,  y = kans op voorkomen per cel
