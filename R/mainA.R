### Main

# Approaches for the analysis

# dataframe df_pp:
# obs = number of observations
# ncells = number of cells (obs >= 1)
# cobs = number of observations of the species class
# ncobs = number of cells (cobs >= 1)

# A. obs ~ year   (neg. binom - loglink )
# B. ncells ~ year  (neg. binom - loglink)  ncells <<<< tot. number of cells
# C. obs ~ year + cobs  (neg. binom - loglink)
# D. ncells ~ year + cobs  (neg. binom - loglink) ncells <<<< tot. number of cells
# E. ncells ~ year + ncobs  (neg. binom - loglink) ncells <<<< tot. number of cells

# df_s:
#### By year + cellID
# obs = number of observations
# pobs = presence / absence
# cobs = number of observations of species class
# ncobs = number of cells (cobs >= 1)

# E. obs ~ year + cobs  (neg. binom. - log)
# F. pobs ~ year + cobs  (binomial - logit link)
# G. pobs ~ year + ncobs (binomial - logit link)

# Three dataframes are generated
# df_s = obs, cobs by taxonKey + year + cellID
# df_pp = obs, ncells, cobs, ncobs by taxonKey + year
# df_spa = data by year and cellID (presence/absence)

# Init
library(tidyverse)
source(file = "./R/2_preproces_data.R")
source(file = "./R/3_select_species.R")
source(file = "./R/3b_incr_times_series.R")
source(file = "./R/4_method_decision_tree.R")
source(file = "./R/5a_method_piecewise_regression.R")
source(file = "./R/5b_method_INLA_AR.R")
source(file = "./R/5c_method_GAM.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

## Load data
df_in <- readRDS(file = "./data/cube_belgium.RDS")
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")

## Preprocessing
firstyear <- 1980
lastyear <- 2017

df_s <- preproc_s(df_in, df_bl, df_xy, spec_names,
                  firstyear = firstyear, lastyear = lastyear)
df_pp <- preproc_pp(df_s)

# Save preprocessed data
saveRDS(df_s, file = "./data/df_s.RDS")
saveRDS(df_pp, file = "./data/df_pp.RDS")


## Some checks (to be deleted)

df_pp %>%
  filter(taxonKey == "3172100") %>%
  ggplot(aes(x = year, y = obs)) + geom_point() + geom_line()

df_s %>%
  filter(taxonKey == "3172100", year == 2005, obs > 0) %>%
  left_join(df_xy, by = "eea_cell_code") %>%
  ggplot(aes(x = x, y = y)) + geom_point() + coord_fixed()


df_s %>%
  filter(taxonKey == "3172100" & year == 2017) %>%
  plot_map()

# Check of there is at least one observation > 0 for each cell + species

temp <- df_s %>%
  group_by(taxonKey, eea_cell_code) %>%
  summarise(maxobs = max(obs)) %>%
  group_by(taxonKey) %>%
  summarise(min_cell = min(maxobs))
summary(temp)

#listlenght = aantal soorten per jaar en per cel  smoother
# enkel cellen en jaren waar minstens 10 soorten gezien zijn

# gamxy

# x= aantal soorten gezien per cel,  y = kans op voorkomen per cel
