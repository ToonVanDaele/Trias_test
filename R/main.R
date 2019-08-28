### Main2

# Approaches for the analysis

# obs = number of observations
# ncells = number of cells (obs >= 1)
# cobs = number of observations of the species class
# ncobs = number of cells (cobs >= 1)

# A. obs ~ year   (neg. binom - loglink )
# B. ncells ~ year  (neg. binom - loglink)  ncells <<<< tot. number of cells
# C. obs ~ year + cobs  (neg. binom - loglink)
# D. ncells ~ year + cobs  (neg. binom - loglink) ncells <<<< tot. number of cells
# E. ncells ~ year + ncobs  (neg. binom - loglink) ncells <<<< tot. number of cells

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
source(file = "./R/5c_method_GAM.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

## Load data
df_in <- readRDS(file = "./data/cube_belgium.RDS")
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

## Preprocessing
fyear <- 1950
lyear <- 2017

df_pp <- preproc(df_in, df_bl, spec_names, firstyear = fyear, lastyear = lyear)
df_s <- preproc_s(df_in, df_bl, spec_names, firstyear = fyear, lastyear = lyear)
#df_spa <- preproc_pa(df_in, df_bl, spec_names, firstyear = fyear, lastyear = lyear)

# Save preprocessed data
saveRDS(df_pp, file = "./data/df_pp.RDS")
saveRDS(df_s, file = "./data/df_s.RDS")

df_pp %>%
  filter(taxonKey == "3172100") %>%
  ggplot(aes(x = year, y = obs)) + geom_point() + geom_line()

df_pp %>%
  filter(taxonKey == "3172100") %>%
  ggplot(aes(x = year, y = cobs)) + geom_point() + geom_line()


df_pp %>%
  filter(taxonKey == "3172100") %>%
  ggplot(aes(x = cobs, y = obs)) + geom_point()


df_pp %>%
  filter(taxonKey == "3172100") %>%
  ggplot(aes(x = ncobs, y = ncells)) + geom_point() + geom_smooth()




#listlenght = aanal soorten per jaar en per cel  smoother

#enkel cellen en jaren waar minstens 10 soorten gezien zijn

# gamxy

# smoother -> listlength
# x= aantal soorten gezien per cel,  y = kans op voorkomen per cel
