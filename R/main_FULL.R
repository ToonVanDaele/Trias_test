### Main full analysis

# This script does the full analysis

# Some overall parameters
firstyear <- 1980
lastyear <- 2017
yb <- 2   # Number of last years to be evaluated
eval_year <- seq(from = lastyear - yb, to = lastyear, by = 1) # years to be evaluated

# Load libraries & source function
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

## Load raw data
# .RDS file are generated with 1_getdata.R
df_in <- readRDS(file = "./data/cube_belgium.RDS")
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")

## Preprocessing
df_s <- preproc_s(df_in, df_bl, df_xy, spec_names,
                  firstyear = firstyear, lastyear = lastyear)
df_pp <- preproc_pp(df_s)

saveRDS(df_s, file = "./data/df_s.RDS")
saveRDS(df_pp, file = "./data/df_pp.RDS")

## Read preprocessed data
#df_s <- readRDS(file = "./data/df_s.RDS")
df_pp <- readRDS(file = "./data/df_pp.RDS")

## Selection of species for testing
df_sp <- selspec(df_pp = df_pp)

# Set group & retrieve group names (i.e. vector with taxonKeys)
df_sp <- df_sp %>% group_by(taxonKey)
taxl <- df_sp %>% group_keys() %>% pull(taxonKey)

## Apply each method on all species

# Decision tree (no statistics)
result_dt <- df_sp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spDT") %>%
  set_names(taxl)
saveRDS(result_dt, file = "./output/result_dt.RDS")

# GAM
result_gam <- df_sp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM") %>%
  set_names(taxl)
saveRDS(result_gam, file = "./output/result_gam.RDS")

# INLA
result_inla <- df_sp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spINLA") %>%
  set_names(taxl)
saveRDS(result_inla, file = "./output/result_inla.RDS")

## Extract 'emerging' information from each method and join
result_dt <- readRDS(file = "./output/result_dt.RDS")
result_gam <- readRDS(file = "./output/result_gam.RDS")
result_inla <- readRDS(file = "./output/result_inla.RDS")

emDT <- result_dt %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM <- result_gam %>%
  map_dfr(~ .x %>% map_dfr("em"))

emINLA <- result_inla %>%
  map_dfr(~ .x %>% map_dfr("em"))

em <- rbind(emDT, emGAM, emINLA) %>%
  spread(key = method_em, value = em)

## Decide the emerging status based on the multiple methods

# For testing use GAM unless GAM = NA, then use DT
em <- em %>%
  mutate(emtot = ifelse(!is.na(GAM), GAM, DT))

## Generate plots
df_em <- df_sp %>%
  left_join(em, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_em <- df_em %>%
  group_split() %>%
  map(plot_incr_em, saveplot = TRUE) %>%
  set_names(taxl)

saveRDS(plot_em, file = "./output/plot_em.RDS")


## Summary statistics
em %>%
  mutate(emtot = as.numeric(emtot)) %>%
  group_by(emtot, eyear) %>%
  summarise(em_nb = n()) %>%
  group_by(eyear) %>%
  spread(key = eyear, value = em_nb) %>%
  arrange(emtot)
