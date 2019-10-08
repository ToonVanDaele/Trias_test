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
source(file = "./R/5c_method_GAM_short.R")
source(file = "./R/5d_method_GAM_pa.R")
source(file = "./R/5e_method_GAM_count.R")
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
df_s <- readRDS(file = "./data/df_s.RDS")
df_pp <- readRDS(file = "./data/df_pp.RDS")

## Selection of species for testing
df_s <- selspec(df = df_s, specs = spec_names$taxonKey[1:20])
df_pp <- selspec(df = df_pp, specs = spec_names$taxonKey[1:20])

# Temporary workaround upscaling to 5x5km cells
df_s <- df_s %>%
  left_join(df_xy %>%
              select(eea_cell_code, x5, y5, cell_code5), by = "eea_cell_code")

df_s5 <- df_s %>%
  group_by(taxonKey, year, classKey, cell_code5) %>%
  summarise(x = first(x5),
            y = first(y5),
            cobs = sum(cobs),
            obs = sum(obs),
            pa_cobs = max(pa_cobs),
            pa_obs = max(pa_obs))

# Set group & retrieve group names (i.e. vector with taxonKeys)
df_pp <- df_pp %>% group_by(taxonKey)
df_s <- df_s %>% group_by(taxonKey)
df_s5 <- df_s5 %>% group_by(taxonKey)

taxl <- df_pp %>% group_keys() %>% pull(taxonKey)

## Apply each method on all species

# Decision tree (no statistics)
result_dt <- df_pp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spDT") %>%
  set_names(taxl)
saveRDS(result_dt, file = "./output/result_dt.RDS")

# GAM_lcount
result_gam_lcount <- df_pp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_lcount") %>%
  set_names(taxl)
saveRDS(result_gam_lcount, file = "./output/result_gam_lcount.RDS")

# GAM_lpa
result_gam_lpa <- df_pp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_lpa") %>%
  set_names(taxl)
saveRDS(result_gam_lpa, file = "./output/result_gam_lpa.RDS")

# GAM_count_ns
result_gam_count_ns <- df_s5 %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_count_ns") %>%
  set_names(taxl)
saveRDS(result_gam_count_ns, file = "./output/result_gam_count_ns.RDS")

# GAM_count
result_gam_count <- df_s5 %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_count") %>%
  set_names(taxl)

# GAM_pa_ns
result_gam_pa_ns <- df_s5 %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_pa_ns") %>%
  set_names(taxl)
saveRDS(result_gam_pa_ns, file = "./output/result_gam_pa_ns.RDS")

# GAM_pa
result_gam_pa <- df_s5 %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spGAM_pa") %>%
  set_names(taxl)
saveRDS(result_gam_pa, file = "./output/result_gam_pa.RDS")






# INLA
result_inla <- df_pp %>%
  group_split() %>%
  map(.f = dfincr, eval_year, "spINLA") %>%
  set_names(taxl)
saveRDS(result_inla, file = "./output/result_inla.RDS")

## Extract 'emerging' information from each method and join
result_dt <- readRDS(file = "./output/result_dt.RDS")
result_gam_lcount <- readRDS(file = "./output/result_gam_lcount.RDS")
result_gam_lpa <- readRDS(file = "./output/result_gam_lpa.RDS")
result_gam_count <- readRDS(file = "./output/result_gam_count.RDS")
result_gam_pa <- readRDS(file = "./output/result_gam_pa.RDS")
result_gam_count_ns <- readRDS(file = "./output/result_gam_count_ns.RDS")
result_gam_pa_ns <- readRDS(file = "./output/result_gam_pa_ns.RDS")

result_inla <- readRDS(file = "./output/result_inla.RDS")

emDT <- result_dt %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_lcount <- result_gam_lcount %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_count <- result_gam_count %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_count_ns <- result_gam_count_ns %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_lpa <- result_gam_lpa %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_pa <- result_gam_pa %>%
  map_dfr(~ .x %>% map_dfr("em"))

emGAM_pa_ns <- result_gam_pa_ns %>%
  map_dfr(~ .x %>% map_dfr("em"))

emINLA <- result_inla %>%
  map_dfr(~ .x %>% map_dfr("em"))

em <- rbind(emDT, emGAM_lcount, emGAM_count_ns, emGAM_count, emGAM_lpa, emGAM_pa_ns, emGAM_pa) %>%
  spread(key = method_em, value = em)

em <- em %>%
  left_join(spec_names %>%
              select(taxonKey, spn),
            by = "taxonKey") %>%
  select('taxonKey', 'spn', 'eyear', 'DT', 'GAM_lcount', 'GAM_count_ns', 'GAM_count',
         'GAM_lpa', 'GAM_pa_ns', 'GAM_pa')


## Decide the emerging status based on the multiple methods

# For testing use GAM unless GAM = NA, then use DT
em <- em %>%
  mutate(emtot = ifelse(!is.na(GAM), GAM, DT))

## Generate plots
df_em <- df_pp %>%
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
