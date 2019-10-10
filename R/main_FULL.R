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

####################################################################
# 1. Preprocessing
# spatial data frame - full
df_s <- preproc_s(df_in, df_bl, df_xy, spec_names,
                  firstyear = firstyear, lastyear = lastyear)
saveRDS(df_s, file = "./data/df_s.RDS")

# spatial data frame - n2k only
df_s_n2k <- preproc_s_N2k(df_s)
saveRDS(df_s_n2k, file = "./data/df_s_n2k.RDS")

# lumped data frame - full
df_pp <- preproc_pp(df_s)
saveRDS(df_pp, file = "./data/df_pp.RDS")

# lumped data_frame - n2k only
df_pp_n2k <- preproc_pp(df_s_n2k)
saveRDS(df_pp_n2k, file = "./data/df_pp_n2k.RDS")

# Aggregate spatial data frames to 5x5km cells (instead of 1x1km)
df_s5 <- aggr_1to5(df_s)
saveRDS(df_s5, file = "./data/df_s5.RDS")

df_s5_n2k <- aggr_1to5(df_s_n2k)
saveRDS(df_s5_n2k, file = "./data/df_s5_n2k.RDS")

########################################################################
# 2. modelling

## Load preprocessed data (if needed)
#df_s <- readRDS(file = "./data/df_s.RDS")   # Don't load 1x1km
#df_s_n2k <- readRDS(file = "./data/df_s_n2k.RDS")   # Don't load 1x1km

df_pp <- readRDS(file = "./data/df_pp.RDS")
df_pp_n2k <- readRDS(file = "./data/df_pp_n2k.RDS")

df_s5 <- readRDS(file = "./data/df_s5.RDS")
df_s5_n2k <- readRDS(file = "./data/df_s5_n2k.RDS")

## Selection of species for testing
#df_s <- selspec(df = df_s, specs = spec_names$taxonKey[1:20])
#df_s_n2k <- selspec(df = df_s_n2k, specs = spec_names$taxonKey[1:20])
df_pp <- selspec(df = df_pp, specs = spec_names$taxonKey[1:100])
df_s5 <- selspec(df = df_s5, specs = spec_names$taxonKey[1:100])
df_pp_n2k <- selspec(df = df_pp_n2k, specs = spec_names$taxonKey[1:100])
df_s5_n2k <- selspec(df = df_s5_n2k, specs = spec_names$taxonKey[1:100])

## Apply each method on all species

## Full data set
apply_method(df_pp, "spDT")   # Decision tree (no statistics)

apply_method(df_pp, "spGAM_lcount")  # GAM occurences on lumped data
apply_method(df_s5, "spGAM_count_ns") # GAM occurences + cobs  5x5km
#apply_method(df_s5, "spGAM_count")    # GAM occurences + cobs + s(x,y)  5x5km

apply_method(df_pp, "spGAM_lpa")    # GAM occupancy on lumped data
apply_method(df_s5, "spGAM_pa_ns")  # GAM occupancy + cobs 5x5km
#apply_method(df_s5, "spGAM_pa")     # GAM occupancy + cobs + s(x,y) 5x5km

## Data in N2000 only
apply_method(df_pp_n2k, "spDT", n2k = TRUE)   # Decision tree (no statistics) n2k

apply_method(df_pp_n2k, "spGAM_lcount", n2k = TRUE)  # GAM occurences on lumped data
apply_method(df_s5_n2k, "spGAM_count_ns", n2k = TRUE) # GAM occurences + cobs  5x5km
#apply_method(df_s5_n2k, "spGAM_count", n2k = TRUE)    # GAM occurences + cobs + s(x,y)  5x5km

apply_method(df_pp_n2k, "spGAM_lpa", n2k = TRUE)    # GAM occupancy on lumped data
apply_method(df_s5_n2k, "spGAM_pa_ns", n2k = TRUE)  # GAM occupancy + cobs 5x5km
#apply_method(df_s5_n2k, "spGAM_pa", n2k = TRUE)     # GAM occupancy + cobs + s(x,y) 5x5km


# 3. Process output

# List out output to be processed
# full data
out_list_full <- list("result_spDT", "result_spGAM_lcount", "result_spGAM_count_ns",
                      "result_spGAM_lpa", "result_spGAM_pa_ns")

# Extract 'emerging' information from each method and join
result_full <- out_list_full %>%
  map_dfr(get_em, path = "./output/")

saveRDS(result_full, file = "./output/result_full.RDS")

# natura2000 only
out_list_n2k <- as.list(paste0(out_list_full, "_n2k"))

result_n2k <- out_list_n2k %>%
  map_dfr(get_em, path = "./output/")

saveRDS(result_n2k, file = "./output/result_n2k.RDS")

result_n2k$method_em <- paste0(result_n2k$method_em, "_n2k")


em <- rbind(result_full, result_n2k) %>%
  spread(key = method_em, value = em)

em <- rbind(emDT, emGAM_lcount, emGAM_count_ns, emGAM_lpa, emGAM_pa_ns) %>%
  spread(key = method_em, value = em)

# Correction of some error because an earlier model run is use
em <- em %>%
  rename(GAM_lcount = GAM)

em <- em %>%
  left_join(spec_names %>%
              select(taxonKey, spn),
            by = "taxonKey") %>%
  select('taxonKey', 'spn', 'eyear',
         'DT_n2k', 'GAM_lcount_n2k', 'GAM_count_ns_n2k', 'GAM_lpa_n2k', 'GAM_pa_ns_n2k',
         'DT', 'GAM_lcount', 'GAM_count_ns', 'GAM_lpa', 'GAM_pa_ns')

# should be corrected too
result_models_df <- em

# Get one value per indicator (if GAM_*: NA -> DT_*)
result_indicator_df <-
  result_models_df %>%
  group_by(taxonKey, spn, eyear) %>%
  summarize(occ = ifelse(is.na(GAM_count_ns),
                         ifelse(is.na(GAM_lcount),
                                DT,
                                GAM_lcount),
                         GAM_count_ns),
            aoo = ifelse(is.na(GAM_pa_ns),
                         ifelse(is.na(GAM_lpa),
                                      DT,
                                      GAM_lpa),
                                GAM_pa_ns),
            occ_n2k = ifelse(is.na(GAM_count_ns_n2k),
                         ifelse(is.na(GAM_lcount_n2k),
                                DT_n2k,
                                GAM_lcount_n2k),
                         GAM_count_ns_n2k),
            aoo_n2k = ifelse(is.na(GAM_pa_ns_n2k),
                         ifelse(is.na(GAM_lpa_n2k),
                                DT_n2k,
                                GAM_lpa_n2k),
                         GAM_pa_ns_n2k)
  ) %>%
  ungroup()

# Ranking
ranking_df_occ <-
  result_indicator_df %>%
  select(taxonKey, spn, eyear, occ) %>%
  spread(key = eyear, value = occ, sep = "occ_")

ranking_df_aoo <-
  result_indicator_df %>%
  select(taxonKey, eyear, aoo) %>%
  spread(key = eyear, value = aoo, sep = "aoo_")

ranking_df_occ_n2k <-
  result_indicator_df %>%
  select(taxonKey, eyear, occ_n2k) %>%
  spread(key = eyear, value = occ_n2k, sep = "occ_n2k_")

ranking_df_aoo_n2k <-
  result_indicator_df %>%
  select(taxonKey, eyear, aoo_n2k) %>%
  spread(key = eyear, value = aoo_n2k, sep = "aoo_n2k_")

ranking_df <- ranking_df_occ %>%
  left_join(ranking_df_aoo, by = "taxonKey") %>%
  left_join(ranking_df_occ_n2k, by = "taxonKey") %>%
  left_join(ranking_df_aoo_n2k, by = "taxonKey") %>%
  rename_at(vars(starts_with("eyear")), ~str_remove(., pattern = "eyear")) %>%
  group_by(taxonKey, spn) %>%
  arrange(
    desc(aoo_n2k_2017),
    desc(occ_n2k_2017),
    desc(aoo_n2k_2016),
    desc(occ_n2k_2016),
    desc(aoo_n2k_2015),
    desc(occ_n2k_2015),
    desc(aoo_2017),
    desc(occ_2017),
    desc(aoo_2016),
    desc(occ_2016),
    desc(aoo_2015),
    desc(occ_2015),
    taxonKey) %>%
  select(taxonKey, spn,
         aoo_n2k_2017, occ_n2k_2017,
         aoo_n2k_2016, occ_n2k_2016,
         aoo_n2k_2015, occ_n2k_2015,
         aoo_2017, occ_2017,
         aoo_2016, occ_2016,
         aoo_2015, occ_2015)


ranking_df

saveRDS(ranking_df, file = "./output/ranking_df.RDS")
write.csv(ranking_df, file = "./output/ranking_df.csv")

## Decide the emerging status based on the multiple methods




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



df_s %>%
  filter(taxonKey == "3189866" & year == "2017") %>%
  left_join(df_xy, by = "eea_cell_code") %>%
  ggplot(aes(x =x, y = y)) + geom_point()
