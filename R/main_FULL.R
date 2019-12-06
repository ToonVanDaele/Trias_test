### Main full analysis

# This full analysis (after 1_getdata.R)

# Some overall parameters
firstyear <- 1980
lastyear <- 2018
nbyear <- 3   # Number of last years to be evaluated

# Load libraries & source function
library(tidyverse)
library(gridExtra)
source(file = "./R/4_method_decision_tree.R")
source(file = "./R/5c_method_GAM_short.R")
source(file = "./R/5d_method_GAM_pa.R")
source(file = "./R/5e_method_GAM_count.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

## Load data
df_ts <- read_tsv(file = "./data/df_timeseries.tsv")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")

# Filter some species - problem with taxonKey
df_ts <- filter(df_ts, !taxonKey %in% c("10173593", "10172912", "10661581"))

# Check for species in spec_names and add new species names if necessary
spec_names <- add_spec(df_ts, spec_names)

# Filter last year
df_ts <- filter(df_ts, year <= lastyear)

# Selection of species for testing
df_ts <- filter(df_ts, taxonKey %in% unique(df_ts$taxonKey)[1:30])

# Aggregate spatial data frames to 5x5km cells (instead of 1x1km)
df_s5 <- aggr_1to5(df_ts)

## Create data set for natura2000
# Identify for each species the first year with an observation
df_begin_n2k <- df_ts %>%
  filter(natura2000 == TRUE & obs > 0) %>%
  group_by(taxonKey) %>%
  summarize(begin_year = min(year))

df_ts_n2k <- df_ts %>%
  filter(natura2000 == TRUE) %>%
  left_join(df_begin_n2k, by = "taxonKey") %>%
  filter(year >= begin_year) %>%
  select(-begin_year)

df_s5_n2k <- aggr_1to5(df_ts_n2k)

# Data sets for the lumped data
df_pp <- df_ts %>%
  group_by(taxonKey, year) %>%
  summarise(obs = sum(obs),
            cobs = sum(native_obs),
            ncell = sum(pa_obs),
            ncobs = sum(pa_native_obs)) %>%
  ungroup()

df_pp_n2k <- df_ts_n2k %>%
  group_by(taxonKey, year) %>%
  summarise(obs = sum(obs),
            cobs = sum(native_obs),
            ncell = sum(pa_obs),
            ncobs = sum(pa_native_obs)) %>%
  ungroup()

########################################################################
# 2. modelling
#
# All locations

#apply_method(df_s, "spGAM_count")    # GAM occurences + cobs + s(x,y)  5x5km
#apply_method(df_s, "spGAM_pa")     # GAM occupancy + cobs + s(x,y) 5x5km

# apply_method(df_ts, "spGAM_count_ns") # GAM occurences + native_obs  1x1km
# apply_method(df_ts, "spGAM_pa_ns")  # GAM occupancy + native_obs 1x1km

apply_method(df_pp, "spGAM_lcount_cobs")  # GAM number of cells with native_obs
apply_method(df_pp, "spGAM_lcount")        # GAM occurences on lumped data

apply_method(df_pp, "spGAM_lpa_cobs")    # GAM occupancy on lumped data
apply_method(df_pp, "spGAM_lpa")    # GAM occupancy on lumped data

apply_method(df_pp, "spDT")                # Decision tree (no statistics)


## In N2000 only

#apply_method(df_s5_n2k, "spGAM_count", n2k = TRUE)    # GAM occurences + cobs + s(x,y)  5x5km
#apply_method(df_s5_n2k, "spGAM_pa", n2k = TRUE)     # GAM occupancy + cobs + s(x,y) 5x5km

# apply_method(df_ts_n2k, "spGAM_count_ns", n2k = TRUE) # GAM occurences + cobs  1x1km
# apply_method(df_ts_n2k, "spGAM_pa_ns", n2k = TRUE)  # GAM occupancy + cobs 1x1km

apply_method(df_pp_n2k, "spGAM_lcount_cobs", n2k = TRUE)  # GAM occurences on lumped data
apply_method(df_pp_n2k, "spGAM_lcount", n2k = TRUE)  # GAM occurences on lumped data

apply_method(df_pp_n2k, "spGAM_lpa_cobs", n2k = TRUE)    # GAM occupancy on lumped data
apply_method(df_pp_n2k, "spGAM_lpa", n2k = TRUE)    # GAM occupancy on lumped data

apply_method(df_pp_n2k, "spDT", n2k = TRUE)   # Decision tree (no statistics) n2k

#######################################################################""
# 3. Process output

# full data
#out_list_full <- list("result_spDT", "result_spGAM_lcount", "result_spGAM_count_ns",
#                      "result_spGAM_lpa", "result_spGAM_pa_ns")

out_list_full <- list("result_spGAM_lpa_cobs", "result_spGAM_lcount_cobs",
                      "result_spGAM_lpa", "result_spGAM_lcount",
                      "result_spDT")

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

df_em <- rbind(result_full, result_n2k) %>%
  select(-lcl) %>%
  spread(key = method_em, value = em)

spec_names$taxonKey <- as.numeric(spec_names$taxonKey)

# selection_list <- c('taxonKey', 'spn', 'eyear',
# 'DT_n2k', 'GAM_lcount_n2k', 'GAM_count_ns_n2k', 'GAM_lpa_n2k', 'GAM_pa_ns_n2k',
# 'DT', 'GAM_lcount', 'GAM_count_ns', 'GAM_lpa', 'GAM_pa_ns')

selection_list <- c('taxonKey', 'spn', 'eyear',
                    'DT_n2k', 'GAM_lcount_cobs_n2k', 'GAM_lcount_n2k',
                    'GAM_lpa_cobs_n2k', 'GAM_lpa_n2k',
                    'DT', 'GAM_lcount_cobs', 'GAM_lcount',
                    'GAM_lpa_cobs', 'GAM_lpa')

df_em <- df_em %>%
  left_join(spec_names %>%
              select(taxonKey, spn),
            by = "taxonKey") %>%
  select(selection_list)

# Add lower confidence level value for GAM_lpa_cobs_n2k
df_em <- df_em %>%
  left_join(result_n2k %>%
              filter(method_em == "GAM_lpa_cobs_n2k") %>%
              select(taxonKey, eyear, lcl),
            by = c("taxonKey", "eyear"))

# Get one value per indicator (if GAM_*: NA -> DT_*)
result_indicator_df <-
  df_em %>%
  group_by(taxonKey, spn, eyear) %>%
  summarize(occ = ifelse("GAM_count_ns" %in% names(df_em) && !is.na(GAM_count_ns),
                         GAM_count_ns,
                         ifelse(is.na(GAM_lcount),
                                DT,
                                GAM_lcount)),
            aoo = ifelse("GAM_pa_ns" %in% names(df_em) && !is.na(GAM_pa_ns),
                         GAM_pa_ns,
                         ifelse(is.na(GAM_lpa),
                                      DT,
                                      GAM_lpa)),
            occ_n2k = ifelse("GAM_count_ns_n2k" %in% names(df_em) && !is.na(GAM_count_ns_n2k),
                             GAM_count_ns_n2k,
                             ifelse(is.na(GAM_lcount_n2k),
                                DT_n2k,
                                GAM_lcount_n2k)),
            aoo_n2k = ifelse("GAM_pa_ns_n2k" %in% names(df_em) && !is.na(GAM_pa_ns_n2k),
                             GAM_pa_ns_n2k,
                             ifelse(is.na(GAM_lpa_n2k),
                                DT_n2k,
                                GAM_lpa_n2k)),
            lcl = lcl) %>%
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

df_mean_lcl <-
  result_indicator_df %>%
  group_by(taxonKey) %>%
  summarise(mean_lcl = round(mean(lcl), 3))

ranking_df <- ranking_df_occ %>%
  left_join(ranking_df_aoo, by = "taxonKey") %>%
  left_join(ranking_df_occ_n2k, by = "taxonKey") %>%
  left_join(ranking_df_aoo_n2k, by = "taxonKey") %>%
  left_join(df_mean_lcl, by = "taxonKey") %>%
  rename_at(vars(starts_with("eyear")), ~str_remove(., pattern = "eyear")) %>%
  group_by(taxonKey, spn) %>%
  arrange(
    desc(aoo_n2k_2018),
    desc(occ_n2k_2018),
    desc(aoo_n2k_2017),
    desc(occ_n2k_2017),
    desc(aoo_n2k_2016),
    desc(occ_n2k_2016),
    desc(aoo_2018),
    desc(occ_2018),
    desc(aoo_2017),
    desc(occ_2017),
    desc(aoo_2016),
    desc(occ_2016),
    desc(mean_lcl),
    taxonKey) %>%
  select(taxonKey, spn,
         mean_lcl,
         aoo_n2k_2018, occ_n2k_2018,
         aoo_n2k_2017, occ_n2k_2017,
         aoo_n2k_2016, occ_n2k_2016,
         aoo_2018, occ_2018,
         aoo_2017, occ_2017,
         aoo_2016, occ_2016)

ranking_df

saveRDS(ranking_df, file = "./output/ranking_df.RDS")
write.csv(ranking_df, file = "./output/ranking_df.csv")


## Generate plots
result <- readRDS("./output/result_spGAM_lcount_cobs.RDS")

dir.create(paste0("./output/plots/"), showWarnings = FALSE)

map(.x = result, .f = output_multiple_plots)




# ## Generate plots based on "aoo_n2k" (result_indicator_df)
# df_plot <- df_pp %>%
#   left_join(result_indicator_df %>%
#               select(taxonKey, eyear, aoo_n2k),
#             by = c("taxonKey" = "taxonKey", "year" = "eyear")) %>%
#   rename(em = aoo_n2k)
#
# df_plot <- group_by(df_plot, taxonKey)
# taxl <- group_keys(df_plot) %>% pull()
#
# plot_em <- df_plot %>%
#   group_split() %>%
#   map(plot_incr_em, saveplot = TRUE) %>%
#   set_names(taxl)
#
# saveRDS(plot_em, file = "./output/plot_em.RDS")
#
