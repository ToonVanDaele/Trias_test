### Main

# Init
library(tidyverse)
source(file = "./R/2_preproces_data.R")
source(file = "./R/3_select_species.R")
source(file = "./R/3b_incr_times_series.R")
source(file = "./R/4_method_decision_tree.R")
source(file = "./R/5a_method_piecewise_regression.R")
source(file = "./R/5c_method_GAM.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

# Get data
df_in <- readRDS(file = "./data/cube_belgium.RDS")

# Do some preprocessing
df_pp <- preproc(df_in)

# Select species
# (Mainly for testing)
df_sp <- selspec(df_pp) %>%
  group_by(taxonKey)

# Get a vector with the group names
taxl <- df_sp %>%
  group_keys() %>%
  pull(taxonKey)

# plot time series
g <- df_sp %>%
  group_split() %>%
  map(.f = plot_ts, printplot = FALSE, saveplot = FALSE) %>%
  set_names(taxl)

#plot(g[[1]])  #plot(g[["2115769"]])

# Create increasing time series
df_spe <- incrts(df_sp, backward = 2) %>%
  group_by(taxonKey, eyear)

# Retrieve a vector with grouping key (used to name the output list later)
taxl_incr <- df_spe %>%
  group_keys() %>%
  mutate(key = paste(taxonKey, eyear, sep = "_")) %>%
  pull(key)


##### Run the different methods

## Decision tree (no statistics)

# on full time series
dt_result <- df_sp %>%
  group_split() %>%
  map(.f = spDT)

# on increasing time series (unique species & evaluation year)
dt_result_incr <- df_spe %>%
  group_split() %>%
  map(.f = spDT)

# Select only the 'em' output for each taxonKey & eyear combination
dt_em <- df_sp %>%
  left_join(dt_result_incr %>%
              map_dfr(c("em")),
            by = c("taxonKey" = "taxonKey", "year" = "eyear"))

dt_em_plot <- dt_em %>%
  group_split() %>%
  map(.f = plot_incr_em, saveplot = FALSE) %>%
  set_names(taxl)

#plot(dt_em_plot[["8542672"]])


## GAM

# GAM on full time series
gam_result <- df_sp %>%
  group_split() %>%
  map(.f = spGAM, saveplot = FALSE) %>%
  set_names(taxl)

# Plot results
for (i in gam_result) {
  if (!is.null(i$plot)) plot(i$plot)
}

# Incrementing time series
gam_result_incr <- df_spe %>%
  group_split() %>%
  map(.f = spGAM, saveplot = FALSE) %>%
  set_names(taxl_incr)

# Plot em status for increasing time series
emGAM <- gam_result_incr %>%  map_dfr(c("em"))
gam_main <- df_sp %>%
  left_join(emGAM, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_emGAM <- gam_main %>%
  group_split() %>%
  map(plot_incr_em5, saveplot = FALSE)



### Join outputs from different methods

# Results from full time series

emDT <- dt_result %>% map_dfr("em")
emGAM <- gam_result %>% map_dfr("em")

main_f <- rbind(emDT, emGAM) %>%
  spread(key = method_em, value = em)


# Results from increasing time series

emDTincr <- dt_result_incr %>% map_dfr("em")
emGAMincr <- gam_result_incr %>% map_dfr("em")
main_incr <- rbind(emDTincr, emGAMincr) %>%
  spread(key = method_em, value = em)

# plot results increasing time series

main_incr_df <- df_sp %>%
  left_join(main_incr, by = c("taxonKey" = "taxonKey", "year" = "eyear")) %>%
  mutate(em = GAM)

main_incr_p <- main_incr_df %>%
  group_split() %>%
  map(plot_incr_em, saveplot = FALSE) %>%
  set_names(taxl)

plot(main_incr_p[["1722299"]])

