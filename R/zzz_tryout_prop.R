### Main for proportional data

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
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

# Do some preprocessing
df_pp <- preproc(df_in)
df_bl <- preprocbl(df_bl)
names(df_bl)[3] <- "occ_bl"

# Select species (for testing)
df_sp <- selspec(df_pp) %>%
  group_by(taxonKey)

# Get kingdomKey and baseline occurence prop = occ / occ_bl
df_sp <- df_sp %>%
  left_join(spec_names, by = "taxonKey") %>%
  left_join(df_bl, by = c("kingdomKey", "year")) %>%
  mutate(prop = occ / occ_bl)

# Get a vector with the group names
taxl <- df_sp %>%
  group_keys() %>%
  pull(taxonKey)

# plot time series proportions
g <- df_sp %>%
  group_split() %>%
  map(.f = plotprop, printplot = TRUE, saveplot = TRUE) %>%
  set_names(taxl)


# Simsons paradox


