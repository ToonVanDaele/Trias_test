### Main  prop cells

# Init
library(tidyverse)
source(file = "./R/2_preproces_data.R")
source(file = "./R/3_select_species.R")
source(file = "./R/3b_incr_times_series.R")
source(file = "./R/5c_method_GAM.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

# Get data
df_in <- readRDS(file = "./data/cube_belgium.RDS")
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

# Preprocessing

# Remove record with incorrect eaa cell code
df_bl <- filter(df_bl, !(eea_cell_code == "1kmE-2305N659"))

# add x, y
df_bl$x <- as.integer(substr(df_bl$eea_cell_code, start = 5, stop = 8))
df_bl$y <- as.integer(substr(df_bl$eea_cell_code, start = 10, stop = 13))



# For each kingdom key a map with the number of observations by cell

df_bl %>%
  filter(kingdomKey == 6, year == 2017) %>%
  ggplot(aes(x = x, y = y, colour = n)) + geom_point()

df_bl %>%
  filter(year == 2017, n <= 50) %>%
  ggplot(aes(x = n)) + geom_histogram() +
  facet_wrap(~kingdomKey)



df_bl <- ppbl(df_bl)
head(df_bl)

df_p <- df_in %>%
  left_join(kingdomkey)
  full_join(df_bl,
            by = c("year", "eea_cell_code"))



df %>%
  filter(year == 2017 & n < 25) %>%
  group_by(n) %>%
  count() %>%
  ggplot(aes(x = n, y = nn)) + geom_line()

# Select species (for testing)
df_sp <- selspec(df_pp) %>%
  group_by(taxonKey)

# Get a vector with the group names
taxl <- df_sp %>%
  group_keys() %>%
  pull(taxonKey)

# plot time series
g <- df_sp %>%
  group_split() %>%
  map(.f = plot_ts, printplot = TRUE, saveplot = TRUE) %>%
  set_names(taxl)





#listlenght = aanal soorten per jaar en per cel  smoother

#enkel cellen en jaren waar minstens 10 soorten gezien zijn

# gamxy

# smoother -> listlength
# x= aantal soorten gezien per cel,  y = kans op voorkomen per cel

