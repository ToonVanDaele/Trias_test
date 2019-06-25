### Main

# Init
library(tidyverse)

# Get data
df_in <- readRDS(file = "./data/cube_belgium.RDS")

# Do some preprocessing
df_pp <- preproc(df_in)

# Select species
# (Mainly for testing)
df_sp <- selspec(df_pp)
length(unique(df_sp$taxonKey))

# plot time series
g <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = plot_ts, printplot = FALSE)

# return van plot_ts aanpassen. Geeft nu NULL terug

# Create increasing time series (for testing)
df_spe <- incrts(df_sp)

# Method: decision tree
# For each unique species & evaluation year in df_sp
df_out <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spDT)

# Select only the em output for each taxonKey & eyear
temp <- df_out %>%  map_dfr(c(2))

df_result <- df_sp %>%
  left_join(temp, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_em <- df_result %>%
  group_by(taxonKey) %>%
  group_map(~tibble(plots = list(plot_incr_em(.x, spec = .y, printplot = TRUE))))


# Run methods
df_dt_out <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spDT)

df_pr_out <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spPR)


# Using GAM

# Complete time series
df_pr_out <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spGAM)

# Incrementing time series
df_outGAM <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spGAM)

# Select only the em output for each taxonKey & eyear
temp <- df_outGAM %>%  map_dfr(c(2))

df_resultGAM <- df_sp %>%
  left_join(temp, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_emGAM <- df_resultGAM %>%
  group_by(taxonKey) %>%
  group_map(~tibble(plots = list(plot_incr_em5(.x, spec = .y, printplot = TRUE))))







#spPR(df_sp)
#spINLA(df_sp)
#spGAM(df_sp)

# Output

#def emerging
#- absoluut aantal
#- relatief aantal (eerste afgeleide)

# Korte termijn eerste afgeleide,
# langere termijn PR lineaire regressie coefficient verschillend van 0

# Voorbeeld enkele tijdreeksen
# Uitwerking in GAM -> + eerste afgeleide voor emerging

# Uitwerking PR -> voorbeeld emerging / not emering + probleem

# Uitwerking in INLA -> RW2 + eerste afgeleide als emerging
