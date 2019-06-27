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
temp <- df_out %>%  map_dfr(c("em"))

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

# Complete time series only
df_pr_out <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spGAM)

# Incrementing time series
df_outGAM <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spGAM)


# ? Why are the times series only up to 2016 and not 2017??


# Retrieve the em output for each taxonKey & eyear combination
temp <- df_outGAM %>%  map_dfr(c("em"))

df_resultGAM <- df_sp %>%
  left_join(temp, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_emGAM <- df_resultGAM %>%
  group_by(taxonKey) %>%
  group_map(~tibble(plots = list(plot_incr_em5(.x, spec = .y, printplot = TRUE))))



t <- spGAM(filter(df_sp, taxonKey == "8193935"))

draw(t$model)
appraise(t$model)
summary(t$model)

draw(t$deriv1)
draw(t$deriv2)

tail(t$deriv1)
tail(t$deriv2)
t$em


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
