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


##### Run the different methods

## Decision tree (no statistics)
df_dt_out <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spDT)


## Piecewise regresion
# Complete time series
PR_result <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spPR)

# Plot results and save
for (i in PR_result) {
  df <- i$df
  df_n <- i$df_n
  ptitle <- paste0("PR/", df[[1,1]], "_", max(df$year))
  if (!is.null(df_n)) {
    plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle, printplot = TRUE)
  }else{
    plot_ts(df = df, ptitle = ptitle, printplot = TRUE)
  }
}

# Incrementing time series
df_outPR <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spPR)

# Plot results and save
for (i in df_outPR) {
  df <- i$df
  df_n <- i$df_n
  ptitle <- paste0("PR/", df[[1,1]], "_", max(df$year))
  if (!is.null(df_n)) {
    plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle, printplot = FALSE)
  }else{
    plot_ts(df = df, ptitle = ptitle, printplot = FALSE)
  }
}


## GAM

# Complete time series only
GAM_result <- df_sp %>%
  group_split(taxonKey) %>%
  map(.f = spGAM)

# Plot results and save
for (i in GAM_result){
  df <- i$df
  df_n <- i$df_n
  df_n$ucl[df_n$ucl > 10000] <- 10000
  df_n$lcl[df_n$lcl > 10000] <- 10000
  ptitle <- paste0("GAM/", df[[1,1]], "_", max(df$year))
  plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle, printplot = TRUE)
}

# Incrementing time series
df_outGAM <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spGAM)

for (i in df_outGAM){
  df <- i$df
  df_n <- i$df_n
  df_n$ucl[df_n$ucl > 10000] <- 10000
  df_n$lcl[df_n$lcl > 10000] <- 10000
  ptitle <- paste0("GAM/", df[[1,1]], "_", max(df$year))
  plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle, printplot = TRUE)
}


# ? Why are the times series plots only up to 2016 and not 2017??

# Retrieve and combine the em result for each method, taxonKey & eyear combination
emGAM <- df_outGAM %>%  map_dfr(c("em")) %>%
emPR <-  df_outPR %>% map_dfr(c("em"))

df_main <- df_sp %>%
  left_join(emGAM %>% rename(gam = em), by = c("taxonKey" = "taxonKey", "year" = "eyear")) %>%
  left_join(emPR %>%
              rename(pr = em) %>%
              dplyr::select(-method_em),  by = c("taxonKey" = "taxonKey", "year" = "eyear"))

df_main <- df_main %>%
  gather(key = method_em, value = em, -taxonKey, -year, -ncells)

# Plot (only GAM for the moment)
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


tpr <- spPR(df = df)
tgam <- spGAM(df = df)

ptitle <- paste0(df[[1,1]])
plot_ribbon_em(df_n = tpr[["df_n"]], df = tpr[["df"]], ptitle = ptitle, printplot = TRUE)
plot_ribbon_em(df_n = tgam[["df_n"]], df = tgam[["df"]], ptitle = ptitle, printplot = TRUE)

# Output

#def emerging
#- absoluut aantal
#- relatief aantal (eerste en tweede afgeleide)

# Korte termijn eerste afgeleide,
# langere termijn PR lineaire regressie coefficient verschillend van 0

# Voorbeeld enkele tijdreeksen
# Uitwerking in GAM -> + eerste afgeleide voor emerging

# Uitwerking PR -> voorbeeld emerging / not emering + probleem

# Uitwerking in INLA -> RW2 + eerste afgeleide als emerging



