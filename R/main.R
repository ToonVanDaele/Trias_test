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
  group_split(taxonKey) %>%
  map(.f = plot_ts, printplot = FALSE, saveplot = FALSE)

#plot(g[[1]])

# Create increasing time series (for testing)
df_spe <- incrts(df_sp, backward = 5) %>%
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
  group_split(taxonKey) %>%
  map(.f = spDT)

# on increasing time series (unique species & evaluation year)
dt_result_incr <- df_spe %>%
  group_split(taxonKey, eyear) %>%
  map(.f = spDT)

# Select only the em output for each taxonKey & eyear
dt_em <- dt_result_incr %>%  map_dfr(c("em"))

dt_result <- df_sp %>%
  left_join(dt_em, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

dt_em_plot <- dt_result %>%
  group_split(taxonKey) %>%
  map(.f = plot_incr_em)

#plot(dt_em_plot[[1]])


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
# get a vector with the group keys
df_spe <- group_by(df_spe, taxonKey, eyear)

gam_result_incr <- df_spe %>%
  group_split() %>%
  map(.f = spGAM, saveplot = FALSE) %>%
  set_names(taxl_incr)

# Plot em status for increasing time series
emGAM <- gam_result_incr %>%  map_dfr(c("em"))
gam_main <- df_sp %>%
  left_join(emGAM, by = c("taxonKey" = "taxonKey", "year" = "eyear"))

plot_emGAM <- gam_main %>%
  group_split(taxonKey) %>%
  map(plot_incr_em5)



### Join outputs from different methods

# Retrieve and combine the em result for each method, taxonKey & eyear combination
emDT <- df_outDT %>%  map_dfr(c("em"))
emGAM <- gam_result_incr %>%  map_dfr(c("em"))
emPR <-  df_outPR %>% map_dfr(c("em"))

df_main <- df_sp %>%
  left_join(emGAM %>%
              rename(gam = em) %>%
              dplyr::select(-method_em), by = c("taxonKey" = "taxonKey", "year" = "eyear")) %>%
  left_join(emPR %>%
              rename(pr = em) %>%
              dplyr::select(-method_em),  by = c("taxonKey" = "taxonKey", "year" = "eyear"))

df_main <- df_main %>%
  gather(key = method_em, value = em, -taxonKey, -year, -ncells)

df_main







df <- filter(df_sp, taxonKey == "7068379")

t <- spGAM(df)

draw(t$model)
appraise(t$model)
summary(t$model)

draw(t$deriv1)
draw(t$deriv2)

tail(t$deriv1)
tail(t$deriv2)
t$em


tpr <- spPR(df = df)

df <- filter(df_pp, taxonKey == "2882849")
plot_ts(df, printplot = TRUE)
tgam <- spGAM(df = df)


tpr <- spPR(df = df)
tpr$df_n

ptitle <- paste0(df[[1,1]])
plot_ribbon_em(df_n = tgam[["df_n"]], df = tgam[["df"]], ptitle = ptitle, printplot = TRUE)

plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle, printplot = TRUE)

plot_ribbon_em(df_n = tpr[["df_n"]], df = tpr[["df"]], ptitle = ptitle, printplot = TRUE)
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


a <- list(x = c(4,5,6),
          y = c(7,8,9))

b <- list(p = c(1,2,3),
          q = c(4,5,6))


t <- list(a, b)

tt <- set_names(x = t, nm = c("first", "second"))


tt["first"]
tt[["first"]]$x

taxl <- df_sp %>%
  group_keys() %>%
  pull(taxonKey)
