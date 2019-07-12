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
df_bl <- readRDS(file = "./data/cube_belgium_baseline.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

# Do some preprocessing
df_pp <- preproc(df_in)
df_bl <- preprocbl(df_bl)

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

#plot(g[["2115769"]])

# Select number of years for evaluation
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

#saveRDS(gam_result_incr, file = "output/models/gam_result_incr.RDS")

# Plot em status for increasing time series
emGAM <- gam_result_incr %>%  map_dfr(c("em"))

# Summary em status based on evaluation of last year
emGAM %>%
  group_by(em) %>%
  count()


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
  mutate(em = GAM) %>%
  left_join(spec_names, by = "taxonKey")


seltax <- c("1031394")

for (tax in taxl){
  df <- main_incr_df %>%
    filter(taxonKey == tax)
  lyear <- max(df$year)
  ptitle <- paste0(df[[1,1]], "_", df[[1,"spn"]], "_DT_", df[df$year == lyear, "DT"],
                   "_GAM_", df[df$year == lyear, "GAM"])
  plot_incr_em(df = df, ptitle = ptitle, printplot = FALSE, saveplot = TRUE)
}




# een van de drie laatste jaren emering -> emerging
# 2 van de laatste drie jaren emerging -> emerging

# ? Gewogen som?


#### Samenvatting


main_incr_df %>%
  group_by(em) %>%
  filter(year == 2017) %>%
  count()


more_one_cells <- main_incr_df %>%
  ungroup() %>%
  group_by(taxonKey) %>%
  summarise(maxcells = max(ncells)) %>%
  filter(maxcells > 1) %>%
  pull(taxonKey)

# 1690429 - 2015-2017 geven allemaal 0. Komt door breed betrouwbaarheidsinterval
# 2225772 - is traag emerging. met 'tp' smoother em(2013), maar niet met 'ts'
# 2226990 - is emerging. Ook goed resultaat met GAM
# 2227000 - hier ontbreekt echt correctie voor observatieinspanning op klasse niveau
# 2287615 - deze zou als duidelijk emerging moeten klasseren
# 2362868 - is niet emerging. Helemaal gestabiliseerd
# 2379089 - zou emerging moeten zijn. GAM geeft niet emerging
# 2426661 - niet emerging, wel gestabiliseerd. GAM afhankelijk van de smoother
# 1718308 - emerging  = ok




