library(tidyverse)
library(INBOtheme)

source(file = "./R/9b_plot_function.R")

df_ts <- read_tsv(file = "./data/df_timeseries.tsv")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")

spec_names$taxonKey <- as.numeric(spec_names$taxonKey)

df_ts <- filter(df_ts, taxonKey == "7501634")

df_ts2 <- df_ts %>%
  left_join(df_xy %>%
              select(eea_cell_code, spa, habitat),
            by = "eea_cell_code") %>%
  select(year, eea_cell_code, obs, pa_obs, natura2000, spa, habitat)

#df_ts_n2k <- filter(df_ts2, natura2000 == TRUE)
df_ts_n2k <- filter(df_ts2, habitat == TRUE)

# Data sets for the lumped data
df_pp <- df_ts2 %>%
  group_by(year) %>%
  summarise(obs = sum(obs),
            ncell = sum(pa_obs)) %>%
  ungroup()

df_pp_n2k <- df_ts_n2k %>%
  group_by(year) %>%
  summarise(obs_n2k = sum(obs),
            ncell_n2k = sum(pa_obs)) %>%
  ungroup()

df_joined <- df_pp %>%
  left_join(df_pp_n2k, by = "year")

write_tsv(df_joined, path = "./rosa_multi.tsv")

df_plot <- df_joined %>%
  select(Jaar = year, All = ncell, Natura2000 = ncell_n2k) %>%
  filter(Jaar > 1990 & Jaar < 2019)

df_plot <- df_plot %>%
  gather(key = locatie, value = occupancy, -Jaar)

df_plot$locatie <- as.factor(df_plot$locatie)
df_plot$locatie <- factor(df_plot$locatie, labels = c("België", "HRL"))

g <- ggplot(df_plot, aes(x = Jaar, y = occupancy, colour = locatie)) +
  geom_point(size = 1.5) +
  geom_line(size = 1) +
  xlim(1990, 2020) +
  theme_set(theme_inbo(transparent = TRUE)) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("Aantal kilometerhokken")

plot(g)

ggsave(filename = "rosa_multi3.eps", plot = g)
ggsave(filename = "rosa_multi3.png", plot = g, dpi = 600)

