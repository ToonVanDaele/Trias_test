library(tidyverse)
library(INBOtheme)

source(file = "./R/9b_plot_function.R")

df_ts <- read_tsv(file = "./data/df_timeseries.tsv")
spec_names <- readRDS(file = "./data/spec_names.RDS")
spec_names$taxonKey <- as.numeric(spec_names$taxonKey)

df_ts <- filter(df_ts, taxonKey == "7501634")

df_ts_n2k <- filter(df_ts, natura2000 == TRUE)

# Data sets for the lumped data
df_pp <- df_ts %>%
  group_by(taxonKey, year) %>%
  summarise(obs = sum(obs),
            ncell = sum(pa_obs)) %>%
  ungroup()

df_pp_n2k <- df_ts_n2k %>%
  group_by(taxonKey, year) %>%
  summarise(obs_n2k = sum(obs),
            ncell_n2k = sum(pa_obs)) %>%
  ungroup()

df_nl <- data.frame(taxonKey = 7501634,
                    soortNL = "Rosa multiflora",
                    stringsAsFactors = FALSE)

df_joined <- df_pp %>%
  left_join(df_pp_n2k, by = c("taxonKey", "year"))

write_tsv(df_joined, path = "./rosa_multi.tsv")

df_plot <- df_joined %>%
  left_join(spec_names %>%
              select(taxonKey, spn), by = c("taxonKey")) %>%
  left_join(df_nl, by = "taxonKey") %>%
  select(Jaar = year, All = ncell, Natura2000 = ncell_n2k, soortnaam = spn) %>%
  filter(Jaar > 1990 & Jaar < 2019)

df_plot <- df_plot %>%
  gather(key = locatie, value = occupancy, -Jaar, -soortnaam)

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

ggsave(filename = "rosa_multi2.eps", plot = g)
ggsave(filename = "rosa_multi2.png", plot = g, dpi = 300)

