library(tidyverse)
library(INBOtheme)

source(file = "./R/9b_plot_function.R")


df_pp <- readRDS(file = "./data/df_pp.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_nl <- data.frame(taxonKey = c("3084015", "3084022"),
                    soortNL = c("Westerse karmozijnbes", "Oosterse karmozijnbes"),
                    stringsAsFactors = FALSE)

df_plot <- df_pp %>%
  filter(taxonKey %in% c("3084015", "3084022")) %>%
  left_join(spec_names %>%
              select(taxonKey, spn), by = c("taxonKey")) %>%
  left_join(df_nl, by = "taxonKey") %>%
  select(Jaar = year, occupancy = ncells, soortNL, soortnaam = spn)


g <- ggplot(df_plot, aes(x = Jaar, y = occupancy, colour = soortNL)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  xlim(1980, 2020) +
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

ggsave(filename = "nf_sp.eps", plot = g)
ggsave(filename = "nf_sp.png", plot = g, dpi = 300)

