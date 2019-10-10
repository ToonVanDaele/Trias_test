### Presentation core group meeting

library(tidyverse)
library(INBOtheme)

source(file = "./R/5d_method_GAM_pa.R")
source(file = "./R/5c_method_GAM_short.R")
source(file = "./R/9_function.R")
source(file = "./R/9b_plot_function.R")

df_pp <- readRDS(file = "./data/df_pp.RDS")
df_s <- readRDS(file = "./data/df_s.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")

spn <- "3084022"

p1 <- df_s %>%
  filter(taxonKey == spn) %>%
  left_join(spec_names %>%
              select(taxonKey, spn), by = c("taxonKey")) %>%
  group_by(taxonKey, year) %>%
  summarise(occurence = sum(obs),
            occupancy = sum(pa_obs)) %>%
  gather(key = type, value = count, -year, - taxonKey) %>%
  ggplot(aes(x = year, y = count, colour = type)) + geom_point() + geom_line() +
  theme_set(theme_inbo(transparent = TRUE)) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())

plot(p1)
ggsave("./plot_1.png", plot = p1)

cbPalette <- c("#cccccc", "#000000")

p2 <- df_s %>%
  filter(taxonKey == spn & year %in% c(2000, 2005, 2010, 2015)) %>%
  left_join(spec_names %>%
              select(taxonKey, spn), by = c("taxonKey")) %>%
  left_join(df_xy, by = "eea_cell_code") %>%
  arrange(pa_obs) %>%
  ggplot(aes(x = x, y = y, colour = as.factor(pa_obs))) + geom_point(size = 1) +
  coord_fixed() +
  facet_wrap(~year, ncol = 2) +
  theme(legend.position = "none") +
  scale_colour_manual(values=cbPalette)


plot(p2)
ggsave("./plot_2.png", plot = p2)



m1 <- df_s %>%
  filter(taxonKey == spn) %>%
  left_join(df_xy %>%
              select(eea_cell_code, x, y),
            by = "eea_cell_code") %>%
  spGAM_pa(savemodel = TRUE)


m2 <- df_pp %>%
  filter(taxonKey == spn & year < 2016) %>%
  spGAM_lcount(savemodel = TRUE)

p3 <- m2$plot
plot(p3)
ggsave("./plot_3.png", plot = p3)



p4a <- draw(m2$model)
plot(p4a)
ggsave("./plot_4a.png", plot = p4a)

p4b <- draw(m2$deriv1)
plot(p4b)
ggsave("./plot_4b.png", plot = p4b)

p4c <- draw(m2$deriv2)
plot(p4c)
ggsave("./plot_4c.png", plot = p4c)

