# MAIN C,D,E

### Main B

# Select species (for testing)
df_sp <- df_pp %>% filter(taxonKey == "5289686") %>%
  group_by(taxonKey)

# Get a vector with the group names
taxl <- df_sp %>%
  group_keys() %>%
  pull(taxonKey)

# plot time series
g <- df_sp %>%
  group_split() %>%
  map(.f = plot_ts, printplot = TRUE, saveplot = FALSE) %>%
  set_names(taxl)


df <- df_sp

g1 <- gam(obs ~ s(year, k = -1, m = 3, bs = "ts"), family = nb(),
          data = df, method = "REML")

summary(g1)
draw(g1)
draw(derivatives(g1, order = 1))
appraise(g1)


df_n1 <- data.frame(year = seq(from = fyear, to = lyear,
                              length.out = (lyear - fyear) * 5))
temp <- predict(object = g1, newdata = df_n1, type = "iterms", interval = "prediction",
                se.fit = TRUE)

intercept <- unname(g1$coefficients[1])
df_n1$fit <- exp(temp$fit[,1] + intercept)
df_n1$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_n1$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_n1$em <- 0

g <- plot_ribbon_em(df_n = df_n1, df = df,
                    printplot = TRUE, saveplot = FALSE)

g2 <- gam(obs ~ s(year, k = -1, m = 3, bs = "ts") +
            s(cobs, k = 3, bs = "ts"), family = nb(),
          data = df, method = "REML")

summary(g2)
draw(g2)
draw(derivatives(g2, order = 1))
appraise(g2)

df$g2fit <- g2$fitted.values


ggplot(df, aes(x = year, y = obs)) + geom_point() +
         geom_line(aes(y = g2fit))

df_n2 <- data.frame(year = seq(from = fyear, to = lyear,
                               length.out = (lyear - fyear) * 5))
temp <- predict(object = g2, newdata = df_n2, type = "iterms", interval = "prediction",
                se.fit = TRUE)

intercept <- unname(g2$coefficients[1])
df_n2$fit <- exp(temp$fit[,1] + intercept)
df_n2$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_n2$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_n2$em <- 0

g <- plot_ribbon_em(df_n = df_n2, df = df,
                    printplot = TRUE, saveplot = FALSE)




# 1690429 - 2015-2017 geven allemaal 0. Komt door breed betrouwbaarheidsinterval
# 2225772 - is traag emerging. met 'tp' smoother em(2013), maar niet met 'ts'
# 2226990 - is emerging. Ook goed resultaat met GAM
# 2227000 - hier ontbreekt echt correctie voor observatieinspanning op klasse niveau
# 2287615 - deze zou als duidelijk emerging moeten klasseren
# 2362868 - is niet emerging. Helemaal gestabiliseerd
# 2379089 - zou emerging moeten zijn. GAM geeft niet emerging
# 2426661 - niet emerging, wel gestabiliseerd. GAM afhankelijk van de smoother
# 1718308 - emerging  = ok




