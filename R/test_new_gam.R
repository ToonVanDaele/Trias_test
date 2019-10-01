### Test gam methods

library(tidyverse)
library(mgcv)
library(gratia)

# Load data
df_s <- readRDS(file = "./data/df_s.RDS")
df_pp <- readRDS(file = "./data/df_pp.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")

# Select species
 spn <- "1718308"
 spn <- "1690429"
spn <- "2219863"  #piek - volledig weg door class correction
# spn <- "2362868"
#spn <- "2439261"  #grote dataset (134000 records). "inner loop 1; can't correct step size
spn <- "2706056"  #goede tijdreeks met niet al te veel zeros

df_ss <- filter(df_s, taxonKey == spn)
df_sp <- filter(df_pp, taxonKey == spn)

df_ss$eea_cell_code <- as.factor(df_ss$eea_cell_code)

head(df_ss)
head(df_sp)

# Check consistency both datasets
ggplot(df_sp, aes(x = year, y = obs)) + geom_point() + geom_line()

df_ss %>%
  group_by(year) %>%
  summarise(yobs = sum(obs)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line()

# Some minor problem with cobs < obs (should never be the case)
# should be solved at the preprocessing stage

temp <- df_ss %>%
  mutate(oops = ifelse(cobs < obs, TRUE, FALSE))
summary(temp)

dd <- temp %>%
  filter(oops == TRUE) %>%
  group_by(eea_cell_code) %>%
  summarise(grid = n())

df_ss <- df_ss %>%
  mutate(cobs = ifelse(cobs < obs, obs, cobs))

# Plot reeks voor één enkele cel
df_grid <- df_ss %>%
  filter(obs > 0) %>%
  group_by(eea_cell_code) %>%
  summarise(minobs = min(obs),
            maxobs = max(obs),
            nb = n())

df_ss %>%
  filter(eea_cell_code == "1kmE3839N3112") %>%
  dplyr::select(year, obs, cobs) %>%
  gather(key = type, value = obs, - year) %>%
  ggplot(aes(x = year, y = obs, colour = type)) + geom_point() + geom_line()

# dikkere lijnen / assen grotere labels

# GAM
# Model A - gaussian
mA <- gam(obs ~ s(year, bs = "tp"), data = df_ss, method = "REML")
summary(mA)
plot(mA)

temp <- predict(object = mA, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mA$coefficients[1])
df_mA <- df_ss
df_mA$fit <- temp$fit[,1] + intercept
df_mA$ucl <- temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96
df_mA$lcl <- temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96

summary(df_mA)
#appraise(mA)

# How does it fit?
df_mA %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)

# Model B - negative binomial
mB <- gam(obs ~ s(year, bs = "tp"), family = nb(),
          data = df_ss, method = "REML")
summary(mB)
plot(mB)
appraise(mB)

# How does it fit?
temp <- predict(object = mB, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mB$coefficients[1])
df_mB <- df_ss
df_mB$fit <- exp(temp$fit[,1] + intercept)
df_mB$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mB$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mB %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)


# Model C - negative binomial - correct for class observations

# How does cobs look like?
df_temp <- df_ss %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            ycobs = sum(cobs))

ggplot(df_temp, aes(x = year, y = ycobs)) + geom_point() + geom_line()
ggplot(df_temp, aes(x = yobs, y = ycobs)) + geom_point()

# C1 -> class observations as smoother
mC1 <- gam(obs ~ s(year, bs = "tp") + s(cobs, k = 3), family = nb(),
          data = df_ss, method = "REML")

summary(mC1)
plot(mC1)
appraise(mC1)

# How does it fit?
temp <- predict(object = mC1, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mC1$coefficients[1])
df_mC1 <- df_ss
df_mC1$fit <- exp(temp$fit[,1] + intercept)
df_mC1$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mC1$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mC1 %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)

# C2 -> class observations as linear predictor
mC2 <- gam(obs ~ s(year, bs = "tp") + cobs, family = nb(),
           data = df_ss, method = "REML")

summary(mC2)
plot(mC2)
appraise(mC2)

# How does it fit?
temp <- predict(object = mC2, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mC2$coefficients[1])
df_mC2 <- df_ss
df_mC2$fit <- exp(temp$fit[,1] + intercept)
df_mC2$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mC2$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mC2 %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  filter(!year == 1990) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)


# C3 -> class observations as offset

mC3 <- gam(obs ~ offset(cobs) + s(year, bs = "tp"), family = nb(),
           data = df_ss, method = "REML")

summary(mC3)
plot(mC3)
appraise(mC3)


# How does it fit?
temp <- predict(object = mC3, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mC3$coefficients[1])
df_mC3 <- df_ss
df_mC3$fit <- exp(temp$fit[,1] + intercept)
df_mC3$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mC3$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mC3 %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)


# Model D - negative binomial - random effect cell ID

length(unique(df_ss$eea_cell_code))

mD <- gam(obs ~ offset(cobs) + s(eea_cell_code, bs = "re") + s(year), family = nb(),
          data = df_ss, method = "REML")

summary(mD)
plot(mD)
appraise(mD)

# How does it fit?
temp <- predict(object = mD, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mD$coefficients[1])
df_mD <- df_ss
df_mD$fit <- exp(temp$fit[,1] + intercept)
df_mD$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mD$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mD %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)



# Model E: with s(x,y) smoother

length(unique(df_ss$eea_cell_code))

df_ss <- df_ss %>% left_join(df_xy %>%
                      dplyr::select(eea_cell_code, x, y, natura2000),
                    by = "eea_cell_code")

df_ss$eea_cell_code <- as.factor(df_ss$eea_cell_code)

df_ss %>%
  filter(obs > 0) %>%
  group_by(eea_cell_code, year) %>%
  summarise(nc = n(),
            x = first(x),
            y = first(y)) %>%
  ggplot(aes(x = x, y = y, colour = nc)) + geom_point() + facet_wrap(~year)

mE <- gam(obs ~ s(year) + s(cobs, k = 3) + s(x, y, bs = "gp", k = 100, m = c(3, 10000)), family = nb(),
          data = df_ss, method = "REML")

summary(mE)
plot(mE)
appraise(mE)

# How does it fit?
temp <- predict(object = mE, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mE$coefficients[1])
df_mE <- df_ss
df_mE$fit <- exp(temp$fit[,1] + intercept)
df_mE$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mE$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mE %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            yfit = sum(fit),
            ylcl = sum(lcl),
            yucl = sum(ucl)) %>%
  ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
  geom_line(aes(y = yfit), colour = "red") +
  geom_ribbon(aes(ymax = yucl, ymin = ylcl),
              fill = grey(0.5),
              alpha = 0.4)






## INLA
library(INLA)

inlaA <- inla(obs ~ s(year, bs = "tp"), data = df_ss, method = "REML")
summary(inlamA)
plot(mA)

temp <- predict(object = mA, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mA$coefficients[1])
df_mA <- df_ss
df_mA$fit <- temp$fit[,1] + intercept
df_mA$ucl <- temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96
df_mA$lcl <- temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96

summary(df_mA)
#appraise(mA)
