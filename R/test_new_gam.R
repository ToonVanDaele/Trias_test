### Test new gam method

library(tidyverse)
library(mgcv)
library(gratia)
## df_s

# random effect for cell

# Load data
df_s <- readRDS(file = "./data/df_s.RDS")
df_pp <- readRDS(file = "./data/df_pp.RDS")
df_xy <- readRDS(file = "./data/df_xy.RDS")
spec_names <- readRDS(file = "./data/spec_names.RDS")


# Select species
spn <- "1718308"

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
df_ss <- df_ss %>%
  mutate(cobs = ifelse(cobs < obs, obs, cobs))

# For a start we only work with data > 1980. (Almost) no zeros.
df_ss <- filter(df_ss, year > 1980)

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
appraise(mA)

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


# Model C - negative binomial - offset class observations

# How does cobs look like?
df_temp <- df_ss %>%
  group_by(year) %>%
  summarise(yobs = sum(obs),
            ycobs = sum(cobs))

ggplot(df_temp, aes(x = year, y = ycobs)) + geom_point() + geom_line()
ggplot(df_temp, aes(x = yobs, y = ycobs)) + geom_point()


mC <- gam(obs ~ s(year, bs = "tp") + s(cobs), family = nb(),
          data = df_ss, method = "REML")

summary(mC)
plot(mC)
appraise(mC)

# How does it fit?
temp <- predict(object = mC, newdata = df_ss, type = "iterms",
                interval = "prediction", se.fit = TRUE)
intercept <- unname(mC$coefficients[1])
df_mC <- df_ss
df_mC$fit <- exp(temp$fit[,1] + intercept)
df_mC$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
df_mC$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

df_mC %>%
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

temp <- df_ss %>%
  group_by(eea_cell_code) %>%
  summarise(omin = min(obs),
            omax = max(obs))

#?? Hoe komt het dat voor bepaalde cellen de maximum obs = 0?

mD <- gam(obs ~ s(eea_cell_code, bs = "re"), family = nb(),
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


