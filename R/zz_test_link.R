### Probleem slope + conf interval bij log link function

df <- data.frame(ryear = seq(from = 1990, to = 2018, by = 1),
                 year = seq(from = 0, to = 28, by = 1))
set.seed(6)
df$f <- 50 + df$year * 1.7
df$ncells <- rpois(n = nrow(df), lambda = df$f)

ggplot(df, aes(x = ryear, y = ncells)) + geom_point() + geom_line()


lm1 <- glm(ncells ~ year, df, family = gaussian)
summary(lm1)
confint(lm1)

ilink <- lm1$family$linkin
temp <- predict(lm1, newdata = df, se.fit = TRUE)
df$lm1_fit <- ilink(temp$fit)
df$lm1_ucl <- ilink(temp$fit + temp$se.fit * 1.96)
df$lm1_lcl <- ilink(temp$fit - temp$se.fit * 1.96)

lm3 <- glm.nb(ncells ~ year, df)
summary(lm3)
confint(lm3)

ilink <- lm3$family$linkin
temp <- predict(lm3, newdata = df, se.fit = TRUE)
df$lm3_fit <- ilink(temp$fit)
df$lm3_ucl <- ilink(temp$fit + temp$se.fit * 1.96)
df$lm3_lcl <- ilink(temp$fit - temp$se.fit * 1.96)



ggplot(df, aes(x = ryear, y = ncells)) + geom_point() + geom_line() +
  geom_point(aes(y = lm1_fit, colour = "lm1")) +
  geom_ribbon(aes(ymin = lm1_lcl, ymax = lm1_ucl, colour = "lm1"),
              fill = "red", alpha = 0.4) +
  geom_point(aes(y = lm3_fit, colour = "lm3")) +
  geom_ribbon(aes(ymin = lm3_lcl, ymax = lm3_ucl, colour = "lm1"),
              fill = "blue", alpha = 0.4)


