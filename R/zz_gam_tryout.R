### GAM tryout

library(gratia)

specid <- "2706501"
specid <- "3026004"
specid <- "8193935"

df <- filter(df_sp, taxonKey == specid)

lyear <- max(df$year)
#df$ncells <- df$ncells / 10000

#g1 <- gam(ncells ~ s(year), family = poisson, data = df)
g1 <- gam(ncells ~ s(year, k = 5), family = nb(), data = df, method = "REML")
#g1 <- gam(ncells ~ s(year, k = 5), family = gaussian, data = df, method = "REML")
# g1 <- gamm(ncells ~ s(year, k = 5), family = gaussian, data = df,
#           correlation = corCAR1(form = ~ year),
#           method = "REML")
#g1$lme

#summary(g1)  #plot(g1)

#Check residuals
#acf(g1$residuals)

temp <- predict(object = g1, newdata = df, se.fit = TRUE)
df$fit <- exp(temp$fit)
df$ucl <- exp(temp$fit + temp$se.fit * 1.96)
df$lcl <- exp(temp$fit - temp$se.fit * 1.96)

# df$fit <- temp$fit
# df$ucl <- temp$fit + temp$se.fit * 1.96
# df$lcl <- temp$fit - temp$se.fit * 1.96

spec <- df[[1,1]]
ptitle <- paste0("GAM/", spec, "_", lyear)
plot_ribbon(df, ptitle, printplot = TRUE)


#simulations
sims <- simulate(g1, nsim = 20, newdata = df, unconditional = TRUE)

colnames(sims) <- paste0("sim", seq_len(20))
sims <- setNames(stack(as.data.frame(sims)), c("simulated", "run"))
sims <- transform(sims, year = rep(df$year, 20),
                  simulated = simulated)
## Plot simulated trends
ggplot(df, aes(x = year, y = fit)) +
  geom_line(data = sims,
            mapping = aes(y = simulated, x = year, group = run),
            colour = "grey80") +
  geom_line(lwd = 1) +
  geom_point(data = df, aes(x = year, y = ncells)) +
  geom_line(data = df, aes(x = year, y = ncells))


confl <- confint(object = g1, parm = "year", newdata = data.frame(df))

ggplot(confl, aes(x = year, y = est)) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper),
              fill = "grey60", inherit.aes = FALSE) +
  geom_line(lwd = 1)

# Derivatives

small.d <- fderiv(g1, n = nrow(df))
small.sint <- with(df,
                   cbind(confint(small.d, nsim = 20,
                                 type = "simultaneous"),
                         year = year))
ggplot(small.sint, aes(x = year, y = est)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,
              fill = "black") +
  geom_line()


#### Answer from
## https://stackoverflow.com/questions/14207250/determining-derivatives-from-gam-smooth-object

library(mgcv)
#The example from predict.gam uses finite differences to approximate the derivatives of the smoothed terms

#Here is an example to do this for a single predictor model. This is more straightforward that the example from the help.

plot_ts(df, ptitle = "gam test", printplot = TRUE)

A <- gam(ncells ~ s(year), data=df, na.action=na.omit)
# new data for prediction
newDF <- with(df, data.frame(year = unique(year)))
# prediction of smoothed estimates at each unique year value
# with standard error
B <- predict(A,  newDF, type="response", se.fit=TRUE)

df <- cbind(df, fit = B$fit, se = B$se.fit)

ggplot(df, aes(x = year, y = fit)) + geom_line() +
  geom_point(aes(y = ncells), colour = "red")

# finite difference approach to derivatives following
# example from ?predict.gam

eps <- 1e-7
X0 <- predict(A, newDF, type = 'lpmatrix')

newDFeps_p <- newDF + eps

X1 <- predict(A, newDFeps_p, type = 'lpmatrix')

# finite difference approximation of first derivative
# the design matrix
Xp <- (X0 - X1) / eps

# first derivative
fd_d1 <- Xp %*% coef(A)
df <- cbind(df, fd = fd_d1[,1])
ggplot(df, aes(x = year, y = fd)) + geom_point() + geom_line()

# second derivative
newDFeps_m <- newDF - eps

X_1 <- predict(A, newDFeps_m, type = 'lpmatrix')
# design matrix for second derivative
Xpp <- (X1 + X_1 - 2*X0)  / eps^2
# second derivative
fd_d2 <- Xpp %*% coef(A)
plot(fd_d2)



