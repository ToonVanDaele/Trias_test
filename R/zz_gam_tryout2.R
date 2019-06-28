### GAM tryout 2

library(gratia)

#specid <- "2706501"
specid <- "3026004"
specid <- "8193935"
specid <- "4314683"
specid <- "3172100"

df <- filter(df_sp, taxonKey == specid)

lyear <- max(df$year)
fyear <- min(df$year)
spec <- df[[1,1]]
# Negative binomial
g_nb <- gam(ncells ~ s(year, k = -1, m = 4), family = nb(), data = df, method = "REML")

#draw(g_nb)  #appraise(g_nb)
df_n <- data.frame(year = seq(from = fyear, to = lyear, length.out = 50))

temp <- predict(object = g_nb, newdata = df_n, se.fit = TRUE)

df_n$method <- "nb"
df_n$fit <- exp(temp$fit)
df_n$ucl <- exp(temp$fit + temp$se.fit * 1.96)
df_n$lcl <- exp(temp$fit - temp$se.fit * 1.96)

ptitle <- paste0("GAM/", spec, "_", lyear)
plot_ribbon(df_n, df, ptitle, printplot = TRUE)

deriv1 <- derivatives(g_nb, type = "central", order = 1, level = 0.9, eps = 1e-5)
draw(deriv1)
deriv2 <- derivatives(g_nb, type = "central", order = 2, level = 0.9, eps = 1e-5)
draw(deriv2)



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



#### derivates with finite differences

finite.differences <- function(x, y) {

  n <- length(x)
  fdx <- vector(length = n)

  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i - 1] <- (y[i - 1] - y[i]) / (x[i - 1] - x[i])
  }

  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)
}



#g_g <- gam(ncells ~ s(year, k = 5), family = gaussian, data = df, method = "REML")
g_g <- gam(ncells ~ s(year, k = 5), family = nb(), data = df, method = "REML")

dfg <- data.frame(year = seq(from = fyear, to = lyear, length.out = 200))

temp <- predict(object = g_g, newdata = dfg, se.fit = TRUE, )

dfg$method <- "g"
dfg$fit <- temp$fit
dfg$ucl <- temp$fit + temp$se.fit * 1.96
dfg$lcl <- temp$fit - temp$se.fit * 1.96

ggplot(dfg, aes(x = year, y = fit)) + geom_line() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, fill = "black") +
  geom_point(data = df, aes(x = year, y = ncells))

dfg$fderiv <- finite.differences(x = dfg$year, y = dfg$fit)

ggplot(dfg, aes(x = year, y = fderiv)) + geom_line()





dat <- gamSim(1, n = 400, dist = "normal", scale = 2, verbose = FALSE)

plot(dat$y, dat$x2)
mod <- gam(y ~ s(x2), data = dat, method = "REML")

## first derivative of all smooths using central finite differenc
draw(mod)
mderiv1 <- derivatives(mod, type = "central", order = 1, level = 0.95)
draw(mderiv1)
mderiv2 <- derivatives(mod, type = "central", order = 2, level = 0.95)
draw(mderiv2)

mderivf1 <- derivatives(mod, type = "forward", order = 1, level = 0.95)
draw(mderivf1)
mderivf2 <- derivatives(mod, type = "forward", order = 2, level = 0.95)
draw(mderivf2)


