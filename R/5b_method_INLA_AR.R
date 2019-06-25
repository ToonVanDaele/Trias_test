## Method: INLA with autocorrelation AR2

spINLA <- function(df){

  require(INLA)
  require(inlatools)

  lyear <- max(df$year)

  #i1 <- inla(ncells ~ year,
  #           family = "nbinomial",
  #           data = df,
  #           control.compute = list(dic = TRUE, waic = TRUE))
  #summary(i1)

  #i1$waic$waic
  #i1$dic$dic

  #i1_betas <- i1$summary.fixed

  #print(i1_betas, digits = 2)

  #i1_betas$mean

  #plot(dispersion_check(i1))
  #plot(fast_distribution_check(i1))
  #plot(simulate_iid(sigma = 2))

  #rw1 <- simulate_rw(sigma = 0.1, length = 40)
  #plot(select_divergence(rw1),link = "log")

  i2 <- inla(ncells ~ f(year,
                        model = "rw2",
                        hyper = list(
                          theta = list(prior = "pc.prec", param = c(0.05, 0.01)))),
             control.compute = list(dic = TRUE),
             family = "nbinomial",
             data = df_sp)

  #summary(i2)


  i2_fit <- i2$summary.fitted.values[,"mean"]
  i2_fit_025 <- i2$summary.fitted.values$"0.025quant"
  i2_fit_975 <- i2$summary.fitted.values$"0.975quant"
  df$fit  <- i2_fit
  df$lcl <- i2_fit_025
  df$ucl <- i2_fit_975

  ptitle <- paste0("/INLA/", df[[1,1]], "_", lyear)
  plot_ribbon(df = df, ptitle = ptitle, printplot = TRUE)
  #plot_inla_rw2(df_sp, printplot = TRUE)


}
