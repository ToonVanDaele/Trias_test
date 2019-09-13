## Method: INLA with autocorrelation AR2

spINLA <- function(df, printplot = FALSE, saveplot = FALSE){

  require(INLA)
  require(inlatools)

  spec <- df[[1,1]]
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  fyear <- min(df$year)
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

  result <- try({
    i2 <- inla(ncells ~ f(year,
                          model = "rw2",
                          hyper = list(
                            theta = list(prior = "pc.prec", param = c(0.05, 0.01)))),
               control.compute = list(dic = TRUE),
               family = "nbinomial",
               data = df)

    #summary(i2)

    df_n <- df
    df_n$fit <- i2$summary.fitted.values[,"mean"]
    df_n$lcl <- i2$summary.fitted.values$"0.025quant"
    df_n$ucl <- i2$summary.fitted.values$"0.975quant"
    df_n$em <- 0  # No emergence status

    ptitle <- paste0("/INLA/", spec, "_", spn, "_", lyear)
    g <- plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle,
                        printplot = printplot, saveplot = saveplot)
    out <- NA   # NO emergence status

  })

  if (class(result) == "try-error") {
    df$fit <- df$ucl <- df$lcl <- out <- NA
    i2 <- df_n <- g <- NULL
  }

  df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "INLA", em = out)
  return(list(df = df, em = df_em, model = i2, df_n = df_n, plot = g, result = result))

}
