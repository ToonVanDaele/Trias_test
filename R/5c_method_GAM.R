### Method GAM

spGAM <- function(df, printplot = FALSE, saveplot = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)
  print(paste0(spec, "_", lyear))

  result <- try({

    maxk <- max(round((lyear - fyear) / 10, 0), 3)  # 1 knot per decade
    mvalue <- min(maxk - 1, 2)
    g1 <- gam(ncells ~ s(year, k = maxk, m = mvalue, bs = "ts"), family = nb(),
              data = df, method = "REML")
    #draw(g1)  #appraise(g1)
    # Predict to new data (200 values between first and last year)
    df_n <- data.frame(year = seq(from = fyear, to = lyear,
                                  length.out = (lyear - fyear) * 5))
    temp <- predict(object = g1, newdata = df_n, type = "iterms", se.fit = TRUE)

    # Calculate confidence intervals in link scale and backtransform to real scale
    intercept <- unname(g1$coefficients[1])
    df_n$fit <- exp(temp$fit[,1] + intercept)
    df_n$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
    df_n$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

    # Calculate first and second derivative + conf. interval
    deriv1 <- derivatives(g1, type = "central", order = 1, level = 0.8,
                          n = nrow(df_n), eps = 1e-4)
    deriv2 <- derivatives(g1, type = "central", order = 2, level = 0.8,
                          n = nrow(df_n), eps = 1e-4)
    #draw(deriv1) #draw(deriv2)

    # Emerging or not, based on last year
    em <- em_level(deriv1, deriv2)
    df_n <- bind_cols(df_n, em)

    # Create plot with conf. interval + colour for emerging status
    ptitle <- paste0("GAM/", spec, "_", lyear)
    g <- plot_ribbon_em(df_n = df_n, df = df, ptitle = ptitle,
                        printplot = printplot, saveplot = saveplot)

    # Save emerging status of the last year
    out <- em %>%
    filter(year == max(year)) %>%
    .$em

  })

  if (!result == 0){
    df$fit <- df$ucl <- df$lcl <- out <- NA
    g1 <- NULL
    df_n <- NULL
    deriv1 <- deriv2 <- NULL
  }

  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear, method_em = "GAM", em = out)
  return(list(df = df, em = df_em, model = g1, df_n = df_n,
              deriv1 = deriv1, deriv2 = deriv2, plot = g))
}


