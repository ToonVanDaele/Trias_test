### Method GAM

spGAM <- function(df) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)

  # Check minimum 5 values > 0
  if (df %>%
    filter(ncells > 0) %>%
    count() %>% pull() >= 5) {

    # Gam model
    g1 <- gam(ncells ~ s(year, k = 5, m =4), family = nb(), data = df, method = "REML")
    #draw(g1)  #appraise(g1)

    # Predict to new data (200 values between first and last year)
    df_n <- data.frame(year = seq(from = fyear, to = lyear, length.out = 200))
    temp <- predict(object = g1, newdata = df_n, se.fit = TRUE)

    # Calculate confidence intervals in link scale and backtransform to real scale
    df_n$fit <- exp(temp$fit)
    df_n$ucl <- exp(temp$fit + temp$se.fit * 1.96)
    df_n$lcl <- exp(temp$fit - temp$se.fit * 1.96)

    # Calculate first and second derivative + conf. interval
    deriv1 <- derivatives(g1, type = "central", order = 1, level = 0.9, eps = 1e-5)
    deriv2 <- derivatives(g1, type = "central", order = 2, level = 0.9, eps = 1e-5)
    #draw(deriv1) #draw(deriv2)

    # Emerging or not, based on last year
    em <- em_level(deriv1, deriv2)
    df_n <- bind_cols(df_n, em)

    # Create plot with conf. interval + colour for emerging status
    ptitle <- paste0("GAM/", spec, "_", lyear)
    plot_ribbon_em(df_n, df, ptitle, printplot = FALSE)

    # Save emerging status of the last year
    out <- em %>%
    filter(year == max(year)) %>%
    .$em

    }else{
    df$fit <- df$ucl <- df$lcl <- out <- NA

    }
  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear, em = out)
  return(list(df = df_n, em = df_em))
}


