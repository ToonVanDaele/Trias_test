### TRIAS - Method GAM

# df dataframe with time serie for one species - grouped by year & location
#
# return list with following elements
# - em dataframe with only one row indicating emerging status for the last year
# - model the whole gam model
# - df_n data frame with prediction data
# - plot ggplot figure with original data and interpretation of emerging status
# - deriv1, deriv2 dataframes with 1st & 2nd derivatives of smoother 'year'
# - result result of try (to capture errors)

spGAM_count <- function(df, printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,"taxonKey"]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  ptitle <- paste0("GAM_count/", spec, "_", spn, "_", lyear)
  print(ptitle)

  # add x,y data
  # df <- df %>% left_join(df_xy %>%
  #                                dplyr::select(eea_cell_code, x, y, natura2000),
  #                              by = "eea_cell_code")
  #
  # df$eea_cell_code <- as.factor(df$eea_cell_code)
  #

  maxk <- max(round((lyear - fyear) / 10, 0), 4)  # 1 knot per decade, min 4

  result <- try({
      g1 <- gam(obs ~ s(year, k = maxk, m = 3, bs = "tp") +
                s(native_obs) + s(x, y, bs = "gp", k = 100, m = c(3, 10)),
                family = nb(),
              data = df, method = "REML")

    #draw(g1)
    df_n <- df
    temp <- predict(object = g1, newdata = df_n, type = "iterms", interval = "prediction",
                    se.fit = TRUE)

    # Calculate confidence intervals in link scale and backtransform to real scale
    intercept <- unname(g1$coefficients[1])
    df_n$fit <- exp(temp$fit[,1] + intercept)
    df_n$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
    df_n$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

    # Create plot with conf. interval + colour for shape status

    g <- df_n %>%
      group_by(year) %>%
      summarise(yobs = sum(obs),
                yfit = sum(fit),
                ylcl = sum(lcl),
                yucl = sum(ucl)) %>%
      ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
      geom_line(aes(y = yfit), colour = "red") +
      geom_ribbon(aes(ymax = yucl, ymin = ylcl),
                  fill = grey(0.5),
                  alpha = 0.4) +
      ggtitle(ptitle)


    # Calculate first and second derivative + conf. interval

    df_new <- data.frame(year = seq(fyear, lyear), cobs = 0, x = 0, y = 0)

    deriv1 <- derivatives(g1, term = "year", type = "central", newdata = df_new,
                          order = 1, level = 0.8, n = nrow(df_new), eps = 1e-4)

    deriv2 <- derivatives(g1, term = "year", type = "central", newdata = df_new,
                          order = 2, level = 0.8, n = nrow(df_new), eps = 1e-4)
    #draw(deriv1) #draw(deriv2)

    # Emerging status based on first and second derivative
    em_level_gam <- em_level(filter(deriv1, smooth == "s(year)"),
                             filter(deriv2, smooth == "s(year)"))
    df_new <- bind_cols(df_new, em_level_gam)

    out <- em_gam2em(em_level_gam) # get emerging status of the last year

    # p-waarde van de smoother te groot -> output NA toevoegen


  })

  if (class(result) == "try-error") {
    out <- NA
    g1 <- g <- em_level_gam <- deriv1 <- deriv2 <- NULL
  }

  df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "GAM_count",
                    em = out)

  if (savemodel == FALSE) g1 <- NULL
  return(list(em = df_em, model = g1, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = result))
}


### GAM count without spatial smooth

spGAM_count_ns <- function(df, printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,"taxonKey"]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  ptitle <- paste0("GAM_count_ns/", spec, "_", spn, "_", lyear)
  print(ptitle)

  maxk <- max(round((lyear - fyear) / 10, 0), 4)  # 1 knot per decade, min 4

  result <- try({
    g1 <- gam(obs ~ s(year, k = maxk, m = 3, bs = "tp") +
                s(native_obs),
              family = nb(),
              data = df, method = "REML")

    #draw(g1)
    df_n <- df
    temp <- predict(object = g1, newdata = df_n, type = "iterms", interval = "prediction",
                    se.fit = TRUE)

    # Calculate confidence intervals in link scale and backtransform to real scale
    intercept <- unname(g1$coefficients[1])
    df_n$fit <- exp(temp$fit[,1] + intercept)
    df_n$ucl <- exp(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
    df_n$lcl <- exp(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

    # Create plot with conf. interval + colour for shape status

    g <- df_n %>%
      group_by(year) %>%
      summarise(yobs = sum(obs),
                yfit = sum(fit),
                ylcl = sum(lcl),
                yucl = sum(ucl)) %>%
      ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
      geom_line(aes(y = yfit), colour = "red") +
      geom_ribbon(aes(ymax = yucl, ymin = ylcl),
                  fill = grey(0.5),
                  alpha = 0.4) +
      ggtitle(ptitle)


    # Calculate first and second derivative + conf. interval

    df_new <- data.frame(year = seq(fyear, lyear), cobs = 0, x = 0, y = 0)

    deriv1 <- derivatives(g1, term = "year", type = "central", newdata = df_new,
                          order = 1, level = 0.8, n = nrow(df_new), eps = 1e-4)

    deriv2 <- derivatives(g1, term = "year", type = "central", newdata = df_new,
                          order = 2, level = 0.8, n = nrow(df_new), eps = 1e-4)
    #draw(deriv1) #draw(deriv2)

    # Emerging status based on first and second derivative
    em_level_gam <- em_level(deriv1, deriv2)
    df_new <- bind_cols(df_new, em_level_gam)

    out <- em_gam2em(em_level_gam) # get emerging status of the last year

    # p-waarde van de smoother te groot -> output NA toevoegen


  })

  if (class(result) == "try-error") {
    out <- NA
    g1 <- g <- em_level_gam <- deriv1 <- deriv2 <- NULL
  }

  df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = "GAM_count_ns",
                  em = out)

  if (savemodel == FALSE) g1 <- NULL
  return(list(em = df_em, model = g1, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = result))
}


