### TRIAS - Method GAM

# df dataframe with time serie for one species - grouped by year
#
# return list with following elements
# - em dataframe with only one row indicating emerging status for the last year
# - model the whole gam model
# - df_n data frame with prediction data
# - plot ggplot figure with original data and interpretation of emerging status
# - deriv1, deriv2 dataframes with 1st & 2nd derivatives of smoother 'year'
# - result result of try (to capture errors)

#########################################################
# GAM lumped data on time series. Only year as covariate

spGAM_lcount <- function(df, native_obs = FALSE,
                         printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn, "_", lyear))

  if (native_obs == TRUE) {
    method_em = "GAM_lcount_cobs"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp") + s(cobs))
    ptitle <- paste0("GAM_lcount_cobs/", spec, "_", spn, "_", lyear)
  }else{
    method_em = "GAM_lcount"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
    ptitle <- paste0("GAM_lcount/", spec, "_", spn, "_", lyear)
  }

  if (nrow(df) > 3 & sum(df$obs[2:nrow(df)]) != 0) {

    result <- try({

      maxk <- max(round((lyear - fyear) / 10, 0), 5)  # 1 knot per decade
      g1 <- gam(formula = fm,
                  family = nb(),
                data = df, method = "REML")

      df_n <- df  # New data for predict
      temp <- predict(object = g1, newdata = df_n, type = "iterms",
                      interval = "prediction",
                      se.fit = TRUE)

      # Calculate confidence intervals & backtransform to real scale
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

      # Emerging status based on first and second derivative
      em_level_gam <- em_level(filter(deriv1, smooth == "s(year)"),
                               filter(deriv2, smooth == "s(year)"))
      df_n <- bind_cols(df_n, em_level_gam)

      # Create plot with conf. interval + colour for status

      g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                          printplot = printplot, saveplot = saveplot)

      out <- em_gam2em(em_level_gam) # get emerging status of the last year

    })

    if (class(result) == "try-error") {
      df$fit <- df$ucl <- df$lcl <- out <- NA
      g1 <- df_n <- g <- em_level_gam <- NULL
      deriv1 <- deriv2 <- NULL
    }

  }else{
    df$fit <- df$ucl <- df$lcl <- out <- NA
    g1 <- df_n <- g <- em_level_gam <- NULL
    deriv1 <- deriv2 <- result <- NULL

  }

  # p-waarde van de smoother te groot -> output NA toevoegen

  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear, method_em = method_em,
                  em = out)
  if (savemodel == FALSE) g1 <- NULL
  return(list(em = df_em, model = g1, df_n = df_n, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = result))
}


### GAM lumped count with ncobs (observations native species)
spGAM_lcount_cobs <- function(df, printplot = FALSE, saveplot = FALSE,
                               savemodel = FALSE){

  spGAM_lcount(df = df, native_obs = TRUE, printplot = printplot,
               saveplot = saveplot, savemodel = savemodel)

}


### GAM lumped presence absence (ncell)

spGAM_lpa <- function(df, native_obs = FALSE, printplot = FALSE,
                      saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn, "_", lyear))

  if (native_obs == TRUE) {
    method_em <- "GAM_lpa_cobs"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp") + s(cobs))
    ptitle <- paste0("GAM_lpa_cobs/", spec, "_", spn, "_", lyear)
  }else{
    method_em <- "GAM_lpa"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
    ptitle <- paste0("GAM_lpa/", spec, "_", spn, "_", lyear)
  }


  if (nrow(df) > 3 & sum(df$obs[2:nrow(df)]) != 0) {

    result <- try({

      maxk <- max(round((lyear - fyear) / 10, 0), 5)  # 1 knot per decade
      g1 <- gam(formula = fm,
                family = nb(),
                data = df, method = "REML")

      # draw(g1)
      # Predict to new data (5 values per year)
      df_n <- df
      temp <- predict(object = g1, newdata = df_n, type = "iterms",
                      interval = "prediction",
                      se.fit = TRUE)

      # Calculate confidence intervals & backtransform to real scale
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

      # Emerging status based on first and second derivative
      em_level_gam <- em_level(filter(deriv1, smooth == "s(year)"),
                               filter(deriv2, smooth == "s(year)"))
      df_n <- bind_cols(df_n, em_level_gam)

      # Create plot with conf. interval + colour for status

      g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                          printplot = printplot, saveplot = saveplot)

      out <- em_gam2em(em_level_gam) # get emerging status of the last year

    })

    if (class(result) == "try-error") {
      df$fit <- df$ucl <- df$lcl <- out <- NA
      g1 <- df_n <- g <- em_level_gam <- NULL
      deriv1 <- deriv2 <- NULL
    }

  }else{
    df$fit <- df$ucl <- df$lcl <- out <- NA
    g1 <- df_n <- g <- em_level_gam <- NULL
    deriv1 <- deriv2 <- result <- NULL

  }

  # p-waarde van de smoother te groot -> output NA toevoegen

  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear, method_em = method_em,
                  em = out)
  if (savemodel == FALSE) g1 <- NULL
  return(list(em = df_em, model = g1, df_n = df_n, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = result))
}


## GAM presence / absence with native observations (cobs)
spGAM_lpa_cobs <- function(df, native_obs = FALSE, printplot = FALSE,
                      saveplot = FALSE, savemodel = FALSE) {

  spGAM_lpa(df = df, native_obs = TRUE, printplot = printplot,
  saveplot = saveplot, savemodel = savemodel)
}
