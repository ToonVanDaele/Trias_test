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

spGAM_lcount <- function(df, native_obs = FALSE, nbyear = 3,
                         printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn))

  if (native_obs == TRUE) {
    method_em = "GAM_lcount_cobs"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp") + s(cobs))
    ptitle <- paste0("GAM_lcount_cobs/", spec, "_", spn)
  }else{
    method_em = "GAM_lcount"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
    ptitle <- paste0("GAM_lcount/", spec, "_", spn)
  }

  # assign NULL values in case something goes wrong later
  g1 <- df_n <- g <- em_level_gam <- NULL
  deriv1 <- deriv2 <- err_result <- NULL
  df_em <- data.frame(taxonKey = spec, eyear = (lyear - nbyear + 1):lyear,
                      method_em = method_em, em = NA, lcl = NA,
                      stringsAsFactors = FALSE)

  if ((lyear - fyear) > 3 & sum(df$obs[2:nrow(df)]) != 0) {

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

      # Mean lower confidence limit from the first derivative
      df_lcl <- get_lcl(deriv1, nbyear)

      # Create plot with conf. interval + colour for status

      g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                          printplot = printplot, saveplot = saveplot)

      out <- em_gam2em(em_level_gam, nbyear) # get emerging status

      df_em <- tibble(taxonKey = spec, eyear = out$year, method_em = method_em,
                      em = out$em_out)

      df_em <- df_em %>%
        left_join(df_lcl %>%
                    select(year, lcl), by = c("eyear" = "year"))

    })

    if (class(result)[1] == "try-error") err_result <- result

  }else{
    err_result <- "Insufficient data"
  }

  if (savemodel == FALSE) g1 <- NULL
  return(list(df_em = df_em, model = g1, df_n = df_n, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = err_result))
}


### GAM lumped count with ncobs (observations native species)
spGAM_lcount_cobs <- function(df, printplot = FALSE, saveplot = FALSE,
                               savemodel = FALSE, nbyear = 3){

  spGAM_lcount(df = df, native_obs = TRUE, printplot = printplot,
               saveplot = saveplot, nbyear = nbyear, savemodel = savemodel)

}


### GAM lumped presence absence (ncell)

spGAM_lpa <- function(df, native_obs = FALSE, printplot = FALSE, nbyear = 3,
                      saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)

  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn))

  if (native_obs == TRUE) {
    method_em <- "GAM_lpa_cobs"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp") + s(cobs))
    ptitle <- paste0("GAM_lpa_cobs/", spec, "_", spn, "_", lyear)
  }else{
    method_em <- "GAM_lpa"
    fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
    ptitle <- paste0("GAM_lpa/", spec, "_", spn, "_", lyear)
  }

  g1 <- df_n <- g <- em_level_gam <- NULL
  deriv1 <- deriv2 <- err_result <- NULL
  df_em <- data.frame(taxonKey = spec, eyear = (lyear - nbyear + 1):lyear,
                      method_em = method_em, em = NA, lcl = NA,
                      stringsAsFactors = FALSE)

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

      # Lower confidence limit of the last three years from the first derivative
      df_lcl <- get_lcl(deriv1, nbyear)

      # Create plot with conf. interval + colour for status

      g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                          printplot = printplot, saveplot = saveplot)

      out <- em_gam2em(em_level_gam, nbyear) # get emerging status

      df_em <- tibble(taxonKey = spec, eyear = out$year, method_em = method_em,
                      em = out$em_out)

      df_em <- df_em %>%
        left_join(df_lcl %>%
                    select(year, lcl), by = c("eyear" = "year"))

    })

    if (class(result)[1] == "try-error") err_result <- result

  }else{
    err_result <- "Insufficient data"
  }

  if (savemodel == FALSE) g1 <- NULL
  return(list(df_em = df_em, model = g1, df_n = df_n, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = err_result))
}


## GAM presence / absence with native observations (cobs)
spGAM_lpa_cobs <- function(df, native_obs = FALSE, printplot = FALSE,
                      saveplot = FALSE, savemodel = FALSE, nbyear = 3) {

  spGAM_lpa(df = df, native_obs = TRUE, printplot = printplot,
  saveplot = saveplot, savemodel = savemodel, nbyear = nbyear)
}
