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

spGAM_lcount <- function(df, method_em = "GAM_lcount", nbyear = 3,
                         printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn))

  fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
  if (method_em == "GAM_lcount_cobs") fm <- update(fm, ~. + s(cobs))

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

      # Check at p-value of least 1 smoother < 0.1
      s_pv <- summary.gam(g1)$s.pv
      p_ok <- ifelse(any(s_pv < 0.1), TRUE, FALSE)

      if (p_ok){
        # Predict in real scale
        df_n <- predict_real_scale(df, g1)

        # Calculate first and second derivative + conf. interval
        deriv1 <- derivatives(g1, term = "s(year)", type = "central", order = 1,
                              level = 0.8, n = nrow(df_n), eps = 1e-4)
        deriv2 <- derivatives(g1, term = "s(year)", type = "central", order = 2,
                              level = 0.8, n = nrow(df_n), eps = 1e-4)

        # Emerging status based on first and second derivative
        em_level_gam <- em_level(deriv1, deriv2)
        df_n <- bind_cols(df_n, em_level_gam)

        # Mean lower confidence limit from the first derivative
        df_lcl <- get_lcl(df_deriv = deriv1, nbyear = nbyear, fam = g1$family)

        out <- em_gam2em(em_level_gam, nbyear) # get emerging status

        df_em <- tibble(taxonKey = spec, eyear = out$year, method_em = method_em,
                        em = out$em_out) %>%
          left_join(df_lcl %>%
                      select(year, lcl), by = c("eyear" = "year"))

        # Create plot with conf. interval + colour for status
        ptitle <- paste0(method_em,"/", spec, "_", spn)
        g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                            printplot = printplot, saveplot = saveplot)
      }
    })

    if (class(result)[1] %in% c("simpleWarning", "simpleError", "try-error"))
      err_result <- result

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

  spGAM_lcount(df = df, method_em = "GAM_lcount_cobs", printplot = printplot,
               saveplot = saveplot, nbyear = nbyear, savemodel = savemodel)
}

### GAM lumped presence absence (ncell)

spGAM_lpa <- function(df, method_em = "GAM_lpa", printplot = FALSE, nbyear = 3,
                      saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,1]]     # species name (taxonKey)

  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  print(paste0(spec, "_", spn))

  fm <- formula(obs ~ s(year, k = maxk, m = 3, bs = "tp"))
  if (method_em == "GAM_lpa_cobs") fm <- update(fm, ~. + s(cobs))

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

      # Check at p-value of least 1 smoother < 0.1
      s_pv <- summary.gam(g1)$s.pv
      p_ok <- ifelse(any(s_pv < 0.1), TRUE, FALSE)

      if (p_ok){
        # Predict in real scale
        df_n <- predict_real_scale(df, g1)

        # Calculate first and second derivative + conf. interval
        deriv1 <- derivatives(g1, term = "s(year)", type = "central", order = 1,
                              level = 0.8, n = nrow(df_n), eps = 1e-4)
        deriv2 <- derivatives(g1, term = "s(year)", type = "central", order = 2,
                              level = 0.8, n = nrow(df_n), eps = 1e-4)

        # Emerging status based on first and second derivative
        em_level_gam <- em_level(deriv1, deriv2)
        df_n <- bind_cols(df_n, em_level_gam)

        # Mean lower confidence limit from the first derivative
        df_lcl <- get_lcl(df_deriv = deriv1, nbyear = nbyear, fam = g1$family)

        # Create plot with conf. interval + colour for status
        ptitle <- paste0(method_em,"/", spec, "_", spn, "_", lyear)
        g <- plot_ribbon_em(df_n = df_n, y_axis = "obs", df = df, ptitle = ptitle,
                            printplot = printplot, saveplot = saveplot)

        out <- em_gam2em(em_level_gam, nbyear) # get emerging status

        df_em <- tibble(taxonKey = spec, eyear = out$year, method_em = method_em,
                        em = out$em_out) %>%
          left_join(df_lcl %>%
                      select(year, lcl), by = c("eyear" = "year"))
      }
    })

    if (class(result)[1] %in% c("simpleWarning", "simpleError", "try-error"))
      err_result <- result

  }else{
    err_result <- "Insufficient data"
  }

  if (savemodel == FALSE) g1 <- NULL
  return(list(df_em = df_em, model = g1, df_n = df_n, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = err_result))
}

## GAM presence / absence with native observations (cobs)
spGAM_lpa_cobs <- function(df, printplot = FALSE,
                      saveplot = FALSE, savemodel = FALSE, nbyear = 3) {

  spGAM_lpa(df = df, method_em = "GAM_lpa_cobs", printplot = printplot,
  saveplot = saveplot, savemodel = savemodel, nbyear = nbyear)
}
