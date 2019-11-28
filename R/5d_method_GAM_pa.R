### TRIAS - Method GAM

# df dataframe with time serie for one species - grouped by year & location -
# as presence / absence data
#
# return list with follwing elements:
# - em dataframe with only one row indicating emerging status for the last year
# - model the whole gam model
# - df_n prediction data frame
# - plot ggplot figure with original data and interpretation of emerging status
# - deriv1, deriv2 dataframes with 1st & 2nd derivatives of smoother 'year'
# - result result of try (to capture errors)

spGAM_pa <- function(df, method_em = "GAM_pa",
                     printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  require(mgcv)
  require(gratia)

  fyear <- min(df$year) # First year
  lyear <- max(df$year) # Last year
  spec <- df[[1,"taxonKey"]]     # species name (taxonKey)
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn) %>% as.character()
  ptitle <- paste0(method_em, "/", spec, "_", spn, "_", lyear)
  print(ptitle)
  start <- date()

  fm <- formula(pa_obs ~ s(year, k = maxk, m = 3, bs = "tp") + s(cobs))
  if (method_em == "GAM_pa")
    fm <- update(fm, ~. + s(x, y, bs = "gp", k = 100, m = c(3, 10)))

  g1 <- df_n <- g <- em_level_gam <- NULL
  deriv1 <- deriv2 <- err_result <- NULL
  df_em <- data.frame(taxonKey = spec, eyear = (lyear - nbyear + 1):lyear,
                      method_em = method_em, em = NA, lcl = NA,
                      stringsAsFactors = FALSE)

  maxk <- max(round((lyear - fyear) / 10, 0), 4)  # 1 knot per decade, min 4

  result <- try({
    g1 <- gam(formula = fm, family = "binomial", data = df, method = "REML")

    # Check at p-value of least 1 smoother < 0.1
    s_pv <- summary.gam(g1)$s.pv
    p_ok <- ifelse(any(s_pv < 0.1), TRUE, FALSE)

    if (p_ok){
      # Predict in real scale
      df_n <- predict_real_scale(df, g1)

      # Calculate first and second derivative + conf. interval
      df_new <- data.frame(year = seq(fyear, lyear), cobs = 0, x = 0, y = 0)

      # Calculate first and second derivative + conf. interval
      deriv1 <- derivatives(g1, term = "s(year)", type = "central", order = 1,
                            level = 0.8, n = nrow(df_n), eps = 1e-4)
      deriv2 <- derivatives(g1, term = "s(year)", type = "central", order = 2,
                            level = 0.8, n = nrow(df_n), eps = 1e-4)

      # Emerging status based on first and second derivative
      em_level_gam <- em_level(deriv1, deriv2)
      df_new <- bind_cols(df_new, em_level_gam)

      out <- em_gam2em(em_level_gam) # get emerging status of the last year

      # Mean lower confidence limit of the last three years from the first derivative
      df_lcl <- get_lcl(df_deriv = deriv1, nbyear = nbyear, fam = g1$family)

      # Create plot with conf. interval + colour for shape status
      g <- df_n %>%
        group_by(year) %>%
        summarise(yobs = sum(pa_obs),
                  yfit = sum(fit),
                  ylcl = sum(lcl),
                  yucl = sum(ucl)) %>%
        ggplot(aes(x = year, y = yobs)) + geom_point() + geom_line() +
        geom_line(aes(y = yfit), colour = "red") +
        geom_ribbon(aes(ymax = yucl, ymin = ylcl),
                    fill = grey(0.5),
                    alpha = 0.4) +
        ggtitle(ptitle)

      df_em <- tibble(taxonKey = spec, eyear = lyear, method_em = method_em,
                      em = out) %>%
        left_join(df_lcl %>%
                    select(year, lcl), by = c("eyear" = "year"))
    }
  })

  if (class(result)[1] %in% c("simpleWarning", "simpleError", "try-error"))
    err_result <- result

  end <- date()
  if (savemodel == FALSE) g1 <- NULL
  return(list(em = df_em, model = g1, em_level_gam = em_level_gam,
              deriv1 = deriv1, deriv2 = deriv2, plot = g, result = result,
              time = list(start, end)))
}


## GAM presence absence no spatial

spGAM_pa_ns <- function(df, printplot = FALSE, saveplot = FALSE, savemodel = FALSE) {

  spGAM_pa(df = df, method_em = "GAM_pa_ns", printplot = printplot,
           saveplot = saveplot, savemodel = savemodel)
}
