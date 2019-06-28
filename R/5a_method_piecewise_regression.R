## Method 1: piecewise regression

spPR <- function(df){

  require(segmented)
  require(MASS)

  fyear <- min(df$year)
  lyear <- max(df$year)
  spec <- df[[1,1]]
  print(paste0(spec, "_", lyear))

  result <- try({
    lm <- glm.nb(ncells ~ year, data = df)

    # The automatic procedure has tendency to overestimate the number of breakpoints
    # (help file segmented - Especially when neg. binomial distribution is used.)

    lm_s1 <- segmented(obj = lm, seg.Z = ~year, ,psi = list(year = NA),
                       control = seg.control(stop.if.error = FALSE,
                                             n.boot = 0, it.max = 50, K = 5))
    #summary(lm_s1)  #plot.segmented(lm_s1)   #draw.history(lm_s1, "year")

    df_n <- data.frame(year = seq(from = fyear, to = lyear, length.out = 200))

    temp <- predict(object = lm_s1, newdata = df_n, se.fit = TRUE)

    ilink <- family(lm)$linkinv

    df_n$fit <- ilink(temp$fit)
    df_n$ucl <- ilink(temp$fit + temp$se.fit * 1.96)
    df_n$lcl <- ilink(temp$fit - temp$se.fit * 1.96)

    #ptitle <- paste0("PR/", spec, "_", lyear)
    #plot_ribbon(df_n = df_n, df = df, ptitle, printplot = TRUE)

    # Get slopes and lcl, ucl for each breakpoint
    if (length(lm_s1$coefficients) > 2){
      # breakpoints (when > 1 slope)
      msm <- as.data.frame(slope(lm_s1, conf.level = 0.9)[[1]][,c(1, 4, 5)]) %>%
        mutate(fyear = c(min(df$year), lm_s1$psi[,2])) %>%
        rename(lcl = 'CI(95%).l') %>%
        rename(ucl = 'CI(95%).u')
    }else{
      # no breakpoint (only one slope from lm)
      conflim <- confint(lm, level = 0.9)
      msm <- data.frame(Est. = lm$coefficients[2],
                        lcl = conflim[2, 1],
                        ucl = conflim[2, 2],
                        fyear = fyear)
    }

    # assign em value (-3 to 3)
    msm <- msm %>%
      mutate(em = case_when(
        lcl > 0 ~ 3,
        ucl < 0 ~ -3,
        lcl * ucl < 0 ~ 0))

    # Assign emerging status to predicted dataset
    df_n$em <- NA
    for (i in 1:nrow(msm)){
      tyear <- msm[i, "fyear"]
      tem <- msm[i, "em"]
      df_n[df_n$year > tyear, "em"] <- tem
    }

    em = msm[nrow(msm), "em"]  # emerging status for the last data point

  })

  if (!result == 0){
    em <- NA
    lm_s1 <- NULL
    df_n <- NULL
  }

  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear,
                  method_em = "PR", em = em)
  return(list(df = df, em = df_em, model = lm_s1, df_n = df_n))

}
