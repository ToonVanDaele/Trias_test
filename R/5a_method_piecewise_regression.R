## Method 1: piecewise regression

spPR <- function(df){


  print(df[[1,1]])
  require(segmented)
  require(MASS)

  maxyear <- max(df$year)

  # Check last 5 years not all zero
  min(5, nrow(df))  # length of time series of <  5
  if (sum(tail(df$ncells, min(5, nrow(df)))) != 0) {

    #lm <- lm(ncells ~ year, data = df)
    #lm <- glm(ncells ~ year, family = "poisson", data = df)
    lm <- glm.nb(ncells ~ year, data = df)

    #one break point
    #lm_s1 <- segmented(obj = lm, seg.Z = ~year, control = seg.control(stop.if.error = TRUE))

    # Automatic breakpoints
    # The automatic procedure has tendency to overestimate the number of breakpoints
    # (help file segmented)
    # (Especially when neg. binomial distribution is used.)

    # Try moet nog verder uitgewerkt, met opvangen van errors.!!
    try({lm_s1 <- segmented(obj = lm, seg.Z = ~year, ,psi = list(year = NA),
                       control = seg.control(stop.if.error = FALSE,
                                             n.boot = 0, it.max = 50, K = 5))})
    #summary(lm_s1)
    #plot.segmented(lm_s1)
    #draw.history(lm_s1, "year")

    #sr_fit <- lm_s1$fitted.values

    #df_new <- data.frame(year = seq(from = minyear, to = maxyear, by = 1))
    #pred <- predict(object = lm_s1, newdata = df_new, se.fit = TRUE)
    #df_new$pred <- pred$fit
    #df_new$se <- pred$se.fit

    temp <- predict(object = lm_s1, newdata = df, se.fit = TRUE)
    df$fit <- exp(temp$fit)
    df$ucl <- exp(temp$fit + temp$se.fit * 1.96)
    df$lcl <- exp(temp$fit - temp$se.fit * 1.96)

    spec <- df[[1,1]]
    ptitle <- paste0(spec, "_PR_", maxyear)
    plot_ribbon(df, ptitle, printplot = TRUE)


    # If there are no breakpoints
    if (length(lm_s1$coefficients == 2)){
          em <- lm_s1$coefficients[2] > 0
    }else{

      msm <- as.data.frame(slope(lm_s1)[[1]][,c(1, 4, 5)]) %>%
        mutate(fyear = c(min(df$year), lm_s1$psi[,2])) %>%
        mutate(em = .$Est. > 0)   ## !! estimate - niet lcl!!

      em <- msm[nrow(msm), 5]
    }
  }else{
    df$fit <- df$ucl <- df$lcl <- em <- NA
  }
  df_em <- tibble(taxonKey = df[[1,1]], eyear = maxyear,
                  methodem = "PR", em = em)
  return(list(df = df, em = em))

}
