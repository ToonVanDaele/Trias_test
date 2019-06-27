## Method 1: piecewise regression

spPR <- function(df){

  require(segmented)
  require(MASS)

  fyear <- min(df$year)
  lyear <- max(df$year)
  spec <- df[[1,1]]
  print(paste0(spec, "_", lyear))

  # Check last 5 years not all zero
  min(5, nrow(df))  # length of time series of <  5
  if (sum(tail(df$ncells, min(5, nrow(df)))) != 0) {

    lm <- glm.nb(ncells ~ year, data = df)

    # Automatic breakpoints
    # The automatic procedure has tendency to overestimate the number of breakpoints
    # (help file segmented - Especially when neg. binomial distribution is used.)

    # Try moet nog verder uitgewerkt, met opvangen van errors.!!
    try({lm_s1 <- segmented(obj = lm, seg.Z = ~year, ,psi = list(year = NA),
                       control = seg.control(stop.if.error = FALSE,
                                             n.boot = 0, it.max = 50, K = 5))})
    #summary(lm_s1)  #plot.segmented(lm_s1)   #draw.history(lm_s1, "year")

    df_n <- data.frame(year = seq(from = fyear, to = lyear, length.out = 200))

    temp <- predict(object = lm_s1, newdata = df_n, se.fit = TRUE)

    ilink <- family(lm)$linkinv

    df_n$fit <- ilink(temp$fit)
    df_n$ucl <- ilink(temp$fit + temp$se.fit * 1.96)
    df_n$lcl <- ilink(temp$fit - temp$se.fit * 1.96)

    ptitle <- paste0("PR/", spec, "_", lyear)
    plot_ribbon(df_n = df_n, df = df, ptitle, printplot = TRUE)

    # Get slopes and lcl, ucl for each breakpoint (only when > 1 slope)
    if (length(lm_s1$coefficients > 2)){

      msm <- as.data.frame(slope(lm_s1)[[1]][,c(1, 4, 5)]) %>%
        mutate(fyear = c(min(df$year), lm_s1$psi[,2])) %>%
        rename(lcl = 'CI(95%).l') %>%
        rename(ucl = 'CI(95%).u') %>%
        mutate(em = case_when(
          lcl > 0 ~ 3,
          ucl < 0 ~ -3,
          lcl * ucl < 0 ~ 0))   ## !! estimate - niet lcl!!

      em <- msm[nrow(msm), 5]
    }else{


    df$fit <- df$ucl <- df$lcl <- em <- NA
  }
  df_em <- tibble(taxonKey = df[[1,1]], eyear = lyear,
                  method_em = "PR", em = em)
  return(list(df = df, df_em = em, model = lm_s1, df_n = df_n))

}
