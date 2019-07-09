# Functions

# Evaluate level of emerging

em_level <- function(df1, df2){

  em1 <- df1 %>%
    mutate(em1 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em1)

  em2 <- df2 %>%
    mutate(em2 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em2)

  # Optie om hier continue score van te maken?
  # Combineren van 1ste en 2de afgeleide

  em <- bind_cols(em1, em2) %>%
    mutate(em = case_when(
      em1 == 1   & em2 == 1  ~  "4",
      em1 == 1   & em2 == 0  ~  "3",
      em1 == 1   & em2 == -1 ~  "2",
      em1 == 0   & em2 == 1  ~  "1",
      em1 == 0   & em2 == 0  ~  "0",
      em1 == 0   & em2 == -1  ~ "-1",
      em1 == -1  & em2 == 1  ~  "-2",
      em1 == -1  & em2 == 0  ~  "-3",
      em1 == -1  & em2 == -1  ~ "-4")) %>%
    dplyr::select(-year1)
  return(em)

}

# Plot segmented regression

plot_sr <- function(msm, df, printplot = FALSE){

  spec <- df[1,"taxonKey"]  # Species name/code
  #minyear <- min(df_sp$year)
  lyear <- max(df$year)
  mncells <- filter(df, year == lyear) %>% .$ncells

  g <- ggplot(df_sp, aes(x = year, y = ncells)) + geom_point()

  for (i in 1:(nrow(msm) - 1)) {
    sgm_xb <- msm[i, 5]
    sgm_xe <- msm[i + 1, 5]
    sgm_yb <- msm[i, 4] + sgm_xb * msm[i,1]
    sgm_ye <- sgm_yb + msm[i,1] * (sgm_xe - sgm_xb)

    if (msm[i,2] > 0) {mycol = "red"}else{mycol = "green"}
    xb <- sgm_xb #* attr(df_sp$year, 'scaled:scale') + attr(df_sp$year, 'scaled:center')
    xe <- sgm_xe #* attr(df_sp$year, 'scaled:scale') + attr(df_sp$year, 'scaled:center')
    g <- g + geom_segment(x = xb, y = exp(sgm_yb), xend = xe, yend = exp(sgm_ye), colour = mycol)
  }

  g <- g + ggtitle(paste0(spec, "_year: ", lyear, "_max_cells: ", mncells))

  ggsave(filename = paste0("./output/figures/sr_", spec, "_", maxyear, ".png"), plot = g)
  if (printplot == TRUE) plot(g)

}

# Plot INLA RW2 time series

plot_inla_rw2 <- function(df_sp, printplot = FALSE){

  spec <- df_sp[1,"taxonKey"]  # Species name/code
  maxyear <- max(df_sp$year)
  maxaantal <- max(df_sp$ncells)

  g <- ggplot(data = df_sp, aes(y = ncells, x = year)) +
    xlab("Year") + ylab("Aantal") +
    theme(text = element_text(size=15)) +
    geom_point(shape = 16, size = 2, col = "black") +
    geom_line(aes(x = year, y = fit)) +
    geom_ribbon(aes(x = year,
                    ymax = fit_975,
                    ymin = fit_025),
                fill = grey(0.5),
                alpha = 0.4) +
    ggtitle(paste0(spec, "_till year: ", maxyear, "_max_cells: ", maxaantal))

  ggsave(filename = paste0("./output/figures/inla_rw2_", spec, "_", maxyear, ".png"), plot = g)
  if (printplot == TRUE) plot(g)


}
