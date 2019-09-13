#### Plot functions


### Plot time series

# input df
# year, ncells

plot_ts <- function(df, y_axis = "ncells", ptitle = NULL,
                    printplot = FALSE, saveplot = FALSE){

  spec <- df[[1,1]]
  lyear <- max(df$year)

  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = get(y_axis))) + geom_line(colour = "grey") +
    geom_point() + ggtitle(ptitle) + ylab(y_axis)

  if (saveplot == TRUE) {
    dir.create("./output/incr_em/", showWarnings = FALSE)
    ggsave(filename = paste0("./output/ts/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}

### Plot proportions

plotprop <- function(df, ptitle = NULL,
                    printplot = FALSE, saveplot = FALSE){

  speckey <- df[[1,1]]
  specname <- df[[1,"spn"]]
  lyear <- max(df$year)

  if (is.null(ptitle)) {ptitle <- paste0(speckey, "_", specname, "_", lyear, "_prop")}

  g <- ggplot(df, aes(x = year, y = prop)) + geom_line(colour = "grey") +
    geom_point() + ggtitle(ptitle)

  if (saveplot == TRUE) {
    dir.create("./output/ts/", showWarnings = FALSE)
    ggsave(filename = paste0("./output/ts/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}

 ### Plot time series with confidence limits

plot_ribbon <- function(df_n, df, ptitle, printplot = FALSE, saveplot = FALSE){

  g <- ggplot(df_n, aes(x = year, y = fit)) + geom_line() +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = ncells)) +
    ggtitle(ptitle)

  if (saveplot == TRUE) {
  dir.create("./output/incr_em/", showWarnings = FALSE)
  ggsave(filename = paste0("./output/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}

### Plot time series with confidence limits + emerging status

plot_ribbon_em <- function(df_n, df, ptitle = NULL,
                           printplot = FALSE, saveplot = FALSE){

  lyear <- max(df$year)
  if (is.null(ptitle)) {
    spec <- df[[1, "taxonKey"]]
    ptitle <- paste0(spec, "_", lyear)
  }

  g <- ggplot(df_n, aes(x = year, y = fit)) +
    geom_line(colour = "grey") +
    geom_point(aes(colour = as.factor(em)), size = 2) +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = ncells)) +
    scale_colour_manual(values = c("4" = "dark red", "3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey50", "-1" = "light yellow",
                                   "-2" = "yellow", "-3" = "green", "-4" = "dark green")) +
    ggtitle(ptitle)

  # color palette to be changed

  if (saveplot == TRUE) {
    dir.create("./output/incr_em/", showWarnings = FALSE)
    ggsave(filename = paste0("./output/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}



### Plot increasing time series emerging evolution

# input df
# year, ncells, em


# Plot incrementing time series with 5 emerging status levels
plot_incr_em <- function(df, ptitle = NULL,
                          printplot = FALSE, saveplot = FALSE) {

  spec <- df[[1,"taxonKey"]]
  lyear <- max(df$year)
  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = ncells, colour = em)) +
    geom_line(colour = "grey") +
    geom_point() +
    scale_colour_manual(values = c("3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey70", "-1" = "yellow",
                                   "-2" = "green", "-3" = "dark green"),
                                   na.value = "grey30") +
    ggtitle(ptitle)

  if (saveplot == TRUE) {
    dir.create("./output/incr/", showWarnings = FALSE)
    ggsave(paste0("./output/incr/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
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


## Plot a map

plot_map <- function(df, year){

  df %>%
  filter(year == year) %>%
  mutate(x = as.integer(substr(eea_cell_code, start = 5, stop = 8)),
         y = as.integer(substr(eea_cell_code, start = 10, stop = 13))) %>%
  plot(aes(x = x, y = y)) + geom_point
}
