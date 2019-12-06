#### Plot functions


### Plot time series

# input df
# year, ncells

plot_ts <- function(df, y_axis = "ncell", ptitle = NULL,
                    printplot = FALSE, saveplot = FALSE){

  spec <- df[[1,"taxonKey"]]
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

plot_ribbon_em <- function(df_n, df, y_axis = "obs", ptitle = NULL,
                           printplot = FALSE, saveplot = FALSE){

  if (is.null(ptitle)) ptitle <- set_title(df_n)

  g <- ggplot(df_n, aes(x = year, y = fit)) +
    geom_line(colour = "dark grey") +
    geom_point(colour = "dark grey", size = 2) +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = get(y_axis))) +
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
  if (is.null(ptitle)) ptitle <- paste0(spec, "_", lyear)

  g <- ggplot(df, aes(x = year, y = obs, colour = as.factor(em))) +
    geom_line(colour = "grey") +
    geom_point() +
    scale_colour_manual(values = c("4" = "dark red", "3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey70", "-1" = "yellow",
                                   "-2" = "green", "-3" = "green", "-4" = "green"),
                                   na.value = "grey30") +
    ggtitle(ptitle)

  if (saveplot == TRUE) {
    dir.create("./output/incr/", showWarnings = FALSE)
    ggsave(paste0("./output/incr/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}


# Plot smoother s_year
plot_smoother <- function(df, ptitle = NULL) {

  dir <- set_title(df)
  if (is.null(ptitle)) ptitle <- paste0(dir, " smoother s_year")

  df$em <- as.factor(df$em)

  g <- ggplot(df, aes(x = year, y = s_year_fit)) +
    geom_line(colour = "black") +
    geom_point(aes(colour = em), size = 2) +
    geom_ribbon(aes(ymax = s_year_ucl, ymin = s_year_lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    scale_colour_manual(values = c("4" = "dark red", "3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey50", "-1" = "light yellow",
                                   "-2" = "yellow", "-3" = "green", "-4" = "dark green")) +
    ggtitle(ptitle)

  return(g)
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

plot_map <- function(df){

  df %>%
  mutate(ob = ifelse(obs > 0, 1, 0)) %>%
  left_join(df_xy, by = "eea_cell_code") %>%
  ggplot(aes(x = x, y = y, colour = ob)) + geom_point() +
    xlim(3780000, 4080000) + ylim(2910000, 3180000) + coord_fixed()
}
