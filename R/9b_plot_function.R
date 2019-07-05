#### Plot functions


### Plot time series

# input df
# year, ncells

plot_ts <- function(df, ptitle = NULL,
                    printplot = FALSE, saveplot = FALSE){

  spec <- df[[1,1]]
  lyear <- max(df$year)

  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = ncells)) + geom_line(colour = "grey") +
    geom_point() + ggtitle(ptitle)

  if (saveplot == TRUE) {ggsave(filename = paste0("./output/ts/", ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}

 ### Plot time series with confidence limits

plot_ribbon <- function(df_n, df, ptitle, printplot = FALSE){

  g <- ggplot(df_n, aes(x = year, y = fit)) + geom_line() +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = ncells)) +
    ggtitle(ptitle)

  ggsave(filename = paste0("./output/", ptitle, ".png"), g)
  if (printplot == TRUE) plot(g)
  return(g)
}

### Plot time series with confidence limits + emerging status

plot_ribbon_em <- function(df_n, df, ptitle = NULL,
                           printplot = FALSE, saveplot = FALSE){

  spec <- df[[1, "taxonKey"]]
  lyear <- max(df$year)
  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df_n, aes(x = year, y = fit)) + geom_point(aes(colour = as.factor(em)), size = 2) +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = ncells)) +
    scale_colour_manual(values = c("3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey50", "-1" = "yellow",
                                   "-2" = "green", "-3" = "dark green")) +
    ggtitle(ptitle)

  if (saveplot == TRUE) ggsave(filename = paste0("./output/", ptitle, ".png"), g)
  if (printplot == TRUE) plot(g)
  return(g)
}



### Plot increasing time series emerging evolution

# input df
# year, ncells, em

plot_incr_em <- function(df, ptitle = NULL, printplot = FALSE, saveplot = FALSE) {

  spec <- df[[1,"taxonKey"]]
  lyear <- max(df$year)
  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = ncells, colour = em)) + geom_point() +
      geom_line(colour = "grey") +
      scale_color_manual(values = c("dark green", "dark red")) +
      ggtitle(ptitle)

  if (saveplot == TRUE) {ggsave(filename = paste0("./output/incr_em/",
                                                  ptitle, ".png"), g)}
  if (printplot == TRUE) plot(g)
  return(g)
}

# Plot incrementing time series with 5 emerging status levels
plot_incr_em5 <- function(df, ptitle = NULL,
                          printplot = FALSE, saveplot = FALSE) {

  spec <- df[[1,"taxonKey"]]
  lyear <- max(df$year)
  if (is.null(ptitle)) {ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = ncells, colour = em)) + geom_point() +
    geom_line(colour = "grey") +
    scale_colour_manual(values = c("3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey60", "-1" = "yellow",
                                   "-2" = "green", "-3" = "dark green"),
                                   na.value = "light grey") +
    ggtitle(ptitle)

  if (saveplot == TRUE) ggsave(paste0("./output/incr/", ptitle, ".png"), g)
  if (printplot == TRUE) plot(g)
  return(g)
}


