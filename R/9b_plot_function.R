#### Plot functions


### Plot time series

# input df
# year, ncells

plot_ts <- function(df, ptitle = NULL, printplot = FALSE){

  spec <- df[[1,1]]
  lyear <- max(df$year)

  if (is.null(ptitle)){ptitle <- paste0(spec, "_", lyear)}

  g <- ggplot(df, aes(x = year, y = ncells)) + geom_line(colour = "grey") +
    geom_point() + ggtitle(ptitle)

  ggsave(filename = paste0("./output/", ptitle, ".png"), g)
  if (printplot == TRUE) plot(g)
  #return(g)
  #return(NULL)
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

plot_ribbon_em <- function(df_n, df, ptitle, printplot = FALSE){

  g <- ggplot(df_n, aes(x = year, y = fit)) + geom_point(aes(colour = as.factor(em)), size = 2) +
    geom_ribbon(aes(ymax = ucl, ymin = lcl),
                fill = grey(0.5),
                alpha = 0.4) +
    geom_point(data = df, aes(x = year, y = ncells)) +
    scale_colour_manual(values = c("3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey50", "-1" = "yellow",
                                   "-2" = "green", "-3" = "dark green")) +
    ggtitle(ptitle)

  ggsave(filename = paste0("./output/", ptitle, ".png"), g)
  if (printplot == TRUE) plot(g)
}



### Plot increasing time series emerging evolution

# input df
# year, ncells, em

plot_incr_em <- function(df, spec, printplot = FALSE) {

  #spec <- as.character(df[1,"taxonKey"])
  lyear <- max(df$year)
  mncells <- filter(df, year == lyear) %>% .$ncells

  g <- ggplot(df, aes(x = year, y = ncells, colour = em)) + geom_point() +
      geom_line(colour = "grey") +
      scale_color_manual(values = c("dark green", "dark red")) +
      ggtitle(paste0(spec, "_", lyear, "_", mncells))

  ggsave(filename = paste0("./output/figures/", spec, "_", lyear, ".png"), g)
  if (printplot == TRUE) plot(g)

}

# Plot incrementing time series with 5 emerging status levels
plot_incr_em5 <- function(df, spec, printplot = FALSE) {

  #spec <- as.character(df[1,"taxonKey"])
  lyear <- max(df$year)
  mncells <- filter(df, year == lyear) %>% .$ncells

  g <- ggplot(df, aes(x = year, y = ncells, colour = em)) + geom_point() +
    geom_line(colour = "grey") +
    scale_colour_manual(values = c("3" = "red", "2" = "orangered", "1" = "orange",
                                   "0" = "grey60", "-1" = "yellow",
                                   "-2" = "green", "-3" = "dark green"),
                                   na.value = "light grey") +
    ggtitle(paste0(spec, "_", lyear, "_", mncells))

  ggsave(filename = paste0("./output/figures/", spec, "_", lyear, ".png"), g)
  if (printplot == TRUE) plot(g)

}


