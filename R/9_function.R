# Functions

# Harmonise levels of emerging status among methods

get_em_levels <- function(){

  df_em_levels <- data.frame(em = c(0, 1, 2, 3),
                   status = c("not emerging", "unclear", "potentially emerging",
                              "emerging"))

}

# Map em_gam out to em
em_gam2em <- function(em_gam, nbyear){

  em <- em_gam %>%
    filter(year >= max(year) - nbyear + 1) %>%
    mutate(em_out = case_when(
      em < 0 ~ 0,
      em == 0 ~ 1,
      em < 3 ~ 2,
      em >= 3 ~ 3)) %>%
    select(year, em_out)

    return(em)
}

# Calculate first and second derivative from smoother 'year'
#
# g gam model
# length length of data frame for prediction
#
# return list with derivatives and emerging evaluation

syear_deriv <- function(g, dfrows){

  deriv1 <- derivatives(g, term = "year", type = "central",
                        order = 1, level = 0.8, n = dfrows, eps = 1e-4)
  deriv2 <- derivatives(g, term = "year", type = "central",
                        order = 2, level = 0.8, n = dfrows, eps = 1e-4)

  em = em_level(deriv1, deriv2)

  return(list(deriv1 = deriv1, deriv2 = deriv2, em = em))

}

# Evaluate emerging based on first and second derivative
#
# df1 result from gratia::derivatives - order 1
# df2 result fomr gratia::derivatives - order 2
#
# Return emerging evaluation

em_level <- function(df1, df2){

  em1 <- df1 %>%
    mutate(em1 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em1) %>%
    mutate(year = round(year, 0))

  em2 <- df2 %>%
    mutate(em2 = case_when(
      .$lower < 0  & .$upper < 0 ~ "-1",
      .$lower < 0  & .$upper > 0 ~ "0",
      .$lower > 0  & .$upper > 0 ~ "1")) %>%
    dplyr::select(year = data, em2) %>%
    mutate(year = round(year, 0))

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


# Apply model and save output

apply_method <- function(df, em_method, n2k = FALSE){

  df <- df %>% group_by(taxonKey)
  taxl <- df %>% group_keys() %>% pull(taxonKey)

  result_list <- df %>%
    group_split() %>%
    map(.f = get(em_method)) %>%
    set_names(taxl)

  filename <- paste0("./output/result_", em_method,
                    ifelse(n2k == TRUE, "_n2k", ""), ".RDS")
  saveRDS(result_list, file = filename)
  return(filename)

}



# Retrieve em information from results

get_em <- function(method_name, path){

  filename = paste0(path, method_name, ".RDS")
  cat(filename, "\n")
  result <- readRDS(file = filename)

  em_result <- map_dfr(result, "df_em")

  remove(result)

  return(em_result)
}

# Calculate the mean lower confidence level of the first derivative of the
# last three years
get_lcl <- function(df_deriv, nbyear, fam){

  lcl <- df_deriv %>%
    filter(var == "year") %>%
    filter(data > max(data) - nbyear + 1) %>%
    select(year = data, lcl = lower) %>%
    mutate(year = round(year, 0),
           lcl = fam$linkinv(lcl))

  return(lcl)
}

add_spec <- function(df_ts, spec_names){

  df_ts_species <- unique(df_ts$taxonKey)
  speclist <- df_ts_species[!df_ts_species %in% spec_names$taxonKey]

  if (!length(speclist) == 0){
    library(rgbif)
    new_spec_names <- data.frame(taxonKey = speclist,
                             spn = map_chr(speclist, ~ name_usage(.)$data$canonicalName),
                             kingdomKey = map_chr(speclist, ~ name_usage(.)$data$kingdomKey),
                             classKey = map_chr(speclist, ~ name_usage(.)$data$classKey),
                             stringsAsFactors = FALSE)

    spec_names <- rbind(spec_names, new_spec_names)
    # write species list with name (canonical), taxon-, kingdom- and classKey
    saveRDS(object = spec_names, file = "./data/spec_names.RDS")
  }
  return(spec_names)
}

# Apply predict and link inverse to real scale
predict_real_scale <- function(df_p, model) {

  pf <- predict(object = model, newdata = df_p, type = "response",
                  interval = "prediction",
                  se.fit = TRUE)

  df_p$fit <- pf$fit
  df_p$ucl <- pf$fit + pf$se.fit * 1.96
  df_p$lcl <- pf$fit - pf$se.fit * 1.96

  p1 <- predict(object = model, newdata = df_p, type = "iterms",
                  interval = "prediction",
                  se.fit = TRUE)

  intercept <- unname(model$coefficients["(Intercept)"])
  ulnk <- model$family$linkinv
  df_p$s_year_fit <- ulnk(p1$fit[,"s(year)"] + intercept)
  df_p$s_year_ucl <- ulnk(p1$fit[,"s(year)"] + intercept + p1$se.fit[,"s(year)"] * 1.96)
  df_p$s_year_lcl <- ulnk(p1$fit[,"s(year)"] + intercept - p1$se.fit[,"s(year)"] * 1.96)

  return(df_p)
}

## Predict to real scale

predict_real_scale2 <- function(df_pred, model) {



  return(df_pred)
}

# Deriv real scale
deriv_to_real <- function(df_deriv, model){

  df_temp <- data.frame(fit = model$family$linkinv(df_deriv$derivative),
                        ucl = model$family$linkinv(df_deriv$derivative + df_deriv$se * 1.96),
                        lcl = model$family$linkinv(df_deriv$derivative - df_deriv$se * 1.96))

}

# Temporary workaround upscaling to 5x5km cells
aggr_1to5 <- function(df){

  df_s5 <- df %>%
    left_join(df_xy %>%
                select(eea_cell_code, x5, y5, cell_code5),
              by = "eea_cell_code") %>%
    filter(!is.na(cell_code5)) %>%
    group_by(taxonKey, year, cell_code5) %>%
    summarise(x = first(x5),
              y = first(y5),
              native_obs = sum(native_obs),
              obs = sum(obs),
              pa_native_obs = max(pa_native_obs),
              pa_obs = max(pa_obs),
              classKey = first(classKey),
              n2k = sum(natura2000),
              nb_1km = n()) %>%
    mutate(natura2000 = ifelse(n2k / nb_1km >= 0.5, TRUE, FALSE)) %>%
    select(-n2k, -nb_1km)
}

## Set title for graph of (used when ptitle = NULL)

set_title <- function(df){

  spec <- df[[1, "taxonKey"]]
  spn <- spec_names %>% filter(taxonKey == spec) %>% pull(spn)
  ptitle <- paste0(spec, "_", spn)
  return(ptitle)
}

# Generate multiple plots for each species in a directory output/species

output_multiple_plots <- function(result_sp){

  if (!is.null(result_sp$df_n)){
    df <- result_sp$df_n
    # Create dir if necessary
    spn_f <- set_title(df)
    print(spn_f)
    dir <- paste0("./output/plots/", spn_f, "/")
    dir.create(dir, showWarnings = FALSE)

    # Plot raw data
    ptitle <- paste0(spn_f, "_raw_data")
    pr <- plot_ts(df = df, y_axis = "obs", ptitle = ptitle)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = pr)

    # Plot raw class observation data
    ptitle <- paste0(spn_f, "_class_obs")
    pc <- plot_ts(df = df, y_axis = "cobs", ptitle = ptitle)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = pc)

    # Plot fit
    ptitle <- paste0(spn_f, "_model_fit")
    pf <- plot_ribbon_em(df_n = df, df = df, y_axis = "obs", ptitle = ptitle)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = pf)

    # Plot smoother s_year
    ptitle <- paste0(spn_f, "_smoother_syear")
    ps <- plot_smoother(df, ptitle)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = ps)

    # Plot first derivative
    ptitle <- paste0(spn_f, "_syear_1st_derivative")
    pd1 <- draw(result_sp$deriv1)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = pd1)

    # Plot second derivative
    ptitle <- paste0(spn_f, "_syear_2nd_derivative")
    pd2 <- draw(result_sp$deriv2)
    ggsave(filename = paste0(dir, ptitle, ".png"), plot = pd2)

    # Plot all 4 plots together
    ptitle <- paste0(spn_f, "_plots_arranged")
    p_ga <- arrangeGrob(pr, pc, pf, ps, pd1, pd2, ncol = 2)
    ggsave(filename = paste0(dir, ptitle, ".png"), width = 10, height = 12, plot = p_ga)
  }
  return(spn_f)
}
