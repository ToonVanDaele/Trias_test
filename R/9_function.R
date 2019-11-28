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
predict_real_scale <- function(df_n, model) {

  temp <- predict(object = model, newdata = df_n, type = "iterms",
                  interval = "prediction",
                  se.fit = TRUE)

  intercept <- unname(model$coefficients[1])
  df_n$fit <- model$family$linkinv(temp$fit[,1] + intercept)
  df_n$ucl <- model$family$linkinv(temp$fit[,1] + intercept + temp$se.fit[,1] * 1.96)
  df_n$lcl <- model$family$linkinv(temp$fit[,1] + intercept - temp$se.fit[,1] * 1.96)

  return(df_n)
}
