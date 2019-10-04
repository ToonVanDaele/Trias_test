# TRIAS - Select species

# The selection of species is mainly for testing purpose
# Finally all species will have to pass the whole procedure

selspec <- function(df, specs){

  set.seed(88)

  #df_sel <- df
  #sp <- sample(x = speclist, size = 1)
  #sp <- c("2927305")
  # sp <- c("3053406")
  # sp <- c("2867614")
  # #sp <- c("2882849")
  sp <- sample(x = specs, size = 20)
  #
   df_sel <- df %>%
     filter(taxonKey %in% sp)

  return(df_sel)
}




