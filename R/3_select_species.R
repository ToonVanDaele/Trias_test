# Select species

# The selection of species is mainly for testing purpose
# Finally all species will have to pass the whole procedure

selspec <- function(df_pp){

  # Get list of species
  speclist <- unique(df_pp$taxonKey)

  # Return a sample selection of the species

  set.seed(55)
  #sp <- sample(x = speclist, size = 1)
  #sp <- c("2927305")
  # sp <- c("3053406")
  # sp <- c("2867614")
  # #sp <- c("2882849")
  sp <- sample(x = speclist, size = 20)

  df_sp <- df_pp %>%
    filter(taxonKey %in% sp)

  #df_sp <- df_pp

  return(df_sp)
}




