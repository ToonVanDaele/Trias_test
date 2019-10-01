## TRIAS - Load data from Github and store as RDS file

# Get the data 'cube_belgium'
df <- read.table(file = "./data/cube_belgium.tsv",
                 header = TRUE, stringsAsFactors = FALSE)

df$taxonKey <- as.character(df$taxonKey)
df <- select(df, -min_coord_uncertainty)

saveRDS(object = df, file = "./data/cube_belgium.RDS")


# Get data 'cube_belgium_baseline'
df_bl <- read.table(file = "./data/cube_belgium_baseline.tsv",
                    header = TRUE, stringsAsFactors =  FALSE)

df_bl$classKey <- as.character(df_bl$classKey)
df_bl <- select(df_bl, -min_coord_uncertainty)

saveRDS(object = df_bl, file = "./data/cube_belgium_baseline.RDS")


# Get data 'eea_cell_code_protected_areas'
df_xy <- read.table(file = "./data/intersect_EEA_ref_grid_protected_areas.tsv",
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df_xy <- df_xy %>%
  dplyr::select(eea_cell_code = CELLCODE, x = EOFORIGIN, y = NOFORIGIN,
                natura2000, spa, habitat) %>%
  semi_join(df %>%
               dplyr::select(eea_cell_code),
             by = "eea_cell_code")

saveRDS(df_xy, file = "./data/df_xy.RDS")


# species - kingdomkey - canonicalName
speclist <- unique(df$taxonKey)
library(rgbif)
spec_names <- data.frame(taxonKey = speclist,
                         spn = map_chr(speclist, ~ name_usage(.)$data$canonicalName),
                         kingdomKey = map_chr(speclist, ~ name_usage(.)$data$kingdomKey),
                         classKey = map_chr(speclist, ~ name_usage(.)$data$classKey),
                         stringsAsFactors = FALSE)

# write species list with name (canonical), taxon-, kingdom- and classKey
saveRDS(object = spec_names, file = "./data/spec_names.RDS")
