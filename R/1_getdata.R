## Load data from Github and store as RDS file

# Data 'cube_belgium'
df <- read.table(file = "./data/cube_belgium.tsv",
                 header = TRUE, stringsAsFactors = FALSE)

df$taxonKey <- as.character(df$taxonKey)
df <- select(df, -min_coord_uncertainty)

saveRDS(object = df, file = "./data/cube_belgium.RDS")


# Data 'cube_belgium_baseline'
df_bl <- read.table(file = "./data/cube_belgium_baseline.tsv",
                    header = TRUE, stringsAsFactors =  FALSE)

df_bl$classKey <- as.character(df_bl$classKey)
df_bl <- select(df_bl, -min_coord_uncertainty)

saveRDS(object = df_bl, file = "./data/cube_belgium_baseline.RDS")

# species - kingdomkey - canonicalName
speclist <- unique(df$taxonKey)
library(rgbif)
spec_names <- data.frame(taxonKey = speclist,
                         spn = map_chr(speclist, ~ name_usage(.)$data$canonicalName),
                         kingdomKey = map_chr(speclist, ~ name_usage(.)$data$kingdomKey),
                         classKey = map_chr(speclist, ~ name_usage(.)$data$classKey),
                         stringsAsFactors = FALSE)

saveRDS(object = spec_names, file = "./data/spec_names.RDS")

