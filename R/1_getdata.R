## Load data from Github and store as RDS file

# Data from the data-cube
df <- read.table(file = "./data/cube_belgium.tsv", header = TRUE)

df$taxonKey <- as.character(df$taxonKey)

head(df)
str(df)
summary(df)

saveRDS(object = df, file = "./data/cube_belgium.RDS")

