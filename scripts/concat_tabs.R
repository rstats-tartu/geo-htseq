
# List files in tmp directory
files <- snakemake@input #list.files("output/tmp", pattern="parsed_suppfiles__", full.names=TRUE)

# Import csv files
dfs <- lapply(files, read.csv)

# Add Set variable to non-imported files
ncols <- unlist(lapply(dfs, ncol))
short <- dfs[ncols == 8]
long <- dfs[ncols == 9]
cols <- colnames(long[[1]])
short1 <- lapply(head(short), function(x) {x["Set"] <- NA; x[,cols]})

# Concatenate data frames
res <- do.call(rbind, c(short1, long))

# Write to csv
write.csv(res, snakemake@output[[1]], row.names = FALSE)
