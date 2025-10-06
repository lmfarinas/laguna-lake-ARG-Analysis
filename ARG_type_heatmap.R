# Load necessary packages
library(pheatmap)

#### Profiling by type; normalized to cell count
# Import args-oap output file: metadata.txt
args_oap_metadata <- read.table('args-oap_output/metadata.txt', header=TRUE, row.names=1)

# Import args-oap output file: normalized_cell.type.txt
normalized_cell.type <- as.data.frame(t(read.table('args-oap_output/normalized_cell.type.txt', row.names=1, header=TRUE)))
rownames(normalized_cell.type) <- sub("_.*", "", rownames(normalized_cell.type)) # Clean rownames

# Add OVERALL abundance row. OVERALL represents the cumulative abundance of ARG types from all samples normalized to the total cell copy number

overall <- c() # declare empty vector
for (i in 1:ncol(normalized_cell.type)) { # Loop over columns of normalized_cell.type
  value <- sum(args_oap_metadata[["nCell"]] * normalized_cell.type[[i]]) / 
    sum(args_oap_metadata[["nCell"]]) # Compute overall for column i
  overall <- c(overall, value) # Append to vector
}
ARG_type_df <- rbind(normalized_cell.type, OVERALL = overall) # add overall vector as new row

# Convert to matrix
ARG_type_mat <- as.matrix(ARG_type_df)

# Get last row index
overall_row <- nrow(ARG_type_mat)

# Reorder columns by values of last row (descending)
ARG_type_mat_ordered <- ARG_type_mat[, order(ARG_type_mat[overall_row, ], decreasing = TRUE)]

# Create heatmap with dendrogram
pheatmap(
  ARG_type_mat_ordered,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red", "darkred"))(100),
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  fontfamily = "sans"
)

# Export ARG table
library(tibble)
ARG_type_transposed <- rownames_to_column(as.data.frame(t(ARG_type_df)), var = "Type") # for export to xlsx file


