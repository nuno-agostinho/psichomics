tablename   <- "RSEM gene expression"   # Name of the created table
filename    <- "RSEM_genes_normalized"   # Name of the file
description <- "Normalised gene expression from RSEM"

# Transpose the data? This is the first step before parsing the information!
# After transposition, a row of the current data equals a column of the original
skip        <- 2     # Rows to skip when parsing file
transpose   <- FALSE

# Format checker information
rowCheck    <- TRUE  # Check format using a row (TRUE) or a column (FALSE)
checkIndex  <- 2     # Index of the row or column used to check the format

# File string to check
check <- c("gene", "raw_counts", "median_length_normalized", "RPKM",
           "raw_counts", "median_length_normalized")

# Parsing information
delim       <- "\t"  # Delimiter used to separate fields
colNames    <- 1     # Row to use for column names
rowNames    <- 1     # Column to use for row names
ignoreCols  <- 1     # Columns to ignore
ignoreRows  <- NULL  # Rows to ignore
commentChar <- NULL  # String to identify comments (these lines will be ignored)

# Other options
unique <- FALSE    # Remove duplicated rows

# Default columns to show (NULL to show all)
show <- NULL
