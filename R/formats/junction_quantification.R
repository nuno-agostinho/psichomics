tablename   <- "Junction quantification"   # Name of the created table
filename    <- "junction_quantification"   # Name of the file
description <- "Quantification of alternatively spliced junctions"

# Transpose the data? This is the first step before parsing the information!
# After transposition, a row of the current data equals a column of the original
transpose   <- FALSE

# Format checker information
rowCheck    <- TRUE  # Check format using a row (TRUE) or a column (FALSE)
checkIndex  <- 2     # Index of the row or column used to check the format

# Parsing information
delim       <- "\t"  # Delimiter used to separate fields
colNames    <- 1     # Row to use for column names
rowNames    <- NULL  # Column to use for row names
ignoreCols  <- NULL  # Columns to ignore
ignoreRows  <- 1:2   # Rows to ignore
commentChar <- NULL  # String to identify comments (these lines will be ignored)

# File string to check
check <- c("junction", "raw_counts", "raw_counts", "raw_counts", "raw_counts")
