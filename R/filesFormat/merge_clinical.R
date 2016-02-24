tablename   <- "Clinical data"      # Name of the created table
filename    <- "clin.merged.txt"    # Name of the file
description <- "Clinical data of the patients"

# Transpose the data? This is the first step before parsing the information!
# After transposition, a row of the current data equals a column of the original
transpose   <- FALSE

# Format checker information
rowCheck    <- FALSE # Check format using a row (TRUE) or a column (FALSE)
checkIndex  <- 1     # Index of the row or column used to check the format

# Parsing information
delim       <- "\t"  # Delimiter used to separate fields
colNames    <- NULL  # Row to use for column names
rowNames    <- 1     # Column to use for row names
ignoreCols  <- 1     # Columns to ignore
ignoreRows  <- NULL  # Rows to ignore
commentChar <- NULL  # String to identify comments (these lines will be ignored)

# File string to check
check <- c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
           "admin.disease_code", "admin.file_uuid",
           "admin.month_of_dcc_upload")