tablename   <- "Clinical data"      # Name of the created table
filename    <- "clin.merged.txt"    # Name of the file
description <- "Clinical data of the patients"

# Transpose the data? This is the first step before parsing the information!
# After transposition, a row of the current data equals a column of the original
transpose   <- TRUE

# Format checker information
rowCheck    <- TRUE  # Check format using a row (TRUE) or a column (FALSE)
checkIndex  <- 1     # Index of the row or column used to check the format

# File string to check
check <- c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
           "admin.disease_code", "admin.file_uuid",
           "admin.month_of_dcc_upload")

# Parsing information
delim       <- "\t"  # Delimiter used to separate fields
colNames    <- 1     # Row to use for column names
rowNames    <- NULL  # Column to use for row names
ignoreCols  <- NULL  # Columns to ignore
ignoreRows  <- 1     # Rows to ignore
commentChar <- NULL  # String to identify comments (these lines will be ignored)

# Other options
unique <- FALSE   # Remove duplicated rows

# Default columns to show in the app (NULL to show all)
show <- c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
          "patient.bcr_patient_barcode")