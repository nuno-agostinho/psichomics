tablename <- "clinical"         # Name of the created table
filename <- "clin.merged.txt"   # Name of the file
delim <- "\t"                   # Delimiter used to separate fields
transpose <- TRUE               # Tranpose data
headerCheck <- 1                # Header row in the file to check format
headerUse <- 1                  # Header row in the file to use as header
commentChar <- NULL # String to identify comments (these lines will be ignored)

# First line of data.frame with content of interest
# If you have a header on line 1, the content probably starts on line 2
contentStart <- 2

# Header of the file
header <- c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
            "admin.disease_code", "admin.file_uuid",
            "admin.month_of_dcc_upload")