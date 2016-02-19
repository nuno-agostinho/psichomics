tablename <- "Junction quantification" # Name of the created table
filename <- "junction_quantification" # Name of the file
delim <- "\t"                   # Delimiter used to separate fields
transpose <- FALSE              # Tranpose data
headerCheck <- 2                # Header row in the file to check format
headerUse <- 1                  # Header row in the file to use as header
commentChar <- NULL # String to identify comments (these lines will be ignored)

# First line of data.frame with content of interest
# If you have a header on line 1, the content probably starts on line 2
contentStart <- 3

# Header of the file
header <- c("junction", "raw_counts", "raw_counts", "raw_counts", "raw_counts")