firehoseJunctionReadsFormat <- function() {
    list(
        tablename   = "Junction quantification",   # Name of the created table
        filename    = "junction_quantification",   # Name of the file
        description = "Read counts of splicing junctions",
        dataType    = "Junction quantification",
        matchName   = TRUE, # Should the file name be matched?
        
        # Transpose the data? This is the first step before parsing the information!
        # After transposition, a row of the current data equals a column of the original
        skip        = 2,     # Rows to skip when parsing file
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using a row (TRUE) or a column (FALSE)
        checkIndex  = 2,     # Index of the row or column used to check the format
        
        # File string to check
        check = c("junction", "raw_counts", "raw_counts", "raw_counts", "raw_counts"),
        
        # Parsing information
        delim       = "\t",  # Delimiter used to separate fields
        colNames    = 1,     # Row to use for column names
        rowNames    = 1,     # Column to use for row names
        ignoreCols  = 1,     # Columns to ignore
        ignoreRows  = NULL,  # Rows to ignore
        commentChar = NULL,  # String to identify comments (these lines will be ignored)
        
        # Other options
        unique = TRUE,    # Remove duplicated rows
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Transform junction positions for easier parsing
            cols <- str_split_fixed(rownames(data), ":|,", 6)
            rownames(data) <- sprintf(
                "%s:%s:%s:%s",
                cols[ , 1], cols[ , 2], cols[ , 5], cols[ , 3])
            return(data)
        }
    )
}

attr(firehoseJunctionReadsFormat, "loader") <- "formats"