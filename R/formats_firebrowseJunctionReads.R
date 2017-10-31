firebrowseJunctionReadsFormat <- function() {
    list(
        tablename   = "Junction quantification",   # Name of the created table
        filename    = "junction_quantification",   # Name of the file
        description = "Read counts of splicing junctions",
        dataType    = "Junction quantification",
        matchName   = TRUE, # Should the file name be matched?
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 2,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE, # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 2,    # Index of row/column to check the format
        
        # File string to check
        check = c("junction", "raw_counts", "raw_counts", "raw_counts", 
                  "raw_counts"),
        
        # Parsing information
        delim       = "\t",  # Delimiter used to separate fields
        colNames    = 1,     # Row to use for column names
        rowNames    = 1,     # Column to use for row names
        ignoreCols  = 1,     # Columns to ignore
        ignoreRows  = NULL,  # Rows to ignore
        commentChar = NULL,  # Ignore lines starting with this string
        
        # Other options
        unique = TRUE,    # Remove duplicated rows
        
        # Identity of rows and columns
        rows    = "splice junctions",
        columns = "samples",
        
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

attr(firebrowseJunctionReadsFormat, "loader") <- "formats"