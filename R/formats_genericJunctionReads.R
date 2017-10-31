genericJunctionReadsFormat <- function() {
    list(
        tablename   = "Junction quantification",
        # filename    = "",
        description = "Read counts of splicing junctions",
        dataType    = "Junction quantification",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = "Junction ID",
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1,    # Columns to ignore
        ignoreRows  = NULL, # Rows to ignore
        commentChar = NULL, # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = TRUE,
        
        # Identity of rows and columns
        rows    = "splice junctions",
        columns = "samples",
        
        process = function(data) {
            strand <- all(grepl("[+-]$", head(rownames(data))))
            if (strand) {
                # Parse the chromosome, start/end positions and strand of the 
                # junction (by this order)
                rownames(data) <- gsub(
                    "^.*?([0-9XY]{1,}).*?([0-9]{1,}).*?([0-9]{1,}).*([+-])$",
                    "chr\\1:\\2:\\3:\\4", rownames(data))
            } else {
                # Parse the chromosome and start/end positions of the junction 
                # (by this order)
                rownames(data) <- gsub(
                    "^.*?([0-9XY]{1,}).*?([0-9]{1,}).*?([0-9]{1,}).*$",
                    "chr\\1:\\2:\\3", rownames(data))
            }
            return(data)
        }
    )
}

attr(genericJunctionReadsFormat, "loader") <- "formats"