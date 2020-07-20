psichomicsInclusionLevelsFormat <- function() {
    list(
        tablename   = "Inclusion levels", # Name of the created table
        description = "PSI values per alternative splicing event",
        dataType    = "Inclusion levels",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = "Inclusion levels",
        
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
        rows    = "alternative splicing events",
        columns = "samples",
        
        process = genericInclusionLevelsFormat()$process
    )
}

attr(psichomicsInclusionLevelsFormat, "loader") <- "formats"