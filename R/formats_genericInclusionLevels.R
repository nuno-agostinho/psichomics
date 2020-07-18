genericInclusionLevelsFormat <- function() {
    list(
        tablename   = "Inclusion levels",
        description = "PSI values per alternative splicing event.",
        dataType    = "Inclusion levels",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = "AS Event ID",
        
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
        
        process = function(data) {
            # Scale data if PSI values are between 0 and 100
            maximum <- max(data, na.rm=TRUE)
            if (maximum > 1 && maximum <= 100) data <- data/100
            
            events <- rownames(data)
            if (!is.null(events)) {
                attr(data, "rowData") <- parseSplicingEvent(events, coords=TRUE)
                data <- preserveAttributes(data)
            }
            return(data)
        }
    )
}

attr(genericInclusionLevelsFormat, "loader") <- "formats"