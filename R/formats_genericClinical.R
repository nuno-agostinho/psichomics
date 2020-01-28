genericClinicalFormat <- function() {
    list(
        tablename   = "Clinical data",
        # filename    = "GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",
        description = "Clinical data",
        dataType    = "Clinical data",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = "Subject ID",
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1,    # Columns to ignore
        ignoreRows  = NULL, # Rows to ignore
        commentChar = NULL, # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = FALSE,
        
        # Identity of rows and columns
        rows    = "subjects",
        columns = "attributes",
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) return(data)
    )
}

attr(genericClinicalFormat, "loader") <- "formats"