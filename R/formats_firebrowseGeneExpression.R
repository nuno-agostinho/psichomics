firebrowseGeneExpressionFormat <- function() {
    list(
        tablename   = "Gene expression",   # Name of the created table
        filename    = "RSEM_genes",   # Name of the file
        description = "Gene expression from RSEM",
        dataType    = "Gene expression",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 2,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 2,     # Index of row/column to check the format
        
        # File string to check
        check = c("gene_id", "raw_count", "scaled_estimate", "transcript_id",
                  "raw_count"),
        
        # Parsing information
        delim       = "\t",  # Delimiter used to separate fields
        colNames    = 1,     # Row to use for column names
        rowNames    = 1,     # Column to use for row names
        ignoreCols  = 1,     # Columns to ignore
        ignoreRows  = NULL,  # Rows to ignore
        commentChar = NULL,  # Ignore lines starting with this string
        
        # Identity of rows and columns
        rows    = "genes",
        columns = "samples",
        
        # Remove duplicated rows
        unique = FALSE,
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            return(data[ , seq(1, ncol(data), 3)])
        }
    )
}

attr(firebrowseGeneExpressionFormat, "loader") <- "formats"