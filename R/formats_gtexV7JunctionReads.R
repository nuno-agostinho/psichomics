gtexV7JunctionReadsFormat <- function() {
    list(
        tablename   = "Junction quantification",
        filename    = "GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct",
        description = "Read counts of splicing junctions",
        dataType    = "Junction quantification",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 3,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE, # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,    # Index of row/column to check the format
        
        # File string to check
        check = c("junction_id", "Description"),
        
        # Parsing information
        delim       = "\t",    # Delimiter used to separate fields
        colNames    = 1,       # Row to use for column names
        rowNames    = 1,       # Column to use for row names
        ignoreCols  = seq(2),  # Columns to ignore
        ignoreRows  = 1,       # Rows to ignore
        commentChar = NULL,    # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = TRUE,
        
        # Identity of rows and columns
        rows    = "splice junctions",
        columns = "samples",
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Clean gtf filenames of columns
            colnames(data) <- gsub("\\.gtf$", "", colnames(data))
            
            # Transform junction positions for easier parsing
            cols <- str_split_fixed(rownames(data), "_", 3)
            rownames(data) <- sprintf(
                "chr%s:%s:%s", cols[ , 1], 
                sprintf("%i", as.numeric(cols[ , 2]) - 1),
                sprintf("%i", as.numeric(cols[ , 3]) + 1))
            return(data)
        }
    )
}

attr(gtexV7JunctionReadsFormat, "loader") <- "formats"
