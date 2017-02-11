gtexJunctionReadsFormat <- function() {
    list(
        tablename   = "Junction quantification", # Name of the created table
        filename    = "GTEx_Analysis_v6_RNA-seq_Flux1.6_junction_reads.txt",
        description = "Read counts of splicing junctions",
        dataType    = "Junction quantification", # General category for the data
        
        # Transpose the data? This is the first step before parsing the information!
        # After transposition, a row of the current data equals a column of the original
        skip        = 1,     # Rows to skip when parsing file
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of the row or column used to check the format
        
        # File string to check
        check = c("TargetID", "Gene_Symbol", "Chr", "Coord"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1:4,    # Columns to ignore
        ignoreRows  = 1,    # Rows to ignore
        commentChar = NULL, # String to identify comments (these lines will be ignored)
        
        # Other options
        unique = TRUE,   # Remove duplicated rows
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Clean gtf filenames of columns
            colnames(data) <- gsub("\\.gtf$", "", colnames(data))
            
            # Transform junction positions for easier parsing
            cols <- str_split_fixed(rownames(data), "_", 3)
            rownames(data) <- sprintf("chr%s:%s:%s",
                                      cols[ , 1], cols[ , 2], cols[ , 3])
            return(data)
        }
    )
}

attr(gtexJunctionReadsFormat, "loader") <- "formats"