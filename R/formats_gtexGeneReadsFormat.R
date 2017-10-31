gtexGeneReadsFormat <- function() {
    list(
        tablename   = "Gene expression",
        filename    = "GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct",
        description = "Read counts of splicing junctions",
        dataType    = "Gene expression",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 3,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = c("Name", "Description"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1:2,  # Columns to ignore
        ignoreRows  = 1,    # Rows to ignore
        commentChar = NULL, # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = TRUE,
        
        # Identity of rows and columns
        rows    = "genes",
        columns = "samples",
        
        # Default columns to show (NULL to show all)
        show = NULL
    )
}

attr(gtexGeneReadsFormat, "loader") <- "formats"