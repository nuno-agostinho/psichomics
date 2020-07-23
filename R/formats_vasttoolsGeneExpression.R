vasttoolsGeneExpressionFormat <- function() {
    list(
        tablename   = "Gene expression",
        description = "VAST-TOOLS' gene expression",
        dataType    = "Gene expression",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = c("ID", "NAME"),
        
        # Parsing information
        delim       = "\t",   # Delimiter used to separate fields
        colNames    = 1,      # Row to use for column names
        rowNames    = 1,      # Column to use for row names
        ignoreCols  = NULL,   # Columns to ignore
        ignoreRows  = NULL,   # Rows to ignore
        commentChar = NULL,   # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = TRUE,
        
        # Identity of rows and columns
        rows    = "genes",
        columns = "samples",
        
        process = function(data) {
            colnames(data) <- gsub("-cRPKM$", "", colnames(data))
            
            # For now, discard count columns (these could be used in the future)
            countCols <- grepl("-Counts$", colnames(data))
            data <- data[ , !countCols]
            
            colnames(data) <- gsub("_[12]{0,1}$", "\\1", colnames(data),
                                   perl=TRUE)
            # Avoid rows containing only missing values
            data <- data[rowSums(!is.na(data[ , -seq(2)])) > 0, ]
            
            # Use unambiguous gene symbols when possible
            genes <- as.character(data$NAME)
            dups  <- names(table(genes)[table(genes) > 1])
            uniq  <- !is.na(genes) & !genes %in% dups
            rownames(data)[uniq] <- genes[uniq]
            return(data[ , -seq(2)])
        }
    )
}

attr(vasttoolsGeneExpressionFormat, "loader") <- "formats"
