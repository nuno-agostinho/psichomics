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
            # Discard rows containing alt, random and Un (unknown chromosome)
            discard <- grepl("random|alt|Un", rownames(data))
            if (sum(discard) > 0) {
                data <- data[!discard, ]
                warning(sprintf(paste(
                    "Discarded %s rows containing 'alt', 'random'",
                    "or 'Un' in chromosome names."), sum(discard)))
            }
            if (is.null(data) || nrow(data) == 0) return(NULL)
            
            # Parse the chromosome and genomic range of the junction
            strand <- all(grepl("[+-]$", head(rownames(data))))
            regex  <- "^.*?([0-9XYZWM]{1,}).*?([0-9]{1,}).*?([0-9]{1,}).*"
            if (strand) {
                # Also parse the strand
                regex <- paste0(regex, "([+-])$")
                rep   <- "chr\\1:\\2:\\3:\\4"
            } else {
                regex <- paste0(regex, "$")
                rep   <- "chr\\1:\\2:\\3"
            }
            rownames(data) <- gsub(regex, rep, rownames(data))
            return(data)
        }
    )
}

attr(genericJunctionReadsFormat, "loader") <- "formats"