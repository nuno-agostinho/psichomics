recountSampleFormat <- function() {
    list(
        tablename   = "Sample metadata",
        description = "Sample metadata for a SRA project",
        dataType    = "Sample metadata",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = c("project", "sample", "experiment", "run",
                  "read_count_as_reported_by_sra", "reads_downloaded"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 4,    # Column to use for row names
        ignoreCols  = 1,    # Columns to ignore
        ignoreRows  = NULL, # Rows to ignore
        commentChar = NULL, # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = FALSE,
        
        # Identity of rows and columns
        rows    = "samples",
        columns = "attributes",
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Unfold data characteristics if available
            if (any(!is.na(data$characteristics))) {
                extra <- gsub(": ", "\"=\"", data$characteristics, fixed=TRUE)
                extra <- gsub("c(\"", "data.frame(\"", extra, fixed=TRUE)
                extra <- gsub("\")", "\", check.names=FALSE)", extra, fixed=TRUE)
                extra <- paste0("plyr::rbind.fill(", 
                                paste0(extra, collapse=", "), ")")
                extra <- eval(parse(text=extra))
                data  <- cbind(data, extra)
                data$characteristics <- NULL
            }
            return(data)
        }
    )
}

attr(recountSampleFormat, "loader") <- "formats"