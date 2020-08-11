vasttoolsInclusionLevelsTidyFormat <- function() {
    list(
        tablename   = "Inclusion levels",
        description = "VAST-TOOLS' PSI values per alternative splicing event",
        dataType    = "Inclusion levels",

        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,

        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format

        # File string to check
        check = c("EVENT"),

        # Parsing information
        delim       = "\t",   # Delimiter used to separate fields
        colNames    = 1,      # Row to use for column names
        rowNames    = 2,      # Column to use for row names
        ignoreCols  = 1,      # Columns to ignore
        ignoreRows  = NULL,   # Rows to ignore
        commentChar = NULL,   # Ignore lines starting with this string

        # Remove duplicated rows
        unique = FALSE,

        # Identity of rows and columns
        rows    = "alternative splicing events",
        columns = "samples",

        process = function(data) {
            colnames(data) <- gsub("_[12]{0,1}$", "\\1", colnames(data),
                                   perl=TRUE)

            # Avoid rows containing only missing values
            psi <- data[rowSums(!is.na(data)) > 0, ]

            # Check if PSI values are all numeric
            if (!all(sapply(psi, is.numeric))) return(NULL)

            # Scale data if PSI values are between 0 and 100
            maximum <- max(psi, na.rm=TRUE)
            if (maximum > 1 && maximum <= 100) {
                attrs <- attributes(psi)
                psi <- psi / 100
                attributes(psi) <- attrs
            }
            return(psi)
        }
    )
}

attr(vasttoolsInclusionLevelsTidyFormat, "loader") <- "formats"
