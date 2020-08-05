vasttoolsGeneExpressionFormat <- function() {
    list(
        tablename   = "Gene expression (cRPKM)",
        description = "VAST-TOOLS' gene expression (cRPKM)",
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
            # Check if data contains gene read counts
            countCols <- which(grepl("-Counts$", colnames(data)))
            countData <- data[ , c(seq(2), countCols)]
            hasReadCounts <- length(countCols) > 0
            if (hasReadCounts) {
                cRPKMdata <- data[ , -countCols]
            } else {
                cRPKMdata <- data
            }

            prepareData <- function(data) {
                colnames(data) <- gsub("(-cRPKM)|(-Counts$)$", "",
                                       colnames(data))
                colnames(data) <- gsub("_[12]{0,1}$", "\\1", colnames(data),
                                       perl=TRUE)
                # Avoid rows containing only missing values
                data <- data[rowSums(!is.na(data[ , -seq(2)])) > 0, ]

                # Use unambiguous gene symbols when possible
                genes <- as.character(data$NAME)
                dups  <- names(table(genes)[table(genes) > 1])
                uniq  <- !is.na(genes) & !genes %in% dups
                rownames(data)[uniq] <- genes[uniq]
                data <- data[ , -seq(2)]
                return(data)
            }
            cRPKMdata <- prepareData(cRPKMdata)
            cRPKMdata <- addObjectAttrs(
                cRPKMdata,
                "tablename"="Gene expression (cRPKM)",
                "description"="VAST-TOOLS' gene expression (cRPKM)")

            if (hasReadCounts) {
                countData <- prepareData(countData)
                countData <- addObjectAttrs(
                    countData,
                    "tablename"="Gene expression (read counts)",
                    "description"="VAST-TOOLS' gene expression (read counts)")
                res <- list(countData, cRPKMdata)
            } else {
                res <- cRPKMdata
            }
            return(res)
        }
    )
}

attr(vasttoolsGeneExpressionFormat, "loader") <- "formats"
