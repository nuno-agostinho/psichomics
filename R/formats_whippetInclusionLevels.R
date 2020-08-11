#' @importFrom stringr str_match
whippetInclusionLevelsFormat <- function() {
    list(
        tablename   = "Inclusion levels",
        description = "Whippet's PSI values per alternative splicing event",
        dataType    = "Inclusion levels",

        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,

        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format

        # File string to check
        check = c("Gene", "Node", "Coord", "Strand", "Type", "Psi"),

        # Parsing information
        delim       = "\t",   # Delimiter used to separate fields
        colNames    = 1,      # Row to use for column names
        rowNames    = NULL,   # Column to use for row names
        ignoreCols  = NULL,   # Columns to ignore
        ignoreRows  = NULL,   # Rows to ignore
        commentChar = NULL,   # Ignore lines starting with this string

        # Remove duplicated rows
        unique = FALSE,

        # Identity of rows and columns
        rows    = "alternative splicing events",
        columns = "samples",

        process = function(data) {
            # Avoid missing PSI values
            data <- data[!is.na(data$Psi), ]
            if (nrow(data) == 0) return(NULL)

            # Convert factor columns to character
            for (col in which(sapply(data, class) == "factor")) {
                data[ , col] <- as.character(data[ , col])
            }

            coordRegex <- "(chr)*(.*):(.*)-(.*)"
            coordInfo  <- str_match(data$Coord, coordRegex)[ , 3:5]

            chrom <- coordInfo[ , 1]
            start <- coordInfo[ , 2]
            end   <- coordInfo[ , 3]
            id    <- paste(data$Type, chrom, data$Strand, start, end, data$Gene,
                           sep="_")
            # Prepare event data
            rowData <- cbind(id                =id,
                             gene              =data$Gene,
                             "full coordinates"=data$Coord,
                             source            ="whippet",
                             type              =data$Type,
                             subtype           =data$Type,
                             chrom             =chrom,
                             start             =start,
                             end               =end,
                             pos               =NA,
                             strand            =data$Strand,
                             constitutive1     =NA,
                             alternative1      =NA,
                             alternative2      =NA,
                             constitutive2     =NA,
                             data[7:14])
            rownames(rowData) <- rowData$id
            class(rowData) <- c("eventData", class(rowData))

            # Set row names
            psi <- data[ , "Psi", drop=FALSE]
            rownames(psi) <- id
            attr(psi, "rowData") <- rowData
            psi <- preserveAttributes(psi)

            # Scale data if PSI values are between 0 and 100
            maximum <- max(psi, na.rm=TRUE)
            if (maximum > 1 && maximum <= 100) {
                attrs <- attributes(psi)
                psi <- psi / 100
                attributes(psi) <- attrs
            }
            return(psi)
        },

        # Join multiple datasets into a single one
        join = function(data) {
            # Identify identical row names among datasets
            allRownames <- unique(unlist(lapply(data, rownames)))

            # Join different datasets
            prep  <- lapply(data, function(k, rows) k[rows, ], rows=allRownames)
            joint <- do.call(cbind, prep)
            rownames(joint) <- allRownames
            df <- data.frame(joint)

            # Name columns based on file names
            samples <- basename(sapply(data, function(i) attr(i, "filename")))
            colnames(df) <- gsub("\\.psi.*", "", samples)

            # Inherit data attributes
            df <- inheritAttrs(df, data[[1]], inheritColnames=FALSE)
            return(df)
        }
    )
}

attr(whippetInclusionLevelsFormat, "loader") <- "formats"
