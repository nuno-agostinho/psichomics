# Process event data based on event type
processEventDataByType <- function(type, eventDataByType) {
    dataType <- eventDataByType[[type]]
    if (type == "RI") { # Retained intron
        C1end   <- sapply(dataType$fullCoords, "[[", 2)
        C2start <- sapply(dataType$fullCoords, "[[", 3)
        exon1   <- lapply(seq(nrow(dataType)), function(i)
            c(dataType$firstCoord[[i]], C1end[[i]]))
        exon2   <- lapply(seq(nrow(dataType)), function(i)
            c(C2start[[i]], dataType$lastCoord[[i]]))
        
        isPosStrand <- dataType$strand == "+"
        dataType$constitutive1 <- ifelse(isPosStrand, exon1, exon2)
        dataType$alternative1  <- NA
        dataType$alternative2  <- NA
        dataType$constitutive2 <- ifelse(isPosStrand, exon2, exon1)
    } else if (type == "SE") { # Skipped exon
        isPosStrand <- dataType$strand == "+"
        dataType$constitutive1 <- ifelse(
            isPosStrand, dataType$firstCoord, dataType$lastCoord)
        dataType$alternative1  <- lapply(
            seq(nrow(dataType)), function(i)
                c(dataType$start[[i]], dataType$end[[i]]))
        dataType$alternative2  <- NA
        dataType$constitutive2 <- ifelse(
            isPosStrand, dataType$lastCoord, dataType$firstCoord)
    }
    return(dataType)
}

vasttoolsInclusionLevelsFormat <- function() {
    list(
        tablename   = "Inclusion levels",
        description = "VAST-TOOLS' PSI values per alternative splicing event.",
        dataType    = "Inclusion levels",
        
        # Transpose data before parsing? If so, a row in the transposed dataset
        # would be a column in the original
        skip        = 1,     # Rows to skip when parsing file (include header)
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 2,    # Column to use for row names
        ignoreCols  = NULL, # Columns to ignore
        ignoreRows  = NULL, # Rows to ignore
        commentChar = NULL, # Ignore lines starting with this string
        
        # Remove duplicated rows
        unique = TRUE,
        
        # Identity of rows and columns
        rows    = "alternative splicing events",
        columns = "samples",
        
        process = function(data) {
            data$EVENT <- NULL
            colnames(data)[seq(5)] <- c("gene", "coordinates", "length",
                                        "full coordinates", "type")
            # Strip "_1/_2" from the sample name if needed
            colnames(data)[seq(6, ncol(data))] <- gsub(
                "_[12](-Q){0,1}$", "\\1", colnames(data)[seq(6, ncol(data))],
                perl=TRUE)
            
            # Prepare AS data information
            rowData <- data[seq(5)]
            for (col in seq(ncol(rowData))) {
                rowData[ , col] <- as.character(rowData[ , col])
            }
            rowData <- cbind(id=rownames(rowData), source="vast-tools", rowData)
            rowData$subtype <- rowData$type
            rowData$chrom   <- NA
            rowData$start   <- NA
            rowData$end     <- NA
            rowData$pos     <- NA
            rowData$strand  <- NA
            
            rowData$constitutive1 <- NA
            rowData$alternative1  <- NA
            rowData$alternative2  <- NA
            rowData$constitutive2 <- NA
            quality <- seq(7, ncol(data), 2)
            rowData <- cbind(rowData, data[quality])
            
            # Prepare PSI data and discard events only containing missing values
            psi <- data[seq(6, ncol(data), 2)]
            attr(psi, "rowData") <- rowData
            psi <- preserveAttributes(psi)
            psi <- psi[rowSums(!is.na(psi)) > 0, ]
            
            # Parse AS data information (after discarding many AS events)
            rowData        <- attr(psi, "rowData")
            coordinates    <- as.character(rowData$coordinates)
            coordRegex     <- "^(chr){0,1}(.*):(.*)-(.*)$"
            rowData$chrom  <- gsub(coordRegex, "\\2", coordinates)
            rowData$start  <- as.numeric(gsub(coordRegex, "\\3", coordinates))
            rowData$end    <- as.numeric(gsub(coordRegex, "\\4", coordinates))
            rowData$pos    <- lapply(seq(nrow(rowData)), function(i, start, end)
                c(start[[i]], end[[i]]), rowData$start, rowData$end)
            
            fullCoords         <- gsub("^.*?:(.*?)(:.*$){0,1}$", "\\1",
                                       rowData$`full coordinates`, perl=TRUE)
            rowData$fullCoords <- strsplit(fullCoords, "[-,+=]")
            rowData$firstCoord <- sapply(rowData$fullCoords, "[[", 1)
            rowData$lastCoord  <- sapply(rowData$fullCoords,
                                         function(i) i[[length(i)]])
            rowData$strand     <- ifelse(
                as.numeric(rowData$firstCoord) < as.numeric(rowData$lastCoord),
                "+", "-")
            
            # Parse coordinates
            SEtypes <- c("S", "C1", "C2", "C3", "ANN", "MIC")
            types <- c(
                setNames(rep("A3SS", 2), c("Alt3", "A_Alt3")),
                setNames(rep("A5SS", 2), c("Alt5", "A_Alt5")),
                setNames(rep("RI", 4), c("IR-C", "IR-S", "A_IR-C", "A_IR-S")),
                setNames(rep("SE", 12),
                         paste0(c("A_", ""), rep(SEtypes, each=2))))
            eventTypes <- types[rowData$subtype]
            rowData$type <- ifelse(!is.na(eventTypes), eventTypes, 
                                   rowData$subtype)
            
            rowDataByType <- split(rowData, rowData$type)
            rowDataByType <- lapply(names(rowDataByType),
                                    processEventDataByType, rowDataByType)
            rowData       <- unsplit(rowDataByType, rowData$type)
            
            rowData$firstCoord <- NULL
            rowData$lastCoord  <- NULL
            rowData$fullCoords <- NULL
            
            class(rowData) <- c("eventData", class(rowData))
            attr(psi, "rowData") <- rowData
            
            # Scale data if PSI values are between 0 and 100
            maximum <- max(psi, na.rm=TRUE)
            if (maximum > 1 && maximum <= 100) {
                attrs <- attributes(psi)
                psi <- psi/100
                attributes(psi) <- attrs
            }
            return(psi)
        }
    )
}

attr(vasttoolsInclusionLevelsFormat, "loader") <- "formats"
