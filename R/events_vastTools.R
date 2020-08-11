#' @rdname parseMisoAnnotation
#' 
#' @param complexEvents Boolean: should complex events in A3SS and A5SS be
#' parsed?
#' 
#' @export
#' @examples 
#' # Load sample files
#' folder <- "extdata/eventsAnnotSample/VASTDB/Hsa/TEMPLATES"
#' vastToolsOutput <- system.file(folder, package="psichomics")
#' 
#' vast <- parseVastToolsAnnotation(vastToolsOutput)
parseVastToolsAnnotation <- function(
    folder,
    types=c("ALT3", "ALT5", "COMBI", "IR", "MERGE3m", "MIC", "EXSK", "MULTI"),
    genome="Hsa",
    complexEvents=FALSE) {
    
    display("Retrieving VAST-TOOLS annotation...")
    typesRegex <- sprintf("(%s)", paste(types, collapse="|"))
    typesFile <- list.files(folder, full.names=TRUE, pattern=sprintf(
        "%s\\.%s\\..*.txt$", genome, typesRegex))
    # Remove first file if available
    typesFile <- grep("\\.1\\.txt", typesFile, invert=TRUE, value=TRUE)
    names(typesFile) <- gsub(sprintf("^.*%s\\.(.*?)\\..*\\.txt$", genome),
                             "\\1", typesFile)
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)
    
    display("Parsing VAST-TOOLS annotation...")
    types <- names(annot)
    skippedExon <- c("COMBI", "MERGE3m", "MIC", "EXSK", "MULTI")
    
    parseEvents <- function(i) {
        type <- types[i]
        display(type)
        a <- annot[[i]]
        if (nrow(a) > 0) {
            parsed <- parseVastToolsEvent(a)
            if (!complexEvents) {
                if (type == "ALT3") {
                    filter <- !is.na(parsed$A1.start)
                    parsed <- parsed[filter, ]
                    parsed <- uniqueBy(parsed, "Chromosome", "Strand",
                                       "C1.end", "A1.start", "C2.start")
                } else if (type == "ALT5") {
                    filter <- !is.na(parsed$A1.end)
                    parsed <- parsed[filter, ]
                    parsed <- uniqueBy(parsed, "Chromosome", "Strand",
                                       "C1.end", "A1.end", "C2.start")
                } else if (type %in% skippedExon) {
                    C1.end   <- vapply(parsed$C1.end,   length, numeric(1)) > 1
                    A1.start <- vapply(parsed$A1.start, length, numeric(1)) > 1
                    A1.end   <- vapply(parsed$A1.end,   length, numeric(1)) > 1
                    C2.start <- vapply(parsed$C2.start, length, numeric(1)) > 1
                    filter   <- !(C1.end | A1.start | A1.end | C2.start)
                    parsed   <- parsed[filter, ]
                }
            }
            return(parsed)
        }
    }
    events <- lapply(seq_along(annot), parseEvents)
    events <- rbind.fill(events)
    events <- unique(events)
    names(events)[match("Gene.symbol", names(events))] <- "Gene"
    
    # Remove duplicated skipped exons from multiple event types
    skipped <- events$Event.type == "SE"
    uniq <- uniqueBy(events[skipped, ], "Chromosome", "Strand", 
                     "C1.end", "A1.start", "A1.end", "C2.start")
    events <- rbind(events[!skipped, ], uniq)
    
    class(events) <- c("ASevents", class(events))
    return(events)
}

#' Parses an alternative splicing event from VAST-TOOLS
#'
#' @details Junctions are parsed from 
#' 
#' @param event Data.frame: VAST-TOOLS event containing gene symbol, event ID,
#' length, junctions coordinates, event type and inclusion levels for both
#' samples
#'
#' @note Only supports to parse one event at a time.
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @keywords internal
#'
#' @examples
#' event <- read.table(text = 
#' "NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2 0 N 0 N"
#' )
#' psichomics:::parseVastToolsEvent(event)
parseVastToolsEvent <- function(event) {
    # Create list with event attributes
    event_attrs <- data.frame("Program" = "VAST-TOOLS",
                              "Gene symbol" = as.character(event[[1]]),
                              "Event ID"    = as.character(event[[2]]),
                              stringsAsFactors = FALSE)
    
    # By default, assumes things may be parsable as an exon skipping
    # TODO (NunoA): make sure this is intended...
    event_type <- as.character(event[1, 6])
    event_type <- switch(event_type,
                         "IR-C" = "RI",   "IR-S" = "RI",
                         "Alt3" = "A3SS", "Alt5" = "A5SS",
                         "SE")
    event_attrs[["Event.type"]] <- event_type
    
    # Split junctions position
    coord <- as.character(event[[5]])
    junctions <- strsplit(coord, ":|,|-|=")
    
    # Split multiple acceptors/donors (separated with +)
    splitJunctions <- function(i) {
        split <- strsplit(i, "+", fixed=TRUE)
        if (length(split) < 4)
            split[[4]] <- character(0)
        
        return(split)
    }
    
    junctions <- lapply(junctions, splitJunctions)
    junctions <- data.matrix(do.call(rbind, junctions))
    
    # Get chromosomes and convert numbers to numeric
    event_attrs[["Chromosome"]] <- junctions[, 1]
    nrowJunctions <- nrow(junctions)
    junctions <- junctions[, 2:ncol(junctions)]
    junctions <- matrix(lapply(junctions, as.numeric), nrow = nrowJunctions)
    
    # Get strand for retained intron
    if (event_type == "RI") {
        len <- nchar(coord)
        strand <- substr(coord, len, len)
        parsed <- parseVastToolsRI(junctions, strand)
    } else {
        parseJunctions <- switch(event_type,
                                 "SE"   = parseVastToolsSE,
                                 "A3SS" = parseVastToolsA3SS,
                                 "A5SS" = parseVastToolsA5SS)
        parsed <- parseJunctions(junctions)
    }
    
    if (ncol(event) > 7) {
        more_attrs <- data.frame("Inclusion level A" = as.numeric(event[[7]]),
                                 "Inclusion level B" = as.numeric(event[[9]]),
                                 stringsAsFactors = FALSE)
        return(cbind(event_attrs, more_attrs, parsed))
    } else {
        return(cbind(event_attrs, parsed))
    }
}

#' Parse junctions of an event from VAST-TOOLS according to event type
#'
#' @param junctions Data.frame or matrix: exon-exon junctions of alternative
#' splicing events (it must have 4 columns)
#'
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (skipped exon)}
#'  \item{\bold{RI} (retained intron)}
#'  \item{\bold{A5SS} (alternative 5' splice site)}
#'  \item{\bold{A3SS} (alternative 3' splice site)}
#' }
#'
#' @seealso \code{\link{parseVastToolsEvent}()}
#'
#' @return List of parsed junctions
#' @keywords internal
#'
#' @examples
#' junctions <- read.table(text = "41040823 41046768 41046903 41051785")
#' psichomics:::parseVastToolsSE(junctions)
#' 
#' # these functions are vectorised!
#' junctions <- read.table(text = "41040823 41046768 41046903 41051785
#'                                 58864658 58864693 58864294 58864563")
#' psichomics:::parseVastToolsSE(junctions)
parseVastToolsSE <- function (junctions) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Strand is plus if the first junction is lower than the last junction
    plus <- sapply(junctions[, 1], "[[", 1) < sapply(junctions[, 4], "[[", 1)
    parsed[["Strand"]] <- ifelse(plus, "+", "-")
    
    parsed[["C1.end"]]   <- junctions[, 1]
    parsed[["C2.start"]] <- junctions[, 4]
    
    # Plus strand
    parsed[plus, ][["A1.start"]] <- junctions[plus, 2]
    parsed[plus, ][["A1.end"]]   <- junctions[plus, 3]
    # Minus strand
    parsed[!plus, ][["A1.start"]] <- junctions[!plus, 3]
    parsed[!plus, ][["A1.end"]]   <- junctions[!plus, 2]
    return(parsed)
}

#' @rdname parseVastToolsSE
#' @param strand Character: positive (+) or negative (-) strand
#'
#' @examples 
#' 
#' junctions <- read.table(text = "58864658 58864693 58864294 58864563")
#' psichomics:::parseVastToolsRI(junctions, strand = "+")
parseVastToolsRI <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    parsed[["Strand"]] <- strand
    # Plus strand
    parsed[plus, ][["C1.start"]] <- junctions[plus, 1]
    parsed[plus, ][["C1.end"]]   <- junctions[plus, 2]
    parsed[plus, ][["C2.start"]] <- junctions[plus, 3]
    parsed[plus, ][["C2.end"]]   <- junctions[plus, 4]
    # Minus strand
    parsed[!plus, ][["C1.start"]] <- junctions[!plus, 2]
    parsed[!plus, ][["C1.end"]]   <- junctions[!plus, 1]
    parsed[!plus, ][["C2.start"]] <- junctions[!plus, 4]
    parsed[!plus, ][["C2.end"]]   <- junctions[!plus, 3]
    return(parsed)
}

#' @rdname parseVastToolsSE
#'
#' @examples 
#' 
#' junctions <- rbind(
#'     c(36276385, list(c(36277798, 36277315)), 36277974),
#'     c(7133604, 7133377, list(c(7133474, 7133456)))
#' )
#' psichomics:::parseVastToolsA3SS(junctions)
parseVastToolsA3SS <- function (junctions) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Check if there aren't junctions missing
    is2Available <- sapply(junctions[,2], length) > 0
    is3Available <- sapply(junctions[,3], length) > 0
    
    # Strand is plus if the first junction is lower than the other junctions
    available <- ifelse(is3Available, junctions[, 3], junctions[, 2])
    plus <- sapply(junctions[, 1], "[[", 1) < sapply(available, "[[", 1)
    parsed[["Strand"]] <- ifelse(plus, "+", "-")
    parsed[["C1.end"]] <- junctions[, 1]
    
    # Plus strand
    plus3 <- plus & is3Available
    bigList <- sapply(junctions[, 2], length) > 2 # filter unrecognised events
    parsed[plus & !bigList, ][c("A1.start", "A2.start")] <-
        ldply(junctions[plus & !bigList, 2])
    parsed[plus & bigList, ][["A2.start"]] <- junctions[plus & bigList, 2]
    parsed[plus3, ][["A2.end"]] <- junctions[plus3, 3]
    
    # Minus strand
    minus2 <- !plus & is2Available
    bigList <- sapply(junctions[, 3], length) > 2 # filter unrecognised events
    parsed[!plus & !bigList, ][c("A1.start", "A2.start")] <-
        ldply(junctions[!plus & !bigList, 3])
    parsed[!plus & bigList, ][["A2.start"]] <- junctions[!plus & bigList, 3]
    parsed[minus2, ][["A2.end"]] <- junctions[minus2, 2]
    return(parsed)
}

#' @rdname parseVastToolsSE
#'
#' @examples 
#' 
#' junctions <- rbind(
#'     c(74650610, list(c(74650654, 74650658)), 74650982),
#'     c(list(c(49557666, 49557642), 49557746, 49557470))
#' )
#' psichomics:::parseVastToolsA5SS(junctions)
parseVastToolsA5SS <- function (junctions) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Check if there aren't junctions missing
    is1Available <- sapply(junctions[,1], length) > 0
    is2Available <- sapply(junctions[,2], length) > 0
    
    # Strand is plus if the first junction is lower than the other junctions
    available <- ifelse(is2Available, junctions[, 2], junctions[, 1])
    plus <- sapply(available, "[[", 1) < sapply(junctions[, 3], "[[", 1)
    parsed[["Strand"]] <- ifelse(plus, "+", "-")
    parsed[["C2.start"]] <- junctions[, 3]
    
    # Plus strand
    plus1 <- plus & is1Available
    bigList <- sapply(junctions[, 2], length) > 2 # filter unrecognised events
    parsed[plus1, ][["A2.start"]] <- junctions[plus1, 1]
    parsed[plus & !bigList, ][c("A1.end", "A2.end")] <-
        ldply(junctions[plus & !bigList, 2])
    parsed[plus & bigList, ][["A2.end"]] <- junctions[plus & bigList, 2]
    
    # Minus strand
    minus2 <- !plus & is2Available
    bigList <- sapply(junctions[, 1], length) > 2 # filter unrecognised events
    parsed[minus2, ][["A2.start"]] <- junctions[minus2, 2]
    parsed[!plus & !bigList, ][c("A1.end", "A2.end")] <-
        ldply(junctions[!plus & !bigList, 1])
    parsed[!plus & bigList, ][["A2.end"]] <- junctions[!plus & bigList, 1]
    return(parsed)
}
