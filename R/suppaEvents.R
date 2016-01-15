#' Parses splicing events of a specific event type from SUPPA
#'
#' @details More information about SUPPA available at
#' \url{https://bitbucket.org/regulatorygenomicsupf/suppa}
#' 
#' @param event Character vector: Splicing event attributes and junction
#' positions
#' @param event_type Character: Alternative splicing event type (see details)
#'
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (exon skipping)}
#'  \item{\bold{RI} (intron retention)}
#'  \item{\bold{MX} (mutually exclusive exons)}
#'  \item{\bold{A5} (alternative 5' splice site)}
#'  \item{\bold{A3} (alternative 3' splice site)}
#'  \item{\bold{AL} (alternative last exon)}
#'  \item{\bold{AF} (alternative first exon)}
#' }
#'
#' @note It only allows to parse one event type at once.
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @export
#'
#' @examples
#' event <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
#' parseSuppaEvent(event)
parseSuppaEvent <- function(event) {
    # Split event ID by semicolon and colon symbols
    event <- stringr::str_split_fixed(event, pattern = ";|:|-", n = 12)
    event <- data.frame(event, stringsAsFactors = FALSE)
    
    # Create list with event attributes
    event_attrs <- data.frame("Program" = "SUPPA",
                              "Gene" = event[[1]],
                              "Chromosome" = event[[3]],
                              stringsAsFactors = FALSE)
    
    event_type <- event[1, 2]
    # Get index of strand (depends on event type)
    strand <- switch(event_type,
                     "SE" = 8, "MX" = 12, "A5" = 8, "A3" = 8,
                     "AF" = 10, "AL" = 10, "RI" = 8)
    
    event_attrs[["Event.type"]] <- switch(event_type,
                                          "SE"="SE",   "MX"="MXE",
                                          "A5"="A5SS", "A3"="A3SS",
                                          "AF"="AFE",  "AL"="ALE",
                                          "RI"="RI")
    
    event_attrs[["Strand"]] <- ifelse(event[[strand]] == "+", "+", "-")
    # Get the junction positions for each exon and parse them
    junctions <- event[4:(strand-1)]
    
    event_type <- event_attrs[["Event.type"]][[1]]
    # Parse junction positions according to event type
    parseJunctions <- switch(event_type,
                             "A3SS" = parseSuppaA3SS,
                             "A5SS" = parseSuppaA5SS,
                             "SE"   = parseSuppaSE,
                             "MXE"  = parseSuppaMXE,
                             "RI"   = parseSuppaRI,
                             "AFE"  = parseSuppaAFE,
                             "ALE"  = parseSuppaALE)
    parsed <- parseJunctions(junctions, event[[strand]])
    return(cbind(event_attrs, parsed))
}

#' Parse junctions of an event from SUPPA according to event type
#'
#' @param junctions List of integers: exon-exon junctions of an event
#' @param strand Character: positive ("+") or negative ("-") strand
#'
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (exon skipping)}
#'  \item{\bold{RI} (intron retention)}
#'  \item{\bold{MXE} (mutually exclusive exons)}
#'  \item{\bold{A5SS} (alternative 5' splice site)}
#'  \item{\bold{A3SS} (alternative 3' splice site)}
#'  \item{\bold{ALE} (alternative last exon)}
#'  \item{\bold{AFE} (alternative first exon)}
#' }
#'
#' @seealso \code{\link{parseSuppaEvent}}
#'
#' @return Data frame of parsed junctions
#' @export
#'
#' @examples
#' junctions <- c(169768099, 169770024, 169770112, 169771762)
#' parseSuppaSE(junctions, "+")
parseSuppaSE <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.end",
                   "A1.start", "A1.end",
                   "C2.start")] <- junctions[plus, ]
    # Minus strand
    parsed[!plus, c("C2.start",
                    "A1.end", "A1.start",
                    "C1.end")] <- junctions[!plus, ]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(196709749, 196709922, 196711005, 196711181)
#' parseSuppaRI(junctions, "+")
parseSuppaRI <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "C2.start", "C2.end")] <- junctions[plus, ]
    # Minus strand
    parsed[!plus, c("C2.end", "C2.start",
                    "C1.end", "C1.start")] <- junctions[!plus, ]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(24790610, 24792494, 24792800, 24790610, 24795476, 24795797)
#' parseSuppaALE(junctions, "+")
parseSuppaALE <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C2.start", "C2.end",
                   "C1.end",
                   "A1.start", "A1.end")] <- junctions[plus, 2:6]
    # Minus strand
    parsed[!plus, c("C2.start", "C2.end",
                    "C1.end",
                    "A1.start", "A1.end")] <- junctions[!plus, 5:1]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169763871, 169764046, 169767998, 169764550, 169765124,
#'                169767998)
#' parseSuppaAFE(junctions, "+")
parseSuppaAFE <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "C2.start",
                   "A1.start", "A1.end")] <- junctions[plus, 1:5]
    # Minus strand
    parsed[!plus, c("C1.start", "C1.end",
                    "C2.start",
                    "A1.start", "A1.end")] <- junctions[!plus, 6:2]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(202060671, 202068453, 202068489, 202073793, 202060671, 
#'                202072798, 202072906, 202073793)
#' parseSuppaMXE(junctions, "+")
parseSuppaMXE <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Note that inclusion values are related to the first alternative isoform
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.end",
                   "A1.start", "A1.end",
                   "A2.start", "A2.end",
                   "C2.start")] <- junctions[plus, -c(4:5)]
    # Minus strand
    parsed[!plus, c("C2.start",
                    "A1.end", "A1.start",
                    "A2.end", "A2.start",
                    "C1.end")] <- junctions[!plus, -c(4:5)]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169772450, 169773216, 169772450, 169773253)
#' parseSuppaA3SS(junctions, "+")
parseSuppaA3SS <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Note that inclusion values are related to the first alternative isoform
    plus <- strand == "+"
    # Plus strand
    parsed[plus, "C1.end"]     <- junctions[plus, 1]
    parsed[plus, ][["C2.start"]] <- listJunctions(junctions, plus, c(2, 4))
    # Minus strand
    parsed[!plus, "C1.end"]       <- junctions[!plus, 4]
    parsed[!plus, ][["C2.start"]] <- listJunctions(junctions, !plus, c(1, 3))
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(99890743, 99891188, 99890743, 99891605)
#' parseSuppaA5SS(junctions, "+")
parseSuppaA5SS <- function (junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    # Note that inclusion values are related to the first alternative isoform
    plus <- strand == "+"
    # Plus strand
    parsed[plus, "C2.start"]   <- junctions[plus, 2]
    parsed[plus, ][["C1.end"]] <- listJunctions(junctions, plus, c(1, 3))
    # Minus strand
    parsed[!plus, "C2.start"]   <- junctions[!plus, 3]
    parsed[!plus, ][["C1.end"]] <- listJunctions(junctions, !plus, c(2, 4))
    return(parsed)
}