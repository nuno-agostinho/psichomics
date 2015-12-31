#' Parse alternative splicing events from MATS
#'
#' @param event Data frame row: MATS splicing event
#' @param event_type Character: Type of event to parse (see details)
#'
#' @details The following event types can be parsed:
#' \itemize{
#'  \item{\strong{SE}: Skipping exon}
#'  \item{\strong{MXE}: Mutually exclusive exons}
#'  \item{\strong{RI}:Retained intron}
#'  \item{\strong{A3SS}: Alternative 3' splice site}
#'  \item{\strong{A5SS}: Alternative 5' splice site}
#' }
#'
#' @return List containing the event attributes and junctions
#' @export
#'
#' @examples
#' # MATS event (alternative 3' splice site)
#' event <- read.table(text = "
#'      2 ENSG00000166012 TAF1D chr11 - 93466515 93466671 93466515 93466563 93467790 93467826
#'      5 ENSG00000166012 TAF1D chr11 - 93466515 93466671 93466515 93466585 93467790 93467826
#'      6 ENSG00000166012 TAF1D chr11 - 93466515 93466585 93466515 93466563 93467790 93467826
#' ")
#' parseMatsEvent(mats_A3SS, "A3SS")
parseMatsEvent <- function(event, event_type) {
    len <- ncol(event)
    # Create list with event attributes
    event_attrs <- data.frame("Program"     = "MATS",
                              "Gene"        = as.character(event[[2]]),
                              "Gene symbol" = as.character(event[[3]]),
                              "Chromosome"  = as.character(event[[4]]),
                              "Strand"      = as.character(event[[5]]),
                              "Event type"  = event_type, 
                              stringsAsFactors = FALSE)
    
    # Parse junction positions according to event type
    strand <- as.character(event[[5]])
    
    # Parse junction positions according to event type
    parseJunctions <- switch(event_type,
                             "SE"   = parseMatsSE,
                             "MXE"  = parseMatsMXE,
                             "RI"   = parseMatsRI,
                             "A3SS" = parseMatsA3SS,
                             "A5SS" = parseMatsA5SS,
                             "AFE"  = parseMatsAFE,
                             "ALE"  = parseMatsALE)
    junctions <- event[6:length(event)]
    parsed <- parseJunctions(junctions, strand)
    
    # Add additional attributes if they exist
    if (len > 13) {
        more_attrs <- data.frame(
            "P value" = as.numeric(event[[len - 4]]),
            "FDR" = as.numeric(event[[len - 3]]),
            "Inclusion level A" = as.numeric(event[[len - 2]]),
            "Inclusion level B" = as.numeric(event[[len - 1]]),
            stringsAsFactors = FALSE)
        return(cbind(event_attrs, parsed, more_attrs))
    } else {
        return(cbind(event_attrs, parsed))
    }
}

#' Parse junctions of an alternative splicing event from MATS according to event 
#' type
#'
#' @param junctions Vector of integers with the event's junctions
#' @param strand Character: strand of the event
#' 
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (exon skipping)}
#'  \item{\bold{MXE} (mutually exclusive exon)}
#'  \item{\bold{RI} (intron retention)}
#'  \item{\bold{A5SS} (alternative 5' splice site)}
#'  \item{\bold{A3SS} (alternative 3' splice site)}
#'  \item{\bold{AFE} (alternative first exon)}
#'  \item{\bold{ALE} (alternative last exon)}
#' }
#' 
#' @seealso \code{\link{parseMisoEvent}}
#' 
#' @return Data frame with parsed junctions
#' @export
#' 
#' @examples 
#' junctions <- read.table(text="79685787 79685910 79685796 79685910 79679566 79679751")
#' parseMatsSE(junctions, strand = "+")
parseMatsSE <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("A1.start", "A1.end",
                   "C1.start", "C1.end",
                   "C2.start", "C2.end")] <- junctions[plus, ]
    # Minus strand
    parsed[!plus, c("A1.end", "A1.start",
                    "C2.end", "C2.start",
                    "C1.end", "C1.start")] <- junctions[!plus, ]
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="158282161 158282276 158282689 158282804 158281047 158281295 158283950 158284199")
#' parseMatsMXE(junctions, strand = "+")
parseMatsMXE <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("A1.start", "A1.end",
                   "A2.start", "A2.end",
                   "C1.start", "C1.end",
                   "C2.start", "C2.end")] <- junctions[plus, ]
    # Minus strand
    parsed[!plus, c("A1.end", "A1.start",
                    "A2.end", "A2.start",
                    "C2.end", "C2.start",
                    "C1.end", "C1.start")] <- junctions[!plus, ]
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="15929853 15932100 15929853 15930016 15930687 15932100")
#' parseMatsRI(junctions, strand = "+")
parseMatsRI <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "C2.start", "C2.end")] <- junctions[plus, ]
    parsed[plus, c("C1.start", "C1.end",
                   "C2.start", "C2.end")] <- junctions[plus, 3:6]
    # Minus strand
    parsed[!plus, c("C2.end", "C2.start",
                    "C1.end", "C1.start")] <- junctions[!plus, ]
    parsed[!plus, c("C1.start", "C1.end",
                    "C2.start", "C2.end")] <- junctions[!plus, 6:3]
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="79685787 79685910 79685796 79685910 79679566 79679751")
#' parseMatsA3SS(junctions, strand = "+")
parseMatsA3SS <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    # Note that inclusion values are related to the first alternative isoform
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "C2.end")] <- junctions[plus, c(5:6, 2)]
    parsed[plus, ][["C2.start"]] <- apply(junctions[plus, c(1, 3)],
                                        1, numericList)
    # Minus strand
    parsed[!plus, c("C1.start", "C1.end",
                    "C2.end")] <- junctions[!plus, c(6:5, 1)]
    parsed[!plus, ][["C2.start"]] <- apply(junctions[!plus, c(2, 4)],
                                         1, numericList)
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="102884421 102884501 102884421 102884489 102884812 102885881")
#' parseMatsA5SS(junctions, strand = "+")
parseMatsA5SS <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    # Note that inclusion values are related to the first alternative isoform
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start",
                   "C2.start", "C2.end")] <- junctions[plus, c(1, 5:6)]
    parsed[plus, ][["C1.end"]] <- apply(junctions[plus, c(2, 4)],
                                        1, numericList)
    # Minus strand
    parsed[!plus, c("C1.start",
                    "C2.start", "C2.end")] <- junctions[!plus, c(2, 6, 5)]
    parsed[!plus, ][["C1.end"]] <- apply(junctions[!plus, c(1, 3)],
                                           1, numericList)
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="16308723 16308879 16308967 16309119 16314269 16314426")
#' parseMatsAFE(junctions, strand = "+")
parseMatsAFE <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "A1.start", "A1.end",
                   "C2.start", "C2.end")] <- junctions[plus, 1:6]
    # Minus strand
    parsed[!plus, c("C1.start", "C1.end",
                    "A1.start", "A1.end",
                    "C2.start", "C2.end")] <- junctions[!plus, c(2:1, 4:3, 6:5)]
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="111858645 111858828 111851063 111851921 111850441 111850543")
#' parseMatsAFE(junctions, strand = "+")
parseMatsALE <- function(junctions, strand) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createFilledJunctions(nrow(junctions))
    
    plus <- strand == "+"
    # Plus strand
    parsed[plus, c("C1.start", "C1.end",
                   "A1.start", "A1.end",
                   "C2.start", "C2.end")] <- junctions[plus, c(5:6, 1:2, 3:4)]
    # Minus strand
    parsed[!plus, c("C1.start", "C1.end",
                    "A1.start", "A1.end",
                    "C2.start", "C2.end")] <- junctions[!plus, c(6:5, 2:1, 4:3)]
    return(parsed)
}