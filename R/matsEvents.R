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
    parsed <- parseMatsJunctions(event, event_type)
    
    # Add more attributes only if the event has them
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


#' Parse MATS alternative splicing junctions
#'
#' @inheritParams parseMatsEvent
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
#' @return List containing the junctions of the event
#' @export
#'
#' @examples
#' # MATS event (alternative 3' splice site)
#' event <- read.table(text = "
#'      3658 ENSG00000067715 SYT1 chr12 + 79685787 79685910 79685796 79685910 79679566 79679751 3658 252 102 73 16 58 56 0.0342916452301 0.274333161841 0.705 0.815 -0.11
#' ")
#' 
#' parseMatsJunctions(mats_A3SS, "A3SS")
parseMatsJunctions <- function(event, event_type) {
    # Fill list of parsed junctions with NAs
    parsed <- as.data.frame(matrix(NA, nrow = nrow(event), ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")
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
    parsed <- parseJunctions(junctions, strand, parsed)
    return(parsed)
}

#' Parse junctions of an alternative splicing event from MATS according to event 
#' type
#'
#' @param junctions Vector of integers with the event's junctions
#' @param strand Character: strand of the event
#' @param parsed Named data frame filled with NAs for faster execution
#' (optional)
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
parseMatsSE <- function(junctions, strand, parsed=data.frame("A1.start"=NA)) {
    ifelse(strand == "+", 
           { # strand is plus
               parsed[c("A1.start", "A1.end",
                        "C1.start", "C1.end",
                        "C2.start", "C2.end")] <- junctions[1:6]},
           { # strand is minus
               parsed[c("A1.end", "A1.start",
                        "C2.end", "C2.start",
                        "C1.end", "C1.start")] <- junctions[1:6]})
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="158282161 158282276 158282689 158282804 158281047 158281295 158283950 158284199")
#' parseMatsMXE(junctions, strand = "+")
parseMatsMXE <- function(junctions, strand, parsed=data.frame("A1.start"=NA)) {
    ifelse (strand == "+",
            # if strand is plus
            parsed[c("A1.start", "A1.end",
                     "A2.start", "A2.end",
                     "C1.start", "C1.end",
                     "C2.start", "C2.end")] <- junctions[1:8],
            # if strand is minus
            parsed[c("A1.end", "A1.start",
                     "A2.end", "A2.start",
                     "C2.end", "C2.start",
                     "C1.end", "C1.start")] <- junctions[1:8]
    )
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="15929853 15932100 15929853 15930016 15930687 15932100")
#' parseMatsRI(junctions, strand = "+")
parseMatsRI <- function(junctions, strand, parsed=data.frame("C1.start"=NA)) {
    ifelse (strand == "+",
            # if strand is plus
            parsed[c("C1.start", "C1.end",
                     "C2.start", "C2.end")] <- junctions[3:6],
            # if strand is minus
            parsed[c("C1.start", "C1.end",
                     "C2.start", "C2.end")] <- junctions[6:3]
    )
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="79685787 79685910 79685796 79685910 79679566 79679751")
#' parseMatsA3SS(junctions, strand = "+")
parseMatsA3SS <- function(junctions, strand, parsed=data.frame(NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C1.start", "C1.end")] <- junctions[5:6]
                parsed[["C2.start"]] <- apply(junctions[c(1, 3)], 1, as.list)
                parsed[["C2.end"]]   <- junctions[[2]]
            }, { # if strand is minus
                parsed[c("C1.start", "C1.end")] <- junctions[6:5]
                parsed[["C2.start"]] <- apply(junctions[c(2, 4)], 1, as.list)
                parsed[["C2.end"]]   <- junctions[[1]]
            }
    )
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="102884421 102884501 102884421 102884489 102884812 102885881")
#' parseMatsA5SS(junctions, strand = "+")
parseMatsA5SS <- function(junctions, strand, parsed=data.frame("C1.start"=NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[["C1.start"]] <- junctions[[1]]
                parsed[["C1.end"]]   <- apply(junctions[c(2, 4)], 1, as.list)
                parsed[c("C2.start", "C2.end")] <- junctions[5:6]
            }, { # if strand is minus
                parsed[["C1.start"]]  <- junctions[[2]]
                parsed[["C1.end"]]    <- apply(junctions[c(1, 3)], 1, as.list)
                parsed[c("C2.start", "C2.end")] <- junctions[6:5]
            })
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="16308723 16308879 16308967 16309119 16314269 16314426")
#' parseMatsAFE(junctions, strand = "+")
parseMatsAFE <- function(junctions, strand, parsed=data.frame("A1.start"=NA)) {
    ifelse (strand == "+", {
        # if strand is plus
        parsed[c("C1.start", "C1.end")] <- junctions[1:2]
        parsed[c("A1.start", "A1.end")] <- junctions[3:4]
        parsed[c("C2.start", "C2.end")] <- junctions[5:6]
    }, { # if strand is minus
        parsed[c("C1.start", "C1.end")] <- junctions[2:1]
        parsed[c("A1.start", "A1.end")] <- junctions[4:3]
        parsed[c("C2.start", "C2.end")] <- junctions[6:5]
    })
    return(parsed)
}

#' @rdname parseMatsSE
#' @examples 
#' 
#' junctions <- read.table(text="111858645 111858828 111851063 111851921 111850441 111850543")
#' parseMatsAFE(junctions, strand = "+")
parseMatsALE <- function(junctions, strand, parsed=data.frame("A1.start"=NA)) {
    ifelse (strand == "+", {
        # if strand is plus
        parsed[c("C1.start", "C1.end")] <- junctions[5:6]
        parsed[c("A1.start", "A1.end")] <- junctions[1:2]
        parsed[c("C2.start", "C2.end")] <- junctions[3:4]
    }, { # if strand is minus
        parsed[c("C1.start", "C1.end")] <- junctions[6:5]
        parsed[c("A1.start", "A1.end")] <- junctions[2:1]
        parsed[c("C2.start", "C2.end")] <- junctions[4:3]
    })
    return(parsed)
}