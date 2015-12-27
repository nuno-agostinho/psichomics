#' Parses splicing event IDs from SUPPA
#' 
#' @details More information about SUPPA available at
#' \url{https://bitbucket.org/regulatorygenomicsupf/suppa}
#' 
#' @param event Character: Splicing event ID
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @export
#'
#' @examples
#' event <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
#' parseSuppaEventID(event)
#' 
#' events <- c("ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-",
#'             "ENSG00000000419;A3:20:49557492-49558568:49557470-49558568:-",
#'             "ENSG00000000003;A5:X:99890743-99891188:99890743-99891605:-")
#' parseSuppaEventID(events)
parseSuppaEventID <- function(event) {
    # Split event ID by semicolon and colon symbols
    event <- strsplit(event, ";|:")
    # Create a list of lists containing event information
    return(lapply(event, parseSuppaEvent))
}

#' Parses splicing event from SUPPA
#'
#' @details More information about SUPPA available at
#' \url{https://bitbucket.org/regulatorygenomicsupf/suppa}
#' 
#' @param event Character vector: Splicing event attributes and junction
#' positions
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @export
#'
#' @examples
#' event <- c("ENSG00000000419", "A3", "20",
#'            "49557492-49557642", "49557470-49557642", "-")
#' parseSuppaEvent(event)
parseSuppaEvent <- function(event) {
    len <- length(event)
    # Create list with event attributes
    event_attrs <- list("Program" = "SUPPA",
                        "Gene" = event[1],
                        "Chromosome" = event[3],
                        "Strand" = event[len])
    
    type <- event[2]
    event_attrs[["Event type"]] <- switch(type,
                                          "MX" = "MXE",
                                          "A5" = "A5SS",
                                          "A3" = "A3SS",
                                          "AF" = "AFE",
                                          "AL" = "ALE",
                                          type)
    
    # Get the junction positions for each exon and parse them
    junctions <- event[4:(len-1)]
    parsed_junctions <- parseSuppaJunctions(event_attrs[["Event type"]],
                                            event_attrs[["Strand"]],
                                            junctions)
    return(c(event_attrs, parsed_junctions))
}

#' Parses splicing junctions from SUPPA
#' 
#' @param event_type Character: Type of the splicing event
#' @param strand Character: Strand (+ or -)
#' @param junctions Character vector: Splicing junctions 
#'
#' @return List of parsed junctions
#' @export
#' 
#' @note In case the -b V (Variable) option is selected, some variability is
#' allowed in some of the boundaries. This is not accounted at the moment.
#'
#' @examples
#' junctionsA5 <- c("99890743-99891188", "99890743-99891605")
#' parseSuppaJunctions(event_type = "A5", strand = "-", junctions = junctionsA5)
parseSuppaJunctions <- function(event_type, strand, junctions) {
    # Split junctions by the hyphen
    junctions <- strsplit(junctions, "-")
    junctions <- as.numeric(unlist(junctions))
    
    # If minus strand, reverse junctions
    if(strand == "-") junctions <- rev(junctions)
    
    # Fill list of parsed junctions with NAs
    parsed = list("C1 start" = NA, "C1 end" = NA,
                  "A1 start" = NA, "A1 end" = NA,
                  "A2 start" = NA, "A2 end" = NA,
                  "C2 start" = NA, "C2 end" = NA)
    
    # Parse junction positions according to event type
    parseJunctions <- switch(event_type,
                             "A3SS" = parseSuppaA3SS,
                             "A5SS" = parseSuppaA5SS,
                             "SE"   = parseSuppaSE,
                             "MXE"  = parseSuppaMXE,
                             "RI"   = parseSuppaRI,
                             "AFE"  = parseSuppaAFE,
                             "ALE"  = parseSuppaALE)
    parsed <- parseJunctions(junctions, strand, parsed)
    return(parsed)
}

#' Parse junctions of an event from SUPPA according to event type
#'
#' @param junctions List of integers: exon-exon junctions of an event
#' @param strand Character: positive ("+") or negative ("-") strand
#' @param parsed Named list filled with NAs for faster execution (optional)
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
#' @return List of parsed junctions
#' @export
#'
#' @examples
#' junctions <- c(169768099, 169770024, 169770112, 169771762)
#' parseSuppaSE(junctions, "+")
parseSuppaSE <- function (junctions, strand, parsed=list()) {
    parsed[c("C1 end",
             "A1 start", "A1 end",
             "C2 start")] <- junctions
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(196709749, 196709922, 196711005, 196711181)
#' parseSuppaRI(junctions, "+")
parseSuppaRI <- function (junctions, strand, parsed=list()) {
    parsed[c("C1 start", "C1 end",
             "C2 start", "C2 end")] <- junctions
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(24790610, 24792494, 24792800, 24790610, 24795476, 24795797)
#' parseSuppaALE(junctions, "+")
parseSuppaALE <- function (junctions, strand, parsed=list()) {
    parsed[c("C2 start", "C2 end",
             "C1 end",
             "A1 start", "A1 end")] <- junctions[2:6]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169763871, 169764046, 169767998, 169764550, 169765124,
#'                169767998)
#' parseSuppaAFE(junctions, "+")
parseSuppaAFE <- function (junctions, strand, parsed=list()) {
    parsed[c("C1 start", "C1 end",
             "C2 start",
             "A1 start", "A1 end")] <- junctions[1:5]
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(202060671, 202068453, 202068489, 202073793, 202060671, 
#'                202072798, 202072906, 202073793)
#' parseSuppaMXE(junctions, "+")
parseSuppaMXE <- function (junctions, strand, parsed=list()) {
    # PSI value is related to the first alternative isoform
    if (strand == "+") {
        parsed[c("C1 end",
                 "A1 start", "A1 end",
                 "A2 start", "A2 end",
                 "C2 start")] <- junctions[-c(4,5)]
    } else if (strand == "-"){
        parsed[c("C1 end",
                 "A2 start", "A2 end",
                 "A1 start", "A1 end",
                 "C2 start")] <- junctions[-c(4,5)]
    }
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169772450, 169773216, 169772450, 169773253)
#' parseSuppaA3SS(junctions, "+")
parseSuppaA3SS <- function (junctions, strand, parsed=list()) {
    parsed[["C1 end"]] <- junctions[1]
    # PSI value is related to the first alternative isoform
    if (strand == "+") {
        parsed[["C2 start"]] <- junctions[c(2, 4)]
    } else if (strand == "-") {
        parsed[["C2 start"]] <- junctions[c(4, 2)]
    }
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(99890743, 99891188, 99890743, 99891605)
#' parseSuppaA5SS(junctions, "+")
parseSuppaA5SS <- function (junctions, strand, parsed=list()) {
    # PSI value is related to the first alternative isoform
    if (strand == "+") {
        parsed[["C1 end"]]   <- junctions[c(1, 3)]
    } else if (strand == "-") {
        parsed[["C1 end"]]   <- junctions[c(3, 1)]
    }
    parsed[["C2 start"]] <- junctions[2]
    return(parsed)
}