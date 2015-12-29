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
    parsed_junctions <- parseSuppaJunctions(event_attrs[["Event.type"]][[1]],
                                            event[[strand]], junctions)
    return(cbind(event_attrs, parsed_junctions))
}

#' Parses splicing junctions from SUPPA
#' 
#' @param event_type Character: Type of the splicing event (see details)
#' @param strand Character: Strand (+ or -)
#' @param junctions Character vector: Splicing junctions 
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
#' @return List of parsed junctions
#' @export
#' 
#' @note In case the -b V (Variable) option is selected, some variability is
#' allowed in some of the boundaries. This is not accounted at the moment.
#'
#' @examples
#' junctions <- c("99890743-99891188", "99890743-99891605")
#' parseSuppaJunctions(event_type = "A5SS", strand = "-", junctions = junctions)
parseSuppaJunctions <- function(event_type, strand, junctions) {
    # Fill list of parsed junctions with NAs
    parsed <- as.data.frame(matrix(NA, nrow = nrow(junctions), ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")
    
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
parseSuppaSE <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C1.end",
                         "A1.start", "A1.end",
                         "C2.start")] <- junctions},
            { # if strand is minus
                parsed[c("C2.start",
                         "A1.end", "A1.start",
                         "C1.end")] <- junctions})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(196709749, 196709922, 196711005, 196711181)
#' parseSuppaRI(junctions, "+")
parseSuppaRI <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C1.start", "C1.end",
                         "C2.start", "C2.end")] <- junctions},
            { # if strand is minus
                parsed[c("C2.end", "C2.start",
                         "C1.end", "C1.start")] <- junctions})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(24790610, 24792494, 24792800, 24790610, 24795476, 24795797)
#' parseSuppaALE(junctions, "+")
parseSuppaALE <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C2.start", "C2.end",
                         "C1.end",
                         "A1.start", "A1.end")] <- junctions[2:6]},
            { # if strand is minus
                parsed[c("C2.start", "C2.end",
                         "C1.end",
                         "A1.start", "A1.end")] <- junctions[5:1]})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169763871, 169764046, 169767998, 169764550, 169765124,
#'                169767998)
#' parseSuppaAFE(junctions, "+")
parseSuppaAFE <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C1.start", "C1.end",
                         "C2.start",
                         "A1.start", "A1.end")] <- junctions[1:5]},
            { # if strand is minus
                parsed[c("C1.start", "C1.end",
                         "C2.start",
                         "A1.start", "A1.end")] <- junctions[6:2]})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(202060671, 202068453, 202068489, 202073793, 202060671, 
#'                202072798, 202072906, 202073793)
#' parseSuppaMXE(junctions, "+")
parseSuppaMXE <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    # PSI value is related to the first alternative isoform
    ifelse (strand == "+",
            { # if strand is plus
                parsed[c("C1.end",
                         "A1.start", "A1.end",
                         "A2.start", "A2.end",
                         "C2.start")] <- junctions[-c(4,5)]},
            { # if strand is minus
                parsed[c("C2.start",
                         "A1.end", "A1.start",
                         "A2.end", "A2.start",
                         "C1.end")] <- junctions[-c(4,5)]})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(169772450, 169773216, 169772450, 169773253)
#' parseSuppaA3SS(junctions, "+")
parseSuppaA3SS <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    # inclusion value is related to the first alternative isoform
    ifelse (strand == "+",
            { # if strand is plus
                parsed[["C1.end"]] <- junctions[[1]]
                parsed[["C2.start"]] <- apply(junctions[c(2, 4)], 1, as.list)},
            { # if strand is minus
                parsed[["C1.end"]] <- junctions[[4]]
                parsed[["C2.start"]] <- apply(junctions[c(1, 3)], 1, as.list)})
    return(parsed)
}

#' @rdname parseSuppaSE
#' @examples 
#' 
#' junctions <- c(99890743, 99891188, 99890743, 99891605)
#' parseSuppaA5SS(junctions, "+")
parseSuppaA5SS <- function (junctions, strand, parsed=data.frame("C1.end"=NA)) {
    # PSI value is related to the first alternative isoform
    ifelse (strand == "+",
            { # if strand is plus
                parsed[["C2.start"]] <- junctions[[2]]
                parsed[["C1.end"]] <- apply(junctions[c(1, 3)], 1, as.list)},
            { # if strand is minus
                parsed[["C2.start"]] <- junctions[[3]]
                parsed[["C1.end"]] <- apply(junctions[c(2, 4)], 1, as.list)})
    return(parsed)
}