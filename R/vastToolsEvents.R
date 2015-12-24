#' Parses multiple alternative splicing events from VAST-TOOLS
#' 
#' Expects the annotation for many events at once and parses each event.
#' 
#' @param events Data.frame: Vast-Tools events containing gene symbol, event ID,
#' length, junctions coordinates, event type and inclusion levels for both 
#' samples
#' @param progress Boolean: show progress?
#' 
#' @seealso \code{\link{parseVastToolsEvent}}
#' 
#' @return List of lists with the event attributes (chromosome, strand, event 
#' type and the position of the exon boundaries)
#' @export
#' events <- read.table(text = "
#' NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2 0 N,N,N,Bn,S@0,0 0 N,N,N,Bn,S@0,0
#' SCYL3   HsaEX0056691    chr1:169839396-169839498        103     chr1:169842837,169839396-169839498,169838180+169838269  S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
#' LAP3    HsaEX0035325    chr4:17587574-17587861  288     chr4:17586759,17587574-17587861,17590442        S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
#' ")
#' parseMultipleVastToolsEvents(events)
parseMultipleVastToolsEvents <- function(events, progress = FALSE) {
    if (progress) pb <- txtProgressBar(min=1, max=nrow(events), style=3)
    parsedEvents <- lapply(1:nrow(events),
                           function(k) {
                               if (progress) setTxtProgressBar(pb, k)
                               parsed <- parseVastToolsEvent(events[k, ])
                               return(parsed)
                           })
    return(parsedEvents)
}

#' Parses an alternative splicing event from VAST-TOOLS
#'
#' @details Junctions are parsed from 
#' 
#' @param event Data.frame row, list or character: Vast-Tools event containing
#' gene symbol, event ID, length, junctions coordinates, event type and 
#' inclusion levels for both samples
#'
#' @note Don't use strings as factors.
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @export
#'
#' @examples
#' event <- c("NFYA", "HsaEX0042823", "chr6:41046768-41046903", "136", 
#'            "chr6:41040823,41046768-41046903,41051785", "C2", "0", 
#'            "N,N,N,Bn,S@0,0", "0", "N,N,N,Bn,S@0,0")
#' parseVastToolsEvent(event)
parseVastToolsEvent <- function(event) {
    # Create list with event attributes
    event_type <- event[[6]]
    event_attrs <- list("Program" = "VAST-TOOLS",
                        "Gene symbol" = as.character(event[[1]]),
                        "Event ID" = as.character(event[[2]]))
    
    if (length(event) > 7) {
        more_attrs <- list("Inclusion level A" = as.numeric(event[[7]]),
                           "Inclusion level B" = as.numeric(event[[9]]))
    } else {
        more_attrs <- list()
    }
    
    # By default, assumes things may be parsable as a skipping exon
    event_type <- switch(event_type,
                         "IR-C" = "RI",
                         "IR-S" = "RI",
                         "Alt3" = "A3SS",
                         "Alt5" = "A5SS",
                         "SE")
    
    event_attrs[["Event type"]] <- event_type
    # Parse junction positions
    junctions <- as.character(event[[5]])
    junctions <- parseVastToolsJunctions(junctions, event_type)
    return(c(event_attrs, more_attrs, junctions))
}

#' Parse junctions from a VAST-TOOLS event
#'
#' @param coord Character: junction coordinates as returned by VAST-TOOLS;
#' something like \emph{chr6:chr6:41040823,41046768-41046903,41051785}
#' @param event_type Chracter: Type of alternative splicing event
#'
#' @details The currently supported types of alternative splicing events are:
#' \itemize{
#'  \item{\bold{S, C1, C2 or C3}: exon skipping events with increasing degrees 
#'  of complexity}
#'  \item{\bold{MIC}: exon skipping events related to micro-exons}
#'  \item{\bold{Alt3 or Alt5}: alternative 3' or 5' splice site (ALTA or ALTD,
#'  respectively)}
#'  \item{\bold{IR-C or IR-S}: intron retention events (IR-C includes
#'  overlapping other AS events)}
#' }
#' 
#' For more information on VAST-TOOLS, please check out the GitHub repository:
#' \url{https://github.com/vastgroup/vast-tools}.
#'
#' @return List of parsed junctions for a given event
#' @export
#'
#' @examples
#' coord <- "chr6:41040823,41046768-41046903,41051785"
#' parseVastToolsJunctions(coord, event_type = "SE")
parseVastToolsJunctions <- function(coord, event_type) {
    # Fill list of parsed junctions with NAs
    parsed <- list("C1 start" = NA, "C1 end" = NA,
                   "A1 start" = NA, "A1 end" = NA,
                   "A2 start" = NA, "A2 end" = NA,
                   "C2 start" = NA, "C2 end" = NA)
    
    # Get strand for intron retention
    if (event_type == "RI") {
        len <- nchar(coord)
        last <- substr(coord, len, len)
        parsed["Strand"] <- last
    }
    
    # Split event information
    coord <- strsplit(coord, ":|,|-|=")[[1]]
    parsed[["Chromosome"]] <- coord[[1]]
    junctions <- coord[2:length(coord)]
    
    # Split junctions by multiple acceptors/donors
    junctions <- strsplit(junctions, "+", fixed = TRUE)
    junctions <- lapply(junctions, as.numeric)
    
    parseJunctions <- switch(event_type,
                             "SE" = parseVastToolsSE,
                             "RI" = parseVastToolsRI,
                             "A3SS" = parseVastToolsA3SS,
                             "A5SS" = parseVastToolsA5SS)
    parsed <- parseJunctions(junctions, parsed)
    return(parsed)
}

#' Parse junctions of an event from VAST-TOOLS according to event type
#'
#' @param junctions List of integers: exon-exon junctions of an event
#' @param parsed Named list filled with NAs for faster execution (optional)
#'
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (exon skipping)}
#'  \item{\bold{RI} (intron retention)}
#'  \item{\bold{A5SS} (alternative 5' splice site)}
#'  \item{\bold{A3SS} (alternative 3' splice site)}
#' }
#'
#' @seealso \code{\link{parseVastToolsEvent}}
#'
#' @return List of parsed junctions
#' @export
#'
#' @examples
#' junctions <- list(41040823, 41046768, 41046903, 41051785)
#' parseVastToolsSE(junctions)
parseVastToolsSE <- function (junctions, parsed=list()) {
    parsed[["C1 end"]]   <- junctions[[1]]
    parsed[["C2 start"]] <- junctions[[4]]
    # Strand is plus if C1 end is lower than C2 start
    if (junctions[[1]][[1]] < junctions[[4]][[1]]) {
        parsed[["Strand"]]   <- "+"
        parsed[["A1 start"]] <- junctions[[2]]
        parsed[["A1 end"]]   <- junctions[[3]]
    } else {
        parsed[["Strand"]]   <- "-"
        parsed[["A1 start"]] <- junctions[[3]]
        parsed[["A1 end"]]   <- junctions[[2]]
    }
    return(parsed)
}

#' @rdname parseVastToolsSE
#' @examples 
#' 
#' junctions <- list(58861736, 58862017, 58858719, 58859006)
#' parseVastToolsRI(junctions)
parseVastToolsRI <- function (junctions, parsed=list()) {
    if (parsed["Strand"] == "+") {
        parsed[["C1 start"]] <- junctions[[1]]
        parsed[["C1 end"]]   <- junctions[[2]]
        parsed[["C2 start"]] <- junctions[[3]]
        parsed[["C2 end"]]   <- junctions[[4]]
    } else if (parsed["Strand"] == "-") {
        parsed[["C1 start"]] <- junctions[[2]]
        parsed[["C1 end"]]   <- junctions[[1]]
        parsed[["C2 start"]] <- junctions[[4]]
        parsed[["C2 end"]]   <- junctions[[3]]
    }
    return(parsed)
}

#' @rdname parseVastToolsSE
#' @examples 
#' 
#' junctions <- list(49558568, 49557402, c(49557492, 49557470))
#' parseVastToolsA3SS(junctions)
parseVastToolsA3SS <- function (junctions, parsed=list()) {
    parsed[["C1 end"]]   <- junctions[[1]]
    # Strand is plus if C1 end is lower than C2 start
    if (length(junctions[[2]]) > 0 && junctions[[1]][[1]] < junctions[[2]][[1]]) {
        parsed["Strand"] <- "+"
        parsed[["C2 start"]] <- junctions[[2]]
        if (length(junctions) == 3) {
            # if there's C2 end information, save it
            parsed[["C2 end"]] <- junctions[[3]]
        }
    } else {
        parsed["Strand"] <- "-"
        parsed[["C2 start"]] <- junctions[[3]]
        if (length(junctions[[2]]) > 0) {
            # if there's C2 end information, save it
            parsed[["C2 end"]] <- junctions[[2]]
        }
    }
    return(parsed)
}

#' @rdname parseVastToolsSE
#' @examples 
#' 
#' junctions <- list(c(99891605, 99891188), 99891686, 99890743)
#' parseVastToolsA5SS(junctions)
parseVastToolsA5SS <- function (junctions, parsed=list()) {
    parsed[["C2 end"]] <- junctions[[3]]
    # Strand is plus if C1 end is lower than C2 start
    if (length(junctions[[2]]) > 0 && junctions[[2]][[1]] < junctions[[3]][[1]]) {
        parsed["Strand"] <- "+"
        if (length(junctions[[1]]) > 0) {
            # if there's C1 start information, save it
            parsed[["C1 start"]] <- junctions[[1]]
        }
        parsed[["C1 end"]] <- junctions[[2]]
    } else {
        parsed["Strand"] <- "-"
        if (length(junctions[[2]]) > 0) {
            # if there's C1 start information, save it
            parsed[["C1 start"]] <- junctions[[2]]
        }
        parsed[["C1 end"]] <- junctions[[1]]
    }
    return(parsed)
}