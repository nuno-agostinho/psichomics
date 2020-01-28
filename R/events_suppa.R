#' @rdname parseMisoAnnotation
#' @export
#' 
#' @examples 
#' # Load sample files
#' folder <- "extdata/eventsAnnotSample/suppa_output/suppaEvents"
#' suppaOutput <- system.file(folder, package="psichomics")
#' 
#' suppa <- parseSuppaAnnotation(suppaOutput)
parseSuppaAnnotation <- function(
    folder,
    types=c("SE", "AF", "AL", "MX", "A5", "A3", "RI"),
    genome="hg19") {
    
    display("Retrieving SUPPA annotation...")
    typesFile <- file.path(folder, paste0(genome, "_", types, ".ioe"))
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)
    
    display("Parsing SUPPA annotation...")
    eventsID <- lapply(annot, "[[", "event_id")
    events <- lapply(eventsID, parseSuppaEvent)
    events <- rbind.fill(events)
    class(events) <- c("ASevents", class(events))
    return(events)
}

#' Parses splicing events of a specific event type from SUPPA
#'
#' @details More information about SUPPA available at
#' \url{https://bitbucket.org/regulatorygenomicsupf/suppa}
#' 
#' @param event Character vector: Splicing event attributes and junction
#' positions
#'
#' @details The following event types are available to be parsed:
#' \itemize{
#'  \item{\bold{SE} (skipped exon)}
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
#' @keywords internal
#'
#' @examples
#' event <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
#' psichomics:::parseSuppaEvent(event)
parseSuppaEvent <- function(event) {
    # Split event ID by semicolon and colon symbols
    id <- event
    event <- stringr::str_split_fixed(event, pattern = ";|:|-", n = 12)
    event <- data.frame(event, stringsAsFactors = FALSE)
    
    # Create list with event attributes
    event_attrs <- data.frame("Program" = "SUPPA",
                              "Gene" = event[[1]],
                              "Event ID" = id,
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

#' Parse junctions of an event from SUPPA
#'
#' @param junctions List of integers: exon-exon junctions of an event
#' @param strand Character: positive-sense (\code{+}) or negative-sense
#' (\code{-}) strand
#' @param coords Character: coordinate positions to fill
#' @param plus_pos Integer: index of the coordinates for a plus strand event
#' @param minus_pos Integer: index of the coordinates for a minus strand event
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
#' @seealso \code{\link{parseSuppaEvent}()}
#'
#' @return Data frame of parsed junctions
#' @keywords internal
#' 
#' @examples
#' # Parse generic event (in this case, an exon skipping event)
#' junctions <- read.table(text = "169768099 169770024 169770112 169771762")
#' coords <- c("C1.end", "A1.start", "A1.end", "C2.start")
#' plus  <- 1:4
#' minus <- 1:4
#' psichomics:::parseSuppaGeneric(junctions, strand = "+", coords, plus, minus)
parseSuppaGeneric <- function(junctions, strand, coords, plus_pos, minus_pos) {
    # Creates a data frame of parsed junctions filled with NAs
    parsed <- createJunctionsTemplate(nrow(junctions))
    
    plus <- strand == "+"
    parsed[plus, coords] <- junctions[plus, plus_pos]
    parsed[!plus, coords] <- junctions[!plus, minus_pos]
    return(parsed)
}

#' @rdname parseSuppaGeneric
#'
#' @examples
#' 
#' junctions <- read.table(text = "169768099 169770024 169770112 169771762")
#' psichomics:::parseSuppaSE(junctions, "+")
parseSuppaSE <- function (junctions, strand) {
    coords <- c("C1.end", 
                "A1.start", "A1.end",
                "C2.start")
    plus_pos  <- seq(4)
    minus_pos <- 4:1
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' 
#' @examples 
#' 
#' junctions <- read.table(text = "196709749 196709922 196711005 196711181")
#' psichomics:::parseSuppaRI(junctions, "+")
parseSuppaRI <- function (junctions, strand) {
    coords <- c("C1.start", "C1.end",
                "C2.start", "C2.end")
    plus_pos  <- seq(4)
    minus_pos <- 4:1
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' @examples 
#' 
#' junctions <- read.table(
#'     text = "24790610 24792494 24792800 24790610 24795476 24795797")
#' psichomics:::parseSuppaALE(junctions, "+")
parseSuppaALE <- function (junctions, strand) {
    coords <- c("C1.end",
                "A1.start", "A1.end",
                "A2.start", "A2.end")
    plus_pos  <- c(seq(3), 5:6)
    minus_pos <- c(6:4, 2:1)
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' @examples 
#' 
#' junctions <- read.table(
#'     text = "169763871 169764046 169767998 169764550 169765124 169767998")
#' psichomics:::parseSuppaAFE(junctions, "+")
parseSuppaAFE <- function (junctions, strand) {
    coords <- c("A2.start", "A2.end",
                "A1.start", "A1.end",
                "C2.start")
    plus_pos  <- c(4:5, seq(3))
    minus_pos <- c(3:2, 6:4)
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' @examples 
#' 
#' junctions <- read.table(
#'     text = "202060671 202068453 202068489 202073793 202060671 202072798 202072906 202073793")
#' psichomics:::parseSuppaMXE(junctions, "+")
parseSuppaMXE <- function (junctions, strand) {
    coords <- c("C1.end",
                "A1.start", "A1.end",
                "A2.start", "A2.end",
                "C2.start")
    plus_pos  <- c(seq(3), 6:8)
    minus_pos <- c(8:6, 3:1)
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' @examples 
#' 
#' junctions <- read.table(text = "169772450 169773216 169772450 169773253")
#' psichomics:::parseSuppaA3SS(junctions, "+")
parseSuppaA3SS <- function (junctions, strand) {
    coords <- c("C1.end", "A1.start", "A2.start")
    plus_pos  <- c(1, 2, 4)
    minus_pos <- c(4, 3, 1)
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}

#' @rdname parseSuppaGeneric
#' @examples 
#' 
#' junctions <- read.table(text = "50193276 50197008 50192997 50197008")
#' psichomics:::parseSuppaA5SS(junctions, "+")
parseSuppaA5SS <- function (junctions, strand) {
    coords <- c("A2.end", "A1.end", "C2.start")
    plus_pos  <- c(3, 1, 4)
    minus_pos <- c(2, 4, 1)
    parseSuppaGeneric(junctions, strand, coords, plus_pos, minus_pos)
}