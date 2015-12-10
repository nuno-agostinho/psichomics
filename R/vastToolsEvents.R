#' Parses an alternative splicing event from VAST-TOOLS
#'
#' @details Junctions are parsed from 
#' 
#' @param event Data.frame row, list or character: Vast-Tools event containing
#' gene symbol, event ID, length, junctions coordinates, event type and inclusion 
#' levels plus  for both samples
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
  event_type <- as.character(event[[6]])
  event_attrs <- list("Program" = "VAST-TOOLS",
                      "Gene symbol" = as.character(event[[1]]),
                      "Event ID" = as.character(event[[2]]),
                      "Event type" = event_type,
                      "Inclusion level A" = as.numeric(event[[7]]),
                      "Inclusion level B" = as.numeric(event[[9]]))
  
  # Parse junction positions
  junctions <- as.character(event[[5]])
  junctions <- parseVastToolsJunctions(junctions, event_type)
  return(c(event_attrs, junctions))
}

#' Parse junctions from a VAST-TOOLS event
#'
#' @param coord Character: junction coordinates as returned by VAST-TOOLS
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
#' parseVastToolsJunctions(coord, event_type = "C2")
parseVastToolsJunctions <- function(coord, event_type) {
  # Fill list of parsed junctions with NAs
  parsed <- list("C1 start" = NA, "C1 end" = NA,
                 "A1 start" = NA, "A1 end" = NA,
                 "A2 start" = NA, "A2 end" = NA,
                 "C2 start" = NA, "C2 end" = NA)
  
  # Get strand for intron retention
  len <- nchar(coord)
  last <- substr(coord, len, len)
  
  # Split event information
  coord <- strsplit(coord, ":|,|-|=")[[1]]
  parsed[["Chromosome"]] <- coord[[1]]
  junctions <- coord[2:length(coord)]
  
  # Split junctions by multiple acceptors/donors
  junctions <- strsplit(junctions, "+", fixed = TRUE)
  junctions <- lapply(junctions, as.numeric)
  
  if (is.element(event_type, c("S", "C1", "C2", "C3", "MIC"))) {
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
  } else if (is.element(event_type, c("IR-C", "IR-S"))) {
    if (last == "+") {
      parsed["Strand"] <- "+"
      parsed[["C1 start"]] <- junctions[[1]]
      parsed[["C1 end"]]   <- junctions[[2]]
      parsed[["C2 start"]] <- junctions[[3]]
      parsed[["C2 end"]]   <- junctions[[4]]
    } else if (last == "-") {
      parsed["Strand"] <- "-"
      parsed[["C1 start"]] <- junctions[[2]]
      parsed[["C1 end"]]   <- junctions[[1]]
      parsed[["C2 start"]] <- junctions[[4]]
      parsed[["C2 end"]]   <- junctions[[3]]
    }
  } else if (event_type == "Alt3") {
    parsed[["C1 end"]]   <- junctions[[1]]
    # Strand is plus if C1 end is lower than C2 start
    if (junctions[[1]][[1]] < junctions[[2]][[1]]) {
      parsed["Strand"] <- "+"
      parsed[["C2 start"]] <- junctions[[2]]
      parsed[["C2 end"]]   <- junctions[[3]]
    } else {
      parsed["Strand"] <- "-"
      parsed[["C2 start"]] <- junctions[[3]]
      parsed[["C2 end"]]   <- junctions[[2]]
    }
  } else if (event_type == "Alt5") {
    parsed[["C2 end"]]   <- junctions[[3]]
    # Strand is plus if C1 end is lower than C2 start
    if (junctions[[1]][[1]] < junctions[[3]][[1]]) {
      parsed["Strand"] <- "+"
      parsed[["C1 start"]] <- junctions[[1]]
      parsed[["C1 end"]]   <- junctions[[2]]
    } else {
      parsed["Strand"] <- "-"
      parsed[["C1 start"]] <- junctions[[2]]
      parsed[["C1 end"]]   <- junctions[[1]]
    }
  }
  return(parsed)
}