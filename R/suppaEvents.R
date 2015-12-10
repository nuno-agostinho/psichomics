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
#' parseSuppaEvent(event)
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
                      "Event type" = event[2],
                      "Chromosome" = event[3],
                      "Strand" = event[len])
  
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
  switch(event_type,
         "A3" = {
           parsed[["C1 end"]]   <- junctions[1]
           # it's important to conserve the order of which exon appears first
           if (strand == "+") {
            parsed[["C2 start"]] <- junctions[c(2, 4)]
           } else if (strand == "-") {
             parsed[["C2 start"]] <- junctions[c(4, 2)]
           }
         }, # end A3
         "A5" = {
           # it's important to conserve the order of which exon appears first
           if (strand == "+") {
             parsed[["C1 end"]]   <- junctions[c(1, 3)]
           } else if (strand == "-") {
             parsed[["C1 end"]]   <- junctions[c(3, 1)]
           }
           parsed[["C2 start"]] <- junctions[2]
         }, # end A5
         "SE" = parsed[c("C1 end",
                         "A1 start", "A1 end",
                         "C2 start")] <- junctions,
         "MX" = {
           # it's important to conserve the order of which exon appears first
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
         }, # end MX
         "RI" = parsed[c("C1 start", "C1 end",
                         "C2 start", "C2 end")] <- junctions,
         "AF" = parsed[c("C1 start", "C1 end",
                         "C2 start",
                         "A1 start", "A1 end")] <- junctions[1:5],
         "AL" = parsed[c("C2 start", "C2 end",
                         "C1 end",
                         "A1 start", "A1 end")] <- junctions[2:6]
  ) # end switch
  return(parsed)
}