parseMultipleMatsEvents <- function(events, event_type) {
    return(lapply(1:nrow(events),
                  function(k) {
                      print(events[k, ])
                      parseMatsEvent(events[k, ], event_type
                      )}
                  ))
}

#' Parse MATS alternative splicing events
#'
#' @param event Data frame: MATS splicing event
#' @param event_type Character: Type of event to be parsed (see details)
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
#' mats_A3SS <- c(ID = "1234", GeneID = "ENSG00000076108", geneSymbol = "BAZ2A", 
#'   chr = "chr12", strand = "-", longExonStart_0base = "57000030", 
#'   longExonEnd = "57000179", shortES = "57000030", shortEE = "57000096",
#'   flankingES = "57000416", flankingEE = "57000517", ID.1 = "1234", 
#'   IJC_SAMPLE_1 = "0", SJC_SAMPLE_1 = "18", IJC_SAMPLE_2 = "2",
#'   SJC_SAMPLE_2 = "0", IncFormLen = "112", SkipFormLen = "56",
#'   PValue = "0.00337999081157", FDR = "0.0540798529851", IncLevel1 = "0",
#'   IncLevel2 = "1", IncLevelDifference = "-1")
#' 
#' parseMatsEvent(mats_A3SS, "A3SS")
parseMatsEvent <- function(event, event_type) {
    len <- length(event)
    # Create list with event attributes
    event_attrs <- list("Program" = "MATS",
                        "Gene" = as.character(event[[2]]),
                        "Gene symbol" = as.character(event[[3]]),
                        "Chromosome" = as.character(event[[4]]),
                        "Strand" = as.character(event[[5]]),
                        "Event type" = event_type)
    
    # Add these attributes only if the event has a given length
    if (len > 13) {
        more_attrs <- list("P value" = as.numeric(event[[len - 4]]),
                           "FDR" = as.numeric(event[[len - 3]]),
                           "Inclusion level A" = as.numeric(event[[len - 2]]),
                           "Inclusion level B" = as.numeric(event[[len - 1]]))
    } else {
        more_attrs <- list()
    }
    
    # Parse junction positions according to event type
    parsed <- parseMatsJunctions(event, event_type)
    return(c(event_attrs, more_attrs, parsed))
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
#' mats_A3SS <- c(ID = "1234", GeneID = "ENSG00000076108", geneSymbol = "BAZ2A", 
#'   chr = "chr12", strand = "-", longExonStart_0base = "57000030", 
#'   longExonEnd = "57000179", shortES = "57000030", shortEE = "57000096",
#'   flankingES = "57000416", flankingEE = "57000517", ID.1 = "1234", 
#'   IJC_SAMPLE_1 = "0", SJC_SAMPLE_1 = "18", IJC_SAMPLE_2 = "2",
#'   SJC_SAMPLE_2 = "0", IncFormLen = "112", SkipFormLen = "56",
#'   PValue = "0.00337999081157", FDR = "0.0540798529851", IncLevel1 = "0",
#'   IncLevel2 = "1", IncLevelDifference = "-1")
#' 
#' parseMatsJunctions(mats_A3SS, "A3SS")
parseMatsJunctions <- function(event, event_type) {
    # Fill list of parsed junctions with NAs
    parsed <- list("C1 start" = NA, "C1 end" = NA,
                   "A1 start" = NA, "A1 end" = NA,
                   "A2 start" = NA, "A2 end" = NA,
                   "C2 start" = NA, "C2 end" = NA)
    
    strand <- as.character(event[[5]])
    
    # Parse junction positions according to event type
    switch(event_type,
           "SE" = {
               if (strand == "+") {
                   parsed[c("A1 start", "A1 end",
                            "C1 start", "C1 end",
                            "C2 start", "C2 end")] <- as.numeric(event[6:11])
               } else if (strand == "-") {
                   parsed[c("A1 end", "A1 start",
                            "C2 end", "C2 start",
                            "C1 end", "C1 start")] <- as.numeric(event[6:11])
               }
           },
           "MXE" = {
               if (strand == "+") {
                   parsed[c("A1 start", "A1 end",
                            "A2 start", "A2 end",
                            "C1 start", "C1 end",
                            "C2 start", "C2 end")] <- as.numeric(event[6:13])
               } else if (strand == "-") {
                   parsed[c("A1 end", "A1 start",
                            "A2 end", "A2 start",
                            "C2 end", "C2 start",
                            "C1 end", "C1 start")] <- as.numeric(event[6:13])
               }
           },
           "RI" = {
               if (strand == "+") {
                   parsed[c("C1 start", "C1 end",
                            "C2 start", "C2 end")] <- as.numeric(event[8:11])
               } else if (strand == "-") {
                   parsed[c("C1 start", "C1 end",
                            "C2 start", "C2 end")] <- as.numeric(event[11:8])
               }
           },
           "A3SS" = {
               junctions <- as.numeric(event[6:11])
               if (strand == "+") {
                   parsed[c("C1 start", "C1 end")] <- junctions[5:6]
                   parsed[["C2 start"]] <- junctions[c(1, 3)]
                   parsed[["C2 end"]]   <- junctions[2]
               } else if (strand == "-") {
                   parsed[c("C1 start", "C1 end")] <- junctions[6:5]
                   parsed[["C2 start"]]   <- junctions[c(2, 4)]
                   parsed[["C2 end"]]  <- junctions[1]
               }
           },
           "A5SS" = {
               junctions <- as.numeric(event[6:11])
               if (strand == "+") {
                   parsed[["C1 start"]] <- junctions[1]
                   parsed[["C1 end"]]   <- junctions[c(2,4)]
                   parsed[c("C2 start", "C2 end")] <- junctions[5:6]
               } else if (strand == "-") {
                   parsed[["C1 start"]]  <- junctions[2]
                   parsed[["C1 end"]]    <- junctions[c(1, 3)]
                   parsed[c("C2 start", "C2 end")] <- junctions[6:5]
               }
           },
           "AFE" = {
               junctions <- as.numeric(event[6:11])
               if (strand == "+") {
                   parsed[c("C1 start", "C1 end")] <- junctions[1:2]
                   parsed[c("A1 start", "A1 end")] <- junctions[3:4]
                   parsed[c("C2 start", "C2 end")] <- junctions[5:6]
               } else if (strand == "-") {
                   parsed[c("C1 start", "C1 end")] <- junctions[2:1]
                   parsed[c("A1 start", "A1 end")] <- junctions[4:3]
                   parsed[c("C2 start", "C2 end")] <- junctions[6:5]
               }
           },
           "ALE" = {
               junctions <- as.numeric(event[6:11])
               if (strand == "+") {
                   parsed[c("C1 start", "C1 end")] <- junctions[5:6]
                   parsed[c("A1 start", "A1 end")] <- junctions[1:2]
                   parsed[c("C2 start", "C2 end")] <- junctions[3:4]
               } else if (strand == "-") {
                   parsed[c("C1 start", "C1 end")] <- junctions[6:5]
                   parsed[c("A1 start", "A1 end")] <- junctions[2:1]
                   parsed[c("C2 start", "C2 end")] <- junctions[4:3]
               }
           }
    ) # end of switch
    return(parsed)
}