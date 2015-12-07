parseMatsEvent <- function(event, event_type) {
  # Create list with event attributes
  event_attrs <- list("Program" = "MATS",
                      "Gene" = as.character(event[[2]]),
                      "Gene symbol" = as.character(event[[3]]),
                      "Chromosome" = as.character(event[[4]]),
                      "Strand" = as.character(event[[5]]),
                      "P value" = as.numeric(event[[len - 4]]),
                      "FDR" = as.numeric(event[[len - 3]]),
                      "Inclusion level A" = as.numeric(event[[len - 2]]),
                      "Inclusion level B" = as.numeric(event[[len - 1]]),
                      "Event type" = event_type)
  
  # Parse junction positions according to event type
  parsed <- parseMatsJunctions(event, event_type)
  return(c(event_attrs, parsed))
}

parseMatsJunctions <- function(event, event_type) {
  # Fill list of parsed junctions with NAs
  parsed <- list("C1 start" = NA, "C1 end" = NA,
                 "A1 start" = NA, "A1 end" = NA,
                 "A2 start" = NA, "A2 end" = NA,
                 "C2 start" = NA, "C2 end" = NA)
  
  # Parse junction positions according to event type
  switch(event_type,
         "SE"   = parsed[c("A1 start", "A1 end",
                           "C1 start", "C1 end",
                           "C2 start", "C2 end")] <- as.numeric(event[6:11]),
         "MXE"  = parsed[c("A1 start", "A1 end",
                           "A2 start", "A2 end",
                           "C1 start", "C1 end",
                           "C2 start", "C2 end")] <- as.numeric(event[6:13]),
         "RI"   = parsed[c("C1 start", "C1 end",
                           "C2 start", "C2 end")] <- as.numeric(event[8:11]),
         "A3SS" = {
           parsed[c("C1 start", "C1 end")] <- as.numeric(event[10:11])
           if (strand == "+") {
             parsed[["C2 start"]] <- as.numeric(event[c(6, 8)])
             parsed[["C2 end"]]   <- as.numeric(event[7])
           } else if (strand == "-") {
             parsed[["C2 start"]]   <- as.numeric(event[c(7, 9)])
             parsed[["C2 end"]]  <- as.numeric(event[6])
           }
         },
         "A5SS" = {
           parsed[c("C2 start", "C2 end")] <- as.numeric(event[10:11])
           if (strand == "+") {
             parsed[["C1 start"]] <- as.numeric(event[6])
             parsed[["C1 end"]]   <- as.numeric(event[c(7,9)])
           } else if (strand == "-") {
             parsed[["C1 start"]]  <- as.numeric(event[7])
             parsed[["C1 end"]]    <- as.numeric(event[c(6, 8)])
           }
         }
  ) # end of switch
  return(parsed)
}