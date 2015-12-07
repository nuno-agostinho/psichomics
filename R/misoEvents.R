#' Match MISO's splicing event IDs with the IDs present in the alternative
#' splicing annotation file and get the event as a data frame
#'
#' @details More information about MISO available at 
#' \url{http://miso.readthedocs.org}
#'
#' @param eventID Character vector: alternative event IDs
#' @param annotation Data.frame: alternative event annotation file
#' @param IDcolumn Integer: index of the column with the event ID's in the
#' alternative event annotation file
#'
#' @return Data frame of the matching event
#' @export
#'
#' @examples
#' eventID <- "2217@uc002poi.1@uc002poe.1"
#' annotation <- read.delim("AFE.hg19.gff3", header=FALSE, comment.char="#")
#' IDcolumn <- 9
#' matchMisoEvent(eventID, annotation, IDcolumn)
matchMisoEventID <- function(eventID, annotation, IDcolumn) {
  # Get first line from annotation matching a given splicing event ID
  event <- match(paste0("ID=", eventID,
                        ";Name=", eventID,
                        ";gid=", eventID),
                 annotation[[IDcolumn]])
  # Get all lines related to a given splicing event
  genes <- which(annotation[["V3"]] == "gene")
  # Get index of the next gene
  next_index <- genes[match(event, genes) + 1]
  ret <- lapply(1:length(event), function(i)
    return(annotation[index[i]:(next_index[i] - 1), ])
  )
  return(ret)
}

#' Parse an alternative splicing event from MISO
#'
#' @details More information about MISO available at
#' \url{http://miso.readthedocs.org}
#'
#' @param event Data.frame containing only one event with at least 7 columns as 
#' retrieved from the annotation files from MISO (GFF3 files)
#'
#' @return List with event attributes and junction positions for the exons
#' (depends on the events)
#' @export
#'
#' @examples
#' # alternative splicing event: retained introns (RI)
#' event <- read.table(text = "
#'   chr1 RI gene 17233 17742 . - .
#'   chr1 RI mRNA 17233 17742 . - .
#'   chr1 RI exon 17233 17742 . - .
#'   chr1 RI mRNA 17233 17742 . - .
#'   chr1 RI exon 17233 17364 . - .
#'   chr1 RI exon 17601 17742 . - .")
#' parseMisoEvent(event)
parseMisoEvent <- function(event) {
  event_attrs <- list("Program" = "MISO",
                      "Chromosome" = as.character(event[1, 1]),
                      "Event type" = as.character(event[1, 2]),
                      "Strand"     = as.character(event[1, 7]))
  
  strand     <- event_attrs[["Strand"]]
  event_type <- event_attrs[["Event type"]]
  
  parsed = list("C1 start" = NA, "C1 end" = NA,
                "A1 start" = NA, "A1 end" = NA,
                "A2 start" = NA, "A2 end" = NA,
                "C2 start" = NA, "C2 end" = NA)
  
  parseEventType <- switch(event_type,
                           "SE"   = parseMisoSE,
                           "MXE"  = parseMisoMXE,
                           "RI"   = parseMisoRI,
                           "A5SS" = parseMisoA5SS,
                           "A3SS" = parseMisoA3SS,
                           "AFE"  = parseMisoAFE,
                           "ALE"  = parseMisoALE,
                           "TandemUTR" = parseMisoTandemUTR)
  parsed <- parseEventType(event, strand, parsed)
  return(c(event_attrs, parsed))
}

parseMisoSE <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[c("C1 start", "C1 end")] <- event[3, 4:5]
    parsed[c("A1 start", "A1 end")] <- event[4, 4:5]
    parsed[c("C2 start", "C2 end")] <- event[5, 4:5]
  } else if (strand == "-") {
    parsed[c("C1 start", "C1 end")] <- event[5, 5:4]
    parsed[c("A1 start", "A1 end")] <- event[4, 5:4]
    parsed[c("C2 start", "C2 end")] <- event[3, 5:4]
  }
  return(parsed)
}

parseMisoMXE <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[c("C1 start", "C1 end")] <- event[3, 4:5]
    parsed[c("A1 start", "A1 end")] <- event[4, 4:5]
    parsed[c("A2 start", "A2 end")] <- event[8, 4:5]
    parsed[c("C2 start", "C2 end")] <- event[5, 4:5]
  } else if (strand == "-") {
    parsed[c("C1 start", "C1 end")] <- event[9, 5:4]
    parsed[c("A1 start", "A1 end")] <- event[8, 5:4]
    parsed[c("A2 start", "A2 end")] <- event[4, 5:4]
    parsed[c("C2 start", "C2 end")] <- event[7, 5:4]
  }
  return(parsed)
}

parseMisoRI <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[c("C1 start", "C1 end")] <- event[5, 4:5]
    parsed[c("C2 start", "C2 end")] <- event[6, 4:5]
  } else if (strand == "-") {
    parsed[c("C1 start", "C1 end")] <- event[6, 5:4]
    parsed[c("C2 start", "C2 end")] <- event[5, 5:4]
  }
  return(parsed)
}

parseMisoA5SS <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[["C1 start"]] <- event[3, 4]
    parsed[["C1 end"]]   <- event[c(3, 6), 5]
    parsed[c("C2 start", "C2 end")] <- event[4, 4:5]
  } else if (strand == "-") {
    parsed[["C1 start"]] <- event[4, 5]
    parsed[["C1 end"]]   <- event[c(4, 7), 4]
    parsed[c("C2 start", "C2 end")] <- event[3, 5:4]
  }
  return(parsed)
}

parseMisoA3SS <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[c("C1 start", "C1 end")] <- event[3, 4:5]
    parsed[["C2 start"]] <- event[c(4, 7), 4]
    parsed[["C2 end"]]   <- event[4, 5]
  } else if (strand == "-") {
    parsed[c("C1 start", "C1 end")] <- event[4, 5:4]
    parsed[["C2 start"]] <- event[c(3, 6), 5]
    parsed[["C2 end"]]   <- event[3, 4]
  }
  return(parsed)
}

#' Get a list of each mRNA and respective exons
list_mRNA <- function(event) {
  mRNA_index <- which(event[ , 3] == "mRNA")
  if (length(mRNA_index) == 1) {
    next_index <- nrow(event)  
  } else {
    next_index <- c(mRNA_index[2:length(mRNA_index)] - 1, nrow(event))
  }
  
  # Get each mRNA and respective exons as a new element of a list
  mRNA <- lapply(1:(length(mRNA_index)), function(i)
    return(event[mRNA_index[i]:next_index[i], 1:8])
  )
  return(mRNA)
}

#' Clear mRNAs with the same exons from a given list of mRNAs
remove_duplicated_mRNA <- function (mRNA) {
  # Get first occurence of each mRNA and remove duplicated index
  uniq <- unique(match(mRNA, mRNA))
  
  # Return a non-redundant list of mRNAs
  return(mRNA[uniq])
}

#' Remove wrong mRNAs of a given event
remove_wrong_mRNA <- function(event) {
  # Clear mRNA and exons with different chromosome identifier
  chr      <- event[1, 1]
  same_chr <- event[ , 1] == chr
  event <- event[same_chr, ]
  
  # Clear mRNA and exons outside event boundaries
  start  <- event[1, 4]
  end    <- event[1, 5]
  inside <- event[ , 4] >= start & event[ , 5] <= end
  event  <- event[inside, ]
  
  return(event)
}

parseMisoAFE <- function(event, strand, parsed) {
  # Remove mRNAs from different chromosomes and mark event
  len <- nrow(event)
  event <- remove_wrong_mRNA(event)
  if (len < nrow(event)) {
    parsed[["MISO condition"]] <- "wrong mRNAs"
  }
  
  # Get a list of each mRNA and respective exons
  mRNA <- list_mRNA(event)
  
  # Remove mRNAs with the same exons and mark event
  len <- length(mRNA)
  mRNA <- remove_duplicated_mRNA(mRNA)
  if (len < length(mRNA)) {
    parsed[["MISO condition"]] <- "duplicated mRNAs"
  }
  
  if (length(mRNA) != 2) {
    # Store (but don't parse) events with more than two mRNAs or only one mRNA
    parsed[["MISO condition"]] <- "unrecognized event"
    parsed[["MISO event"]] <- event
  } else if (strand == "+") {
    mRNA1 <- mRNA[[1]]
    exon1 <- mRNA1[nrow(mRNA1), 4:5]
    mRNA2 <- mRNA[[2]]
    exon2 <- mRNA2[nrow(mRNA2), 4:5]
    # Check if the most downstream exons are equal in both mRNAs
    if (all(exon1 == exon2)) {
      # Save positions of shared (constitutive) exon and alternative
      # first exons
      parsed[c("C1 start", "C1 end")] <- mRNA1[nrow(mRNA2) - 1, 4:5]
      parsed[c("A1 start", "A1 end")] <- mRNA1[nrow(mRNA2) - 1, 4:5]
      parsed[c("C2 start", "C2 end")] <- exon1
    } else {
      # Save the positions of the downstream exons
      parsed[c("C1 start", "C1 end")] <- exon1
      parsed[c("A1 start", "A1 end")] <- exon2
    }
  } else if (strand == "-") {
    mRNA1 <- mRNA[[1]]
    exon1 <- mRNA1[2, 5:4]
    mRNA2 <- mRNA[[2]]
    exon2 <- mRNA2[2, 5:4]
    # Check if the most downstream exons are equal in both mRNAs
    if (all(exon1 == exon2)) {
      # Save positions of shared (constitutive) exon and alternative
      # first exons
      parsed[c("C1 start", "C1 end")] <- mRNA1[3, 5:4]
      parsed[c("A1 start", "A1 end")] <- mRNA2[3, 5:4]
      parsed[c("C2 start", "C2 end")] <- exon1
    } else {
      # Save the positions of the downstream exons
      parsed[c("C1 start", "C1 end")] <- exon1
      parsed[c("A1 start", "A1 end")] <- exon2
    }
  }
  return(parsed)
}

parseMisoALE <- function(event, strand, parsed) {
  len <- nrow(event)
  if (len == 1) {
    parsed[["MISO condition"]] <- "zero mRNAs"
    parsed[["MISO event"]] <- event
  } else if (len == 5) {
    # Most ALE events have length of 5, so let's avoid wasting time
    if (strand == "+") {
      parsed[c("A1 start", "A1 end")] <- event[3, 4:5]
      parsed[c("C2 start", "C2 end")] <- event[5, 4:5]
    } else if (strand == "-") {
      parsed[c("A1 start", "A1 end")] <- event[5, 5:4]
      parsed[c("C2 start", "C2 end")] <- event[3, 5:4]
    }
  } else {
    event <- remove_wrong_mRNA(event)
    if (len < nrow(event)) {
      parsed[["MISO condition"]] <- "wrong mRNAs"
    }
    
    # Get a list of each mRNA and respective exons
    mRNA <- list_mRNA(event)
    
    # Remove mRNAs with the same exons and mark event
    len <- length(mRNA)
    mRNA <- remove_duplicated_mRNA(mRNA)
    if (len < length(mRNA)) {
      parsed[["MISO condition"]] <- "duplicated mRNAs"
    }
    
    if (length(mRNA) != 2) {
      # Store (but don't parse) events with more than two mRNAs or only one mRNA
      parsed[["MISO condition"]] <- "unrecognized event"
      parsed[["MISO event"]] <- event
    } else if (strand == "+") {
      mRNA1 <- mRNA[[1]]
      exon1 <- mRNA1[2, 4:5]
      mRNA2 <- mRNA[[2]]
      exon2 <- mRNA2[2, 4:5]
      # Check if the most downstream exons are equal in both mRNAs
      if (all(exon1 == exon2)) {
        parsed[c("C1 start", "C1 end")] <- exon1
        parsed[c("A1 start", "A1 end")] <- mRNA1[3, 4:5]
        parsed[c("C2 start", "C2 end")] <- mRNA2[3, 4:5]
      } else {
        parsed[c("A1 start", "A1 end")] <- exon1
        parsed[c("C2 start", "C2 end")] <- exon2
      }
    } else if (strand == "-") {
      mRNA1 <- mRNA[[1]]
      exon1 <- mRNA1[nrow(mRNA1), 5:4]
      mRNA2 <- mRNA[[2]]
      exon2 <- mRNA2[nrow(mRNA2), 5:4]
      # Check if the most downstream exons are equal in both mRNAs
      if (all(exon1 == exon2)) {
        parsed[c("C1 start", "C1 end")] <- exon1
        parsed[c("A1 start", "A1 end")] <- mRNA1[nrow(mRNA2) - 1, 5:4]
        parsed[c("C2 start", "C2 end")] <- mRNA2[nrow(mRNA2) - 1, 5:4]
      } else {
        parsed[c("A1 start", "A1 end")] <- exon1
        parsed[c("C2 start", "C2 end")] <- exon2
      }
    }
  }
  return(parsed)
}

parseMisoTandemUTR <- function(event, strand, parsed) {
  if (strand == "+") {
    parsed[["C2 start"]] <- event[3, 4]
    parsed[["C2 end"]]   <- event[c(3, 5), 5]
  } else if (strand == "-") {
    parsed[["C2 start"]] <- event[3, 5]
    parsed[["C2 end"]]   <- event[c(3, 5), 4]
  }
  return(parsed)
}