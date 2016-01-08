#' Get rows of a data frame between two row indexes
#'
#' @details For a given iteration i, returns data from first_row[i] to
#' last_row[i]
#'
#' @param i Integer: current iteration
#' @param data Data.frame: contains the data of interest
#' @param first_row Vector of integers: First row index of interest; value must
#' be less than the respective last row index and less than the number of rows
#' in the data frame
#' @param last_row Vector of integers: Last row index of interest; value must be 
#' higher than the respective first row index and less than the number of rows
#' in the data frame
#' 
#' @return Data frame subset givne two row indexes (returns NA if the first row 
#' index is NA)
#' @export
getDataRows <- function(i, data, first_row, last_row) {
    first <- first_row[i]
    last  <- last_row[i]
    if (is.na(first)) {
        # if there is no first position, there is no item match
        return(NA)
    } else if (is.na(last)) {
        # if there's no last position, use the last row of the data frame
        last <- nrow(data)
    }
    seq <- seq(first, last)
    return(data[seq, 1:8])
}

#' Match MISO's splicing event IDs with the IDs present in the alternative
#' splicing annotation file and get events in a data frame
#'
#' @details For faster execution times, provide a vector of event IDs.
#' 
#' For more information about MISO, see \url{http://miso.readthedocs.org}.
#'
#' @param eventID Character: alternative event IDs
#' @param annotation Data.frame: alternative event annotation file
#' @param IDcolumn Integer: index of the column with the event ID's in the
#' alternative event annotation file
#'
#' @note If possible, it's recommend to use smaller subsets of the alternative
#' events' annotation instead of all data for faster runs. For example, when
#' trying to match only skipped exons event IDs, only use the annotation of
#' skipped exons instead of using a mega annotation with all event types.
#'
#' @importFrom fastmatch fmatch
#'
#' @return Data frame of the matching events (or NA when nothing is matched)
#' @export
#'
#' @examples
#' eventID <- c("2217@uc002poi.1@uc002poe.1", "57705@uc009xob.1@uc001jgy.2")
#' # the annotation is one of the GFF3 files needed to run MISO
#' annotation <- read.delim("AFE.hg19.gff3", header=FALSE, comment.char="#")
#' IDcolumn <- 9
#' parseMisoEventID(eventID, annotation, IDcolumn)
parseMisoEventID <- function(eventID, annotation, IDcolumn) {
    # Get first row from annotation matching a given splicing event ID
    index <- fmatch(paste0("ID=", eventID,
                           ";Name=", eventID,
                           ";gid=", eventID),
                    annotation[[IDcolumn]])
    # Get every index of splicing events present in the annotation
    events <- which(annotation[["V3"]] == "gene")
    # Get the index of the next gene
    next_index <- events[fmatch(index, events) + 1]
    
    # Get the rows relative to each event
    rows <- lapply(1:(length(index)), getDataRows, annotation, index,
                   next_index - 1)
    return(rows)
}

#' Parses multiple alternative splicing events from MISO
#' 
#' Expects the annotation for many events at once and parses each event.
#' 
#' @param events Data.frame containing one or more events with at least 7
#' columns as  retrieved from the alternative splicing annotation files from
#' MISO (GFF3 files)
#' @param progress Boolean: show progress?
#' 
#' @seealso \code{\link{parseMisoEvent}}
#' 
#' @export
#' @examples 
#' # example of alternative splicing events: skipped exon (SE)
#' events <- read.table(text = "
#'   chr1 SE gene 16854	18061	. - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17055 . - .
#'   chr1 SE exon 17233 17742 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17955 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr2 SE gene 6854 8061	. + .
#'   chr2 SE mRNA 6854 8061 . + .
#'   chr2 SE exon 6854 7055 . + .
#'   chr2 SE exon 7233 7742 . + .
#'   chr2 SE exon 7915 8061 . + .
#'   chr2 SE mRNA 6854 8061 . + .
#'   chr2 SE exon 6854 7955 . + .
#'   chr2 SE exon 7915 8061 . + .")
#' parseMultipleMisoEvents(events)
# parseMultipleMisoEvents <- function(events, progress=FALSE) {
#     firstIndex <- which(events[ , 3] == "gene")
#     if (length(firstIndex) > 1) {
#         # if there are many events, get each event's last index
#         lastIndex <- c(firstIndex[2:length(firstIndex)] - 1, nrow(events))
#     } else {
#         # if there is only one event, get last index
#         lastIndex <- nrow(events)
#     }
#     
#     if (progress)
#         pb <- txtProgressBar(min=1, max=length(firstIndex), style=3)
#     
#     res <- lapply(1:length(firstIndex), function(i, data) {
#         if (progress)
#             setTxtProgressBar(pb, i)
#         event <- data[firstIndex[i]:lastIndex[i], ]
#         parseMisoEvent(event)
#     }, events)
#     return(res)
# }

#' Filters the events with valid elements according to the given validator
#'
#' @param events Data.frame: events annotation from MISO  
#' @param validator Character: valid elements for each event
#' @param isRunLengthValid Boolean: consider runs of equal elements as valid
#' when comparing with the validator? Default if FALSE (see details)
#'
#' @details \code{isRunLengthValid} allows to consider runs of equal elements
#' (i.e. sequences where the same element consecutively occurs) as valid when
#' comparing with a given validator. For example, if the validator is
#' \code{c("gene", "mRNA", "exon")} and \code{isRunLengthValid = FALSE}, this
#' function will only considerate events as valid if they have the exact same 
#' elements. If \code{isRunLengthValid = TRUE}, a valid events could include the
#' elements \code{c("gene", "mRNA", "exon", "exon", "exon")}.
#'
#' @return Data.frame with valid events
#' @export
#'
#' @examples
#' event <- read.table(text = "
#'  chr1 SE gene 17233 18061  .  -  .
#'  chr1 SE dkfd 00000 30000  .  -  .
#'  chr1 SE mRNA 17233 18061  .  -  .
#'  chr1 SE exon 17233 17368  .  -  .
#'  chr1 SE exon 17526 17742  .  -  .
#'  chr1 SE exon 17915 18061  .  -  .
#'  chr1 SE mRNA 17233 18061  .  -  .
#'  chr1 SE exon 17233 17368  .  -  .
#'  chr1 SE exon 17915 18061  .  -  .
#'  chr1 SE gene 17233 18061  .  -  .
#'  chr1 SE mRNA 17233 18061  .  -  .
#'  chr1 SE exon 17233 17368  .  -  .
#'  chr1 SE exon 17606 17742  .  -  .
#'  chr1 SE exon 17915 18061  .  -  .
#'  chr1 SE mRNA 17233 18061  .  -  .
#'  chr1 SE exon 17233 17368  .  -  .
#'  chr1 SE exon 17915 18061  .  -  .
#' ")
#' validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
#' getValidEvents(event, validator)
getValidEvents <- function(events, validator, isRunLengthValid = FALSE) {
    elem <- events[[3]]
    
    # If run length is valid, consider runs of equal elements as one element
    if (isRunLengthValid)
        # TODO(NunoA): rle should only consider exons (if there are two
        # consecutive mRNA), there could be problems
        elem <- rle(as.character(elem))$values
    
    # Get starting position of each event
    index <- which(elem == "gene")
    
    # Get starting position of the next event
    next_index <- length(elem) + 1
    if (length(index) > 1)
        next_index <- c(index[2:length(index)], next_index)
    
    # Get the elements composing each event as elements of a list
    diff <- next_index - index
    groups <- unlist(lapply(1:length(diff), function(i, diff)
        rep.int(i, times = diff[i]), diff))
    splitElems <- split(elem, groups)
    
    # Check if those elements are valid
    valid <- vapply(splitElems, function(e)
        identical(as.character(e), validator), logical(1))
    
    if (sum(!valid) == 0) {
        # If all events are valid, return all events
        return(events)
    } else if (sum(valid) > 0) {
        # If there are invalid events, return only valid events
        if (isRunLengthValid) {
            # If run length is valid, consider runs of equal elements as one
            
            # Get starting position of each event
            index <- which(events[[3]] == "gene")
            
            # Get starting position of the next event
            next_index <- nrow(events) + 1
            if (length(index) > 1)
                next_index <- c(index[2:length(index)], next_index)
        }
        
        invalid_index <- unlist(
            lapply(1:sum(!valid),
                   function(i) index[!valid][i]:(next_index[!valid][i]-1)
            )
        )
        return(events[-invalid_index, ])
    }
}

#' Parse an alternative splicing event from MISO
#'
#' @details More information about MISO available at
#' \url{http://miso.readthedocs.org}
#'
#' @param event Data.frame containing only one event with at least 7 columns as 
#' retrieved from the alternative splicing annotation files from MISO (GFF3
#' files)
#'
#' @return List with event attributes and junction positions for the exons
#' (depends on the events)
#' @export
#'
#' @examples
#' # example of alternative splicing event: skipped exon (SE)
#' event <- read.table(text = "
#'   chr1 SE gene 16854	18061	. - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17055 . - .
#'   chr1 SE exon 17233 17742 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17955 . - .
#'   chr1 SE exon 17915 18061 . - .")
#' parseMisoEvent(event)
parseMisoEvent <- function(event) {
    event_type <- as.character(event[1, 2])
    parseEvent <- switch(event_type,
                         "SE"   = parseMisoSE,
                         "MXE"  = parseMisoMXE,
                         "RI"   = parseMisoRI,
                         "A5SS" = parseMisoA5SS,
                         "A3SS" = parseMisoA3SS,
                         "AFE"  = parseMisoAFE,
                         "ALE"  = parseMisoALE,
                         "TandemUTR" = parseMisoTandemUTR)
    parsed <- parseEvent(event)
    return(parsed)
}

#' Parse junctions of an event from MISO according to event type
#'
#' @inheritParams parseMisoEvent
#' @param strand Character: strand of the event
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
#'  \item{\bold{Tandem UTR}}
#' }
#' 
#' @seealso \code{\link{parseMisoEvent}}
#'
#' @return List of parsed junctions
#' @export
#'
#' @examples
#' # skipping exon event (SE)
#' event <- read.table(text = "
#'   chr1 SE gene 16854 18061 . - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17055 . - .
#'   chr1 SE exon 17233 17742 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17955 . - .
#'   chr1 SE exon 17915 18061 . - .")
#' parseMisoSE(event)
parseMisoSE <- function(event) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "SE",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("C1.start", 
                             "C1.end")] <- event[index[plus] + 2, 4:5]
            parsed[plus, ][c("A1.start", 
                             "A1.end")] <- event[index[plus] + 3, 4:5]
            parsed[plus, ][c("C2.start", 
                             "C2.end")] <- event[index[plus] + 4, 4:5]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("C1.start", 
                              "C1.end")] <- event[index[minus] + 4, 5:4]
            parsed[minus, ][c("A1.start", 
                              "A1.end")] <- event[index[minus] + 3, 5:4]
            parsed[minus, ][c("C2.start", 
                              "C2.end")] <- event[index[minus] + 2, 5:4]
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#' 
#' # mutually exclusive exon (MXE) event
#' event <- read.table(text = "
#'  chr1 MXE gene 764383 788090 . + .
#'  chr1 MXE mRNA 764383 788090 . + .
#'  chr1 MXE exon 764383 764484 . + .
#'  chr1 MXE exon 776580 776753 . + .
#'  chr1 MXE exon 787307 788090 . + .
#'  chr1 MXE mRNA 764383 788090 . + .
#'  chr1 MXE exon 764383 764484 . + .
#'  chr1 MXE exon 783034 783186 . + .
#'  chr1 MXE exon 787307 788090 . + .")
#' parseMisoMXE(event)
parseMisoMXE <- function(event) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 3))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "MXE",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("C1.start", 
                             "C1.end")] <- event[index[plus] + 2, 4:5]
            parsed[plus, ][c("A1.start", 
                             "A1.end")] <- event[index[plus] + 3, 4:5]
            parsed[plus, ][c("A2.start", 
                             "A2.end")] <- event[index[plus] + 7, 4:5]
            parsed[plus, ][c("C2.start", 
                             "C2.end")] <- event[index[plus] + 4, 4:5]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("C1.start",
                              "C1.end")] <- event[index[minus] + 8, 5:4]
            parsed[minus, ][c("A1.start", 
                              "A1.end")] <- event[index[minus] + 3, 5:4]
            parsed[minus, ][c("A2.start", 
                              "A2.end")] <- event[index[minus] + 7, 5:4]
            parsed[minus, ][c("C2.start", 
                              "C2.end")] <- event[index[minus] + 6, 5:4]
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#'
#' # intron retention (RI) event
#' event <- read.table(text = "
#'  chr1 RI gene 17233 17742 . - .
#'  chr1 RI mRNA 17233 17742 . - .
#'  chr1 RI exon 17233 17742 . - .
#'  chr1 RI mRNA 17233 17742 . - .
#'  chr1 RI exon 17233 17364 . - .
#'  chr1 RI exon 17601 17742 . - .")
#' parseMisoRI(event)
parseMisoRI <- function(event, strand) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", "exon", "mRNA", rep("exon", 2))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "RI",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("C1.start", 
                             "C1.end")] <- event[index[plus] + 4, 4:5]
            parsed[plus, ][c("C2.start",
                             "C2.end")] <- event[index[plus] + 5, 4:5]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("C1.start",
                              "C1.end")] <- event[index[minus] + 5, 5:4]
            parsed[minus, ][c("C2.start",
                              "C2.end")] <- event[index[minus] + 4, 5:4]
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#'
#' # alternative 5' splice site (A5SS) event
#' event <- read.table(text = "
#'  chr1 A5SS gene 17233 17742 . - .
#'  chr1 A5SS mRNA 17233 17742 . - .
#'  chr1 A5SS exon 17233 17368 . - .
#'  chr1 A5SS exon 17526 17742 . - .
#'  chr1 A5SS mRNA 17233 17742 . - .
#'  chr1 A5SS exon 17233 17368 . - .
#'  chr1 A5SS exon 17606 17742 . - .")
#' parseMisoA5SS(event)
parseMisoA5SS <- function(event) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", rep("exon", 2), "mRNA", rep("exon", 2))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get the first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "A5SS",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ]["C1.start"] <- event[index[plus] + 2, 4]
            parsed[plus, ][["C1.end"]] <-
                lapply(index[plus], function(i) event[i + c(2, 5), 5])
            parsed[plus, ][c("C2.start", "C2.end")] <- 
                event[index[plus] + 3, 4:5]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ]["C1.start"] <- event[index[minus] + 3, 5]
            parsed[minus, ][["C1.end"]] <-
                lapply(index[minus], function(i) event[i + c(3, 6), 4])
            parsed[minus, ][c("C2.start", "C2.end")] <- 
                event[index[minus]+2, 5:4]
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#'
#' # alternative 3' splice site (A3SS) event
#' event <- read.table(text = "
#'  chr1 A3SS gene 15796 16765 . - .
#'  chr1 A3SS mRNA 15796 16765 . - .
#'  chr1 A3SS exon 15796 15947 . - .
#'  chr1 A3SS exon 16607 16765 . - .
#'  chr1 A3SS mRNA 15796 16765 . - .
#'  chr1 A3SS exon 15796 15942 . - .
#'  chr1 A3SS exon 16607 16765 . - .")
#' parseMisoA3SS(event)
parseMisoA3SS <- function(event, strand) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", rep("exon", 2), "mRNA", rep("exon", 2))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get the first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "A3SS",
                                        chromosome = chr,
                                        strand = strand)
        
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("C1.start", "C1.end")] <- 
                event[index[plus] + 2, 4:5]
            parsed[plus, ][["C2.start"]] <-
                lapply(index[plus], function(i) event[i + c(3, 6), 4])
            parsed[plus, ]["C2.end"] <- event[index[plus] + 3, 5]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("C1.start", "C1.end")] <- 
                event[index[minus] + 3, 5:4]
            parsed[minus, ][["C2.start"]] <-
                lapply(index[minus], function(i) event[i + c(2, 5), 5])
            parsed[minus, ]["C2.end"] <- event[index[minus] + 2, 4]
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#'
#' # Tandem UTR event
#' event <- read.table(text = "
#'  chr19 TandemUTR gene  10663759  10664625  .  -  .
#'  chr19 TandemUTR mRNA  10663759  10664625  .  -  .
#'  chr19 TandemUTR exon  10663759  10664625  .  -  .
#'  chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#'  chr19 TandemUTR exon  10664223  10664625  .  -  .")
#' parseMisoTandemUTR(event)
parseMisoTandemUTR <- function(event, strand) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", "exon", "mRNA", rep("exon", 1))
    event <- getValidEvents(event, validator)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get the first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "TandemUTR",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ]["C2.start"] <- event[index[plus] + 2, 4]
            parsed[plus, ][["C2.end"]] <-
                lapply(index[plus], function(i) event[i + c(2, 4), 5])
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ]["C2.start"] <- event[index[minus] + 2, 5]
            parsed[minus, ][["C2.end"]] <-
                lapply(index[minus], function(i) event[i + c(2, 4), 4])
        }
        return(parsed)
    }
}

#' @rdname parseMisoSE
#' @examples
#'
#' # alternative first exon (AFE) event
#' event <- read.table(text = "
#'  chr12 AFE gene 57916659 57920171  .  +  .
#'  chr12 AFE mRNA 57919131 57920171  .  +  .
#'  chr12 AFE exon 57919131 57920171  .  +  .
#'  chr12 AFE mRNA 57916659 57918199  .  +  .
#'  chr12 AFE exon 57916659 57916794  .  +  .
#'  chr12 AFE exon 57917812 57917875  .  +  .
#'  chr12 AFE exon 57918063 57918199  .  +  .")
#' parseMisoAFE(event)
parseMisoAFE <- function(event) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", "exon", "mRNA", "exon")
    event <- getValidEvents(event, validator, isRunLengthValid = TRUE)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get the first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        next_index <- nrow(event) + 1
        if (length(index) > 1)
            next_index <- c(index[2:length(index)], next_index)
        
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Get mRNAs index
        mRNA <- which(event[[3]] == "mRNA")
        splitting <- split(mRNA, c(TRUE, FALSE))
        mRNA1 <- splitting$`TRUE`
        mRNA2 <- splitting$`FALSE`
        
        # Create a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "AFE",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("C1.start", 
                             "C1.end")] <- event[mRNA2-1, 4:5][plus, ]
            parsed[plus, ][c("A1.start",
                             "A1.end")] <- event[next_index-1, 4:5][plus, ]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("C1.start",
                              "C1.end")] <- event[mRNA1+1, 5:4][minus, ]
            parsed[minus, ][c("A1.start", 
                              "A1.end")] <- event[mRNA2+1, 5:4][minus, ]
        }
        return(parsed)
    }
    
#     # Remove mRNAs from different chromosomes and mark event
#     len <- nrow(event)
#     event <- remove_wrong_mRNA(event)
#     if (len < nrow(event))
#         parsed[["comment"]] <- "wrong mRNAs"
#     
#     # Get a list of each mRNA and respective exons
#     mRNA <- list_mRNA(event)
#     
#     # Remove mRNAs with the same exons and mark event
#     len <- length(mRNA)
#     mRNA <- remove_duplicated_mRNA(mRNA)
#     if (len < length(mRNA))
#         parsed[["comment"]] <- "duplicated mRNAs"
#     
#     if (length(mRNA) != 2) {
#         # Store (but don't parse) events with more than two mRNAs or only one mRNA
#         parsed[["comment"]] <- "unrecognized event"
#     } else if (strand == "+") {
#         mRNA1 <- mRNA[[1]]
#         exon1 <- mRNA1[nrow(mRNA1), 4:5]
#         mRNA2 <- mRNA[[2]]
#         exon2 <- mRNA2[nrow(mRNA2), 4:5]
#         # Check if the most downstream exons are equal in both mRNAs
#         if (all(exon1 == exon2)) {
#             # Save positions of shared (constitutive) exon and alternative
#             # first exons
#             parsed[c("C1.start", "C1.end")] <- mRNA1[nrow(mRNA2) - 1, 4:5]
#             parsed[c("A1.start", "A1.end")] <- mRNA1[nrow(mRNA2) - 1, 4:5]
#             parsed[c("C2.start", "C2.end")] <- exon1
#         } else {
#             # Save the positions of the downstream exons
#             parsed[c("C1.start", "C1.end")] <- exon1
#             parsed[c("A1.start", "A1.end")] <- exon2
#         }
#     } else if (strand == "-") {
#         mRNA1 <- mRNA[[1]]
#         exon1 <- mRNA1[2, 5:4]
#         mRNA2 <- mRNA[[2]]
#         exon2 <- mRNA2[2, 5:4]
#         # Check if the most downstream exons are equal in both mRNAs
#         if (all(exon1 == exon2)) {
#             # Save positions of shared (constitutive) exon and alternative
#             # first exons
#             parsed[c("C1.start", "C1.end")] <- mRNA1[3, 5:4]
#             parsed[c("A1.start", "A1.end")] <- mRNA2[3, 5:4]
#             parsed[c("C2.start", "C2.end")] <- exon1
#         } else {
#             # Save the positions of the downstream exons
#             parsed[c("C1.start", "C1.end")] <- exon1
#             parsed[c("A1.start", "A1.end")] <- exon2
#         }
#     }
#     return(parsed)
}

#' @rdname parseMisoSE
#' @examples
#'
#' # alternative last exon (ALE) event
#' event <- read.table(text = "
#'  chr6 ALE gene 30620579 30822593  .  +  .
#'  chr6 ALE mRNA 30822190 30822593  .  +  .
#'  chr6 ALE exon 30822190 30822593  .  +  .
#'  chr6 ALE mRNA 30620579 30620982  .  +  .
#'  chr6 ALE exon 30620579 30620982  .  +  .")
#' parseMisoALE(event)
parseMisoALE <- function(event, strand) {
    # Filter out events that aren't valid
    validator <- c("gene", "mRNA", "exon", "mRNA", "exon")
    event <- getValidEvents(event, validator, isRunLengthValid = TRUE)
    
    # If there are valid events
    if (!is.null(event)) {
        # Get the first index, chromosome and strand of valids events
        index <- which(event[[3]] == "gene")
        next_index <- nrow(event) + 1
        if (length(index) > 1)
            next_index <- c(index[2:length(index)], next_index)
        
        chr <- as.character(event[index, 1])
        strand <- as.character(event[index, 7])
        
        # Get mRNAs index
        mRNA <- which(event[[3]] == "mRNA")
        splitting <- split(mRNA, c(TRUE, FALSE))
        mRNA1 <- splitting$`TRUE`
        mRNA2 <- splitting$`FALSE`
        
        # Creates a data frame of parsed junctions filled with NAs
        parsed <- createFilledJunctions(length(index),
                                        program = "MISO",
                                        event.type = "ALE",
                                        chromosome = chr,
                                        strand = strand)
        plus <- strand == "+"
        # Plus strand
        if (nrow(event[plus, ]) > 0) {
            parsed[plus, ][c("A1.start",
                             "A1.end")] <- event[mRNA1+1, 4:5][plus, ]
            parsed[plus, ][c("C2.start",
                             "C2.end")] <- event[mRNA2+1, 4:5][plus, ]
        }
        # Minus strand
        minus <- !plus
        if (nrow(event[minus, ]) > 0) {
            parsed[minus, ][c("A1.start", 
                              "A1.end")] <- event[mRNA2-1, 5:4][minus, ]
            parsed[minus, ][c("C2.start", 
                              "C2.end")] <- event[next_index-1, 5:4][minus, ]
        }
        return(parsed)
    }
    
    
#     # Creates a data frame of parsed junctions filled with NAs
#     parsed <- createFilledJunctions(1)
#     
#     len <- nrow(event)
#     if (len == 5) {
#         # Most ALE events have length of 5, so let's avoid wasting time
#         if (strand == "+") { 
#             parsed[c("A1.start", "A1.end")] <- event[3, 4:5]
#             parsed[c("C2.start", "C2.end")] <- event[5, 4:5]
#         } else if (strand == "-") {
#             parsed[c("A1.start", "A1.end")] <- event[5, 5:4]
#             parsed[c("C2.start", "C2.end")] <- event[3, 5:4]
#         }   
#     } else if (len < 5) {
#         # A GFF3 valid event needs at least 5 lines to be described
#         parsed[["comment"]] <- "unrecognized event"
#     } else {
#         event <- remove_wrong_mRNA(event)
#         if (len < nrow(event)) {
#             parsed[["comment"]] <- "wrong mRNAs"
#         }
#         
#         # Get a list of each mRNA and respective exons
#         mRNA <- list_mRNA(event)
#         
#         # Remove mRNAs with the same exons and mark event
#         len <- length(mRNA)
#         mRNA <- remove_duplicated_mRNA(mRNA)
#         if (len < length(mRNA))
#             parsed[["comment"]] <- "duplicated mRNAs"
#         
#         if (length(mRNA) != 2) {
#             # Store (but don't parse) events with more than two mRNAs or only one mRNA
#             parsed[["comment"]] <- "unrecognized event"
#         } else if (strand == "+") {
#             mRNA1 <- mRNA[[1]]
#             exon1 <- mRNA1[2, 4:5]
#             mRNA2 <- mRNA[[2]]
#             exon2 <- mRNA2[2, 4:5]
#             # Check if the most downstream exons are equal in both mRNAs
#             if (all(exon1 == exon2)) {
#                 parsed[c("C1.start", "C1.end")] <- exon1
#                 parsed[c("A1.start", "A1.end")] <- mRNA1[3, 4:5]
#                 parsed[c("C2.start", "C2.end")] <- mRNA2[3, 4:5]
#             } else {
#                 parsed[c("A1.start", "A1.end")] <- exon1
#                 parsed[c("C2.start", "C2.end")] <- exon2
#             }
#         } else if (strand == "-") {
#             mRNA1 <- mRNA[[1]]
#             exon1 <- mRNA1[nrow(mRNA1), 5:4]
#             mRNA2 <- mRNA[[2]]
#             exon2 <- mRNA2[nrow(mRNA2), 5:4]
#             # Check if the most downstream exons are equal in both mRNAs
#             if (all(exon1 == exon2)) {
#                 parsed[c("C1.start", "C1.end")] <- exon1
#                 parsed[c("A1.start", "A1.end")] <- mRNA1[nrow(mRNA2) - 1, 5:4]
#                 parsed[c("C2.start", "C2.end")] <- mRNA2[nrow(mRNA2) - 1, 5:4]
#             } else {
#                 parsed[c("A1.start", "A1.end")] <- exon1
#                 parsed[c("C2.start", "C2.end")] <- exon2
#             }
#         }
#     }
#     return(parsed)
}

#' Get a list of each mRNA and respective exons
#' 
#' From a given alternative splicing event, put each mRNA and respective exons
#' in a list to be returned.
#' 
#' @inheritParams parseMisoEvent
#' 
#' @return List of mRNAs and respective exons
#' @export
list_mRNA <- function(event) {
    mRNA_index <- which(event[ , 3] == "mRNA")
    
    next_index <- nrow(event)
    if (length(mRNA_index) > 1)
        next_index <- c(mRNA_index[2:length(mRNA_index)], next_index + 1)
    
    # Get each mRNA and respective exons as a new element of a list
    mRNA <- lapply(1:(length(mRNA_index)), getDataRows,
                   event, mRNA_index, next_index - 1)
    return(mRNA)
}

#' Clear mRNAs with the same exons from a given list of mRNAs and respective
#' exons
#' 
#' For a given list where each element contains one mRNA and respective exons,
#' compare if there are identical elements and only return each unique element
#' once.
#' 
#' @param mRNA List of mRNAs and respective exons
#' 
#' @importFrom fastmatch fmatch
#' 
#' @return Non-redundant list of mRNAs and respective exons
#' @export
remove_duplicated_mRNA <- function (mRNA) {
    # Get first occurence of each mRNA and remove duplicated index
    uniq <- unique(fmatch(mRNA, mRNA))
    
    # Return a non-redundant list of mRNAs
    return(mRNA[uniq])
}

#' Remove wrong mRNAs of a given alternative splicing event
#'
#' For a given event, checks if each mRNA and respective exons are located in 
#' the same chromosome and in the same region as the alternative splicing event
#' (this information is in the first row of the data frame)
#'
#' @inheritParams parseMisoEvent
#' 
#' @return Data frame with incorrect rows removed
#' @export
#' 
#' @examples
#' # skipping exon (SE) event with incorrect mRNAs
#' event <- read.table(text = "
#'   chr1 SE gene 16854	18061	. - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17055 . - .
#'   chr1 SE exon 17233 17742 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr1 SE mRNA 16854 18061 . - .
#'   chr1 SE exon 16854 17955 . - .
#'   chr1 SE exon 17915 18061 . - .
#'   chr1 SE mRNA 10 25 . - .
#'   chr1 SE exon 10 16 . - .
#'   chr1 SE exon 20 25 . - .
#'   chr7 SE mRNA 16854 18061 . - .
#'   chr7 SE exon 16854 17955 . - .")
#' remove_wrong_mRNA(event)
remove_wrong_mRNA <- function(event) {
    # Clear mRNA and exons with different chromosome identifier
    chr      <- event[1, 1]
    same_chr <- event[ , 1] == chr
    
    # Clear mRNA and exons outside event boundaries
    start  <- event[1, 4]
    end    <- event[1, 5]
    inside <- event[ , 4] >= start & event[ , 5] <= end
    
    return(event[same_chr & inside, ])
}