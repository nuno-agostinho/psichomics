library(plyr)
library(dplyr)

source("R/matsEvents.R")
source("R/misoEvents.R")
source("R/vastToolsEvents.R")
source("R/suppaEvents.R")

#' Creates a template of alternative splicing junctions
#' 
#' @param nrow Integer: Number of rows
#' @param program Character: Program used to get the junctions
#' @param event.type Character: Event type of the respective events
#' @param chromosome Character: Chromosome of the junctions
#' @param strand Strand: positive ("+") or negative ("-") strand of the event
#' 
#' @return A data frame with the junctions coordinate names pre-filled with NAs
#' @export
#' 
#' @examples
#' createJunctionsTemplate(nrow = 8)
createJunctionsTemplate <- function(nrow, program = character(0),
                                    event.type = character(0),
                                    chromosome = character(0),
                                    strand = character(0)) {
    parsed <- as.data.frame(matrix(NA, nrow = nrow, ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")
    
    if (length(program) > 0)    parsed[["Program"]] <- "MISO"
    if (length(event.type) > 0) parsed[["Event.type"]] <- event.type
    if (length(chromosome) > 0) parsed[["Chromosome"]] <- chromosome
    if (length(strand) > 0)     parsed[["Strand"]] <- strand
    return(parsed)
}

getMisoAnnotation <- function() {
    types <- c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI", "TandemUTR")
    typesFile <- paste0("/genedata/Resources/Annotations/MISO/hg19/", types,
                        ".hg19.gff3")
    miso.hg19 <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                        comment.char="#", header=FALSE)
    ## TODO: ALE events are baldy formatted, they have two consecutive gene
    ## lines... remove them for now
    miso.hg19[[3]] <- miso.hg19[[3]][-c(49507, 49508), ]
    
    misoEvents <- lapply(miso.hg19, parseMisoEvent)
    misoEvents <- plyr::rbind.fill(misoEvents)
    return(misoEvents)
}

getSuppaAnnotation <- function() {
    types <- c("SE", "AF", "AL", "MX", "A5", "A3", "RI")
    typesFile <- paste0("~/Documents/psi_calculation/suppa/suppaEvents/hg19_", 
                        types, ".ioe")
    suppa.hg19 <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                         comment.char="#", header=TRUE)
    
    eventsID <- lapply(suppa.hg19, "[[", "event_id")
    suppaEvents <- lapply(eventsID, parseSuppaEvent)
    suppaEvents <- plyr::rbind.fill(suppaEvents)
    return(suppaEvents)
}

getMatsAnnotation <- function() {
    types <- c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI")
    typesFile <- paste("~/Documents/psi_calculation/mats_out/ASEvents/fromGTF",
                       c(types, paste0("novelEvents.", types)), "txt",
                       sep = ".")
    mats.hg19 <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                        comment.char="#", header=TRUE)
    
    matsEvents <- lapply(seq_along(mats.hg19),
                         function(i) {
                             type <- rep(types, 2)[i]
                             annotation <- mats.hg19[[i]]
                             if (nrow(annotation) > 0)
                                 return(parseMatsEvent(annotation, type))
                         })
    matsEvents <- plyr::rbind.fill(matsEvents)
    
    # Sum 1 position to the start/end of MATS events (depending on the strand)
    matsNames <- names(matsEvents)
    plus <- matsEvents$Strand == "+"
    # Plus
    start <- matsNames[grep(".start", matsNames)]
    matsEvents[plus, start] <- matsEvents[plus, start] + 1
    # Minus
    end <- matsNames[grep(".end", matsNames)]
    matsEvents[!plus, end] <- matsEvents[!plus, end] + 1
    
    return(matsEvents)
}

getVastToolsAnnotation <- function() {
    types <- c("ALT3", "ALT5", "COMBI", "IR", "MERGE3m", "MIC",
               rep(c("EXSK", "MULTI"), 1))
    typesFile <- sprintf(
        "/genedata/Resources/Software/vast-tools/VASTDB/Hsa/TEMPLATES/Hsa.%s.Template%s.txt",
        types, c(rep("", 6), rep(".2", 2))#, rep(".2", 2))
    )
    
    vastTools.hg19 <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                             comment.char="#", header=TRUE)
    
    vastToolsEvents <- lapply(seq_along(vastTools.hg19),
                              function(i) {
                                  type <- types[i]
                                  print(type)
                                  annotation <- vastTools.hg19[[i]]
                                  if (nrow(annotation) > 0)
                                      return(parseVastToolsEvent(annotation))
                              })
    vastToolsEvents <- plyr::rbind.fill(vastToolsEvents)
    return(vastToolsEvents)
}

#' Returns the coordinates of interest for a given event type
#' @export
getCoordinates <- function(type) {
    switch(type,
           "SE"   = c("C1.end", "A1.start", "A1.end", "C2.start"),
           "A3SS" = c("C1.end", "C2.start", "A1.start"),
           "A5SS" = c("C1.end", "C2.start", "A1.end"),
           "AFE"  = c("C1.start", "C1.end", "A1.start", "A1.end", "C2.start"),
           "ALE"  = c("A1.start", "A1.end", "C2.start", "C2.end"),
           "RI"   = c("C1.start", "C1.end", "C2.start", "C2.end"),
           "MXE"  = c("C1.end", "A1.start", "A1.end",
                      "A2.start", "A2.end", "C2.start"), 
           "TandemUTR" = c("C2.start", "C2.end", "A1.end"))
}

#' Get the annotation for all event types
#' @export
getAnnotation <- function() {
    events <- list("miso" = getMisoAnnotation(),
                   "mats" = getMatsAnnotation(),
                   "vast-tools" = getVastToolsAnnotation(),
                   "suppa" = getSuppaAnnotation())
    print("Annotation retrieved")
    
    # Remove the "chr" prefix from the chromosome field
    for (each in seq_along(events)) {
        chr <- grepl("chr", events[[each]]$Chromosome)
        events[[each]]$Chromosome[chr] <-
            gsub("chr", "", events[[each]]$Chromosome[chr])
    }
    
    events <- rbind.fill(events)
    events <- dlply(events, .(Event.type))
    events <- lapply(events, dlply, .(Program))
    return(events)
}

#' Convert a column to numeric if possible and ignore given columns composed
#' of lists
#' 
#' @export
#' 
#' @examples
#' getNumerics(vast, by = c("Strand", "C1.end", "A1.end", "A1.start"),
#'     toNumeric = c(F,T,T,T))
getNumerics <- function(table, by = NULL, toNumeric = FALSE) {
    # Check which elements are lists of specified length
    bool <- TRUE
    for (each in by)
        bool <- bool & vapply(table[[each]], length, integer(1)) == 1
    
    table <- table[bool, ]
    # Convert elements to numeric
    conv <- by[toNumeric]
    table[conv] <- as.numeric(as.character(unlist(table[conv])))
    return(table)
}

#' Full outer join all given annotation based on select columns
#' @export
joinAnnotation <- function(annotation) {
    types <- names(annotation)
    joint <- lapply(types, function(type, annotation) {
        print(type)
        # Create vector with comparable columns
        id <- c("Strand", "Chromosome", "Event.type")
        by <- c(id, getCoordinates(type))
        toNumeric <- !by %in% id
        
        # Convert given columns to numeric if possible
        tables <- lapply(annotation[[type]], getNumerics, by, toNumeric)
        
        # Make the names of non-comparable columns distinct
        cols <- lapply(names(tables), function(k) {
            ns <- names(tables[[k]])
            inBy <- ns %in% by
            ifelse(inBy, ns, paste(k, ns, sep="."))
        })
        
        # Full join all the tables
        res <- Reduce(function(x, y) full_join(x, y, by), tables)
        names(res) <- unique(unlist(cols))
        
        # Remove equal rows
        res <- distinct_(res)
        return(res)
    }, annotation)
    names(joint) <- types
    return(joint)
}

#' Write the annotation of an event type to a file
#' 
#' @param join List of lists of data frame
#' @param eventType Character: type of event
#' @param filename Character: path to the annotation file
#' 
#' @return Invisible TRUE if everything's okay
#' @export
writeAnnotation <- function(jointEvents, eventType,
                            filename = paste0("data/annotation_",
                                              eventType, ".txt")) {
    res <- jointEvents[[eventType]]
    # Show only the columns Chromosome, Strand and coordinates of interest
    by <- c("Chromosome", "Strand", getCoordinates(eventType))
    res <- subset(res, select = by)
    
    # Order by chromosome and coordinates
    orderBy <- lapply(c(1, 3:ncol(res)), function(x) return(res[[x]]))
    res <- res[do.call(order, orderBy), ]
    
    write.table(res, file = filename, quote = FALSE, row.names = FALSE)
    return(invisible(TRUE))
}

#' Read the annotation of an event type from a file
#' 
#' @inheritParams writeAnnotation
#' 
#' @return Data frame with the annotation
#' @export
readAnnotation <- function(eventType, filename = paste0("data/annotation_",
                                                        eventType, ".txt")) {
    res <- read.table(filename, header = TRUE)
    return(res)
}

#' Compare the number of events from the different programs in a Venn diagram
#' 
#' @param join List of lists of data frame
#' @param eventType Character: type of event
#' 
#' @return Venn diagram
#' 
#' @export
vennEvents <- function(join, eventType) {
    join <- join[[eventType]]
    
    programs <- join[grep("Program", names(join))]
    nas <- !is.na(programs)
    nas <- ifelse(nas, row(nas), NA)
    p <- lapply(1:ncol(nas), function(col) nas[!is.na(nas[ , col]), col])
    names(p) <- sapply(programs, function(x) unique(x[!is.na(x)]))
    gplots::venn(p)
}

#' String used to search for matches in a junction quantification file
#' @param chr Character: chromosome
#' @param strand Character: strand
#' @param junc5 Integer: 5' end junction
#' @param junc3 Integer: 3' end junction
#' 
#' @return Formatted character string
junctionString <- function(chr, strand, junc5, junc3) {
    plus <- strand == "+"
    res <- sprintf("chr%s:%s:%s,chr%s:%s:%s",
                   chr, ifelse(plus, junc5, junc3), strand,
                   chr, ifelse(plus, junc3, junc5), strand)
    return(res)
}

#' Calculate inclusion levels using alternative splicing event annotation and
#' junction quantification for many samples
#' 
#' @param eventType Character: type of the alternative event to calculate
#' @param junctionQuant Data.frame: junction quantification with samples as
#' columns and junctions as rows
#' @param annotation Data.frame: alternative splicing annotation related to
#' event type
#' 
#' @importFrom fastmatch fmatch
#' @return Matrix with inclusion levels
#' @export
calculateInclusionLevels <- function(eventType, junctionQuant, annotation) {
    # TODO(NunoA): Really? This should be unique by now
    annotation <- unique(annotation)
    chr <- annotation$Chromosome
    strand <- annotation$Strand
    
    if (eventType == "SE") {
        c1e <- annotation$C1.end
        a1s <- annotation$A1.start
        a1e <- annotation$A1.end
        c2s <- annotation$C2.start
        
        # Convert junction quantification files to matrix
        ## TODO(NunoA): shouldn't this be processed already?
        junctionQuant <- unique(junctionQuant)
        row.names(junctionQuant) <- junctionQuant$Hybridization.REF    
        junctions <- data.matrix(junctionQuant[, 2:ncol(junctionQuant)])
        
        # Get specific junctions
        incAstr <- junctionString(chr, strand, c1e, a1s)
        incBstr <- junctionString(chr, strand, a1e, c2s)
        exclstr <- junctionString(chr, strand, c1e, c2s)
        incA <- junctions[fmatch(incAstr, rownames(junctions)), ]
        incB <- junctions[fmatch(incBstr, rownames(junctions)), ]
        excl <- junctions[fmatch(exclstr, rownames(junctions)), ]
        
        # Calculate inclusion levels
        inc <- (incA + incB) / 2
        psi <- inc/(excl + inc)
        
        # Clear rows with nothing but NAs
        psi <- psi[complete.cases(psi), ]
        return(psi)
    }
}