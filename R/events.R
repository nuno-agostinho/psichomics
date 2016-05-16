#' @include matsEvents.R
#' @include misoEvents.R
#' @include vastToolsEvents.R
#' @include suppaEvents.R
NULL

#' Creates a template of alternative splicing junctions
#' 
#' @param nrow Integer: Number of rows
#' @param program Character: Program used to get the junctions
#' @param event.type Character: Event type of the respective events
#' @param chromosome Character: Chromosome of the junctions
#' @param strand Character: positive ("+") or negative ("-") strand of the event
#' @param id Character: events' ID
#' 
#' @return A data frame with the junctions coordinate names pre-filled with NAs
#' @export
#' 
#' @examples
#' createJunctionsTemplate(nrow = 8)
createJunctionsTemplate <- function(nrow, program = character(0),
                                    event.type = character(0),
                                    chromosome = character(0),
                                    strand = character(0),
                                    id = character(0)) {
    ## TODO(NunoA): only accept a "+" or a "-" strand
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
    if (length(id) > 0)         parsed[["Event.ID"]] <- id
    return(parsed)
}

#' @importFrom utils read.delim
getMisoAnnotation <- function() {
    types <- c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI", "TandemUTR")
    typesFile <- paste0("/genedata/Resources/Annotations/MISO/hg19/", types,
                        ".hg19.gff3")
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=FALSE)
    
    ## TODO: ALE events are baldy formatted, they have two consecutive gene
    ## lines... remove them for now
    annot[[3]] <- annot[[3]][-c(49507, 49508), ]
    return(annot)
}

#' @importFrom plyr rbind.fill
parseMisoAnnotation <- function(annot) {
    events <- lapply(annot, parseMisoEvent)
    events <- rbind.fill(events)
    return(events)
}

#' @importFrom utils read.delim
getSuppaAnnotation <- function() {
    types <- c("SE", "AF", "AL", "MX", "A5", "A3", "RI")
    typesFile <- paste0("~/Documents/psi_calculation/suppa/suppaEvents/hg19_", 
                        types, ".ioe")
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)
    return(annot)
}

#' @importFrom plyr rbind.fill
parseSuppaAnnotation <- function(annot) {
    eventsID <- lapply(annot, "[[", "event_id")
    events <- lapply(eventsID, parseSuppaEvent)
    events <- rbind.fill(events)
    return(events)
}

#' @importFrom utils read.delim
getMatsAnnotation <- function() {
    types <- c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI")
    typesFile <- paste("~/Documents/psi_calculation/mats_out/ASEvents/fromGTF",
                       c(types, paste0("novelEvents.", types)), "txt",
                       sep = ".")
    names(typesFile) <- rep(types, 2)
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)
    return(annot)
}

#' @importFrom plyr rbind.fill
parseMatsAnnotation <- function(annot) {
    types <- names(annot)
    events <- lapply(seq_along(annot), function(i)
        if (nrow(annot[[i]]) > 0) 
            return(parseMatsEvent(annot[[i]], types[[i]])))
    events <- rbind.fill(events)
    
    # Sum 1 position to the start/end of MATS events (depending on the strand)
    matsNames <- names(events)
    plus <- events$Strand == "+"
    # Plus
    start <- matsNames[grep(".start", matsNames)]
    events[plus, start] <- events[plus, start] + 1
    # Minus
    end <- matsNames[grep(".end", matsNames)]
    events[!plus, end] <- events[!plus, end] + 1
    
    return(events)
}

#' @importFrom utils read.delim
getVastToolsAnnotation <- function() {
    types <- c("ALT3", "ALT5", "COMBI", "IR", "MERGE3m", "MIC",
               rep(c("EXSK", "MULTI"), 1))
    typesFile <- sprintf(
        "/genedata/Resources/Software/vast-tools/VASTDB/Hsa/TEMPLATES/Hsa.%s.Template%s.txt",
        types, c(rep("", 6), rep(".2", 2))#, rep(".2", 2))
    )
    names(typesFile) <- types
    
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)
    return(annot)
}

#' @importFrom plyr rbind.fill
parseVastToolsAnnotation <- function(annot) {
    types <- names(annot)
    events <- lapply(seq_along(annot),
                     function(i) {
                         print(types[i])
                         a <- annot[[i]]
                         if (nrow(a) > 0)
                             return(parseVastToolsEvent(a))
                     })
    events <- rbind.fill(events)
    return(events)
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
           "TandemUTR" = c("C1.start", "C1.end", "A1.end"))
}

#' Get the annotation for all event types
#' @export
getParsedAnnotation <- function() {
    print("Retrieving MISO annotation...")
    annot <- getMisoAnnotation()
    print("Parsing MISO annotation...")
    miso <- parseMisoAnnotation(annot)
    
    print("Retrieving SUPPA annotation...")
    annot <- getSuppaAnnotation()
    print("Parsing SUPPA annotation...")
    suppa <- parseSuppaAnnotation(annot)
    
    print("Retrieving VAST-TOOLS annotation...")
    annot <- getVastToolsAnnotation()
    print("Parsing VAST-TOOLS annotation...")
    vast <- parseVastToolsAnnotation(annot)
    
    print("Retrieving MATS annotation...")
    annot <- getMatsAnnotation()
    print("Parsing MATS annotation...")
    mats <- parseMatsAnnotation(annot)
    
    events <- list(
        "miso" = miso, "mats" = mats, "vast-tools" = vast, "suppa" = suppa)
    
    # Remove the "chr" prefix from the chromosome field
    print("Standarising chromosome field")
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
#' @examples
#' event <- read.table(text = "ABC123 + 250 300 350
#'                             DEF456 - 900 800 700")
#' names(event) <- c("Event ID", "Strand", "C1.end", "A1.end", "A1.start")
#' 
#' # Let's change one column to character
#' event[ , "C1.end"] <- as.character(event[ , "C1.end"])
#' is.character(event[ , "C1.end"])
#' 
#' event <- getNumerics(event, by = c("Strand", "C1.end", "A1.end", "A1.start"),
#'                      toNumeric = c(FALSE, TRUE, TRUE, TRUE))
#' # Let's check if the same column is now integer
#' is.numeric(event[ , "C1.end"])
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
#' @importFrom dplyr full_join
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
        return(unique(res))
    }, annotation)
    names(joint) <- types
    return(joint)
}

#' Write the annotation of an event type to a file
#' 
#' @param jointEvents List of lists of data frame
#' @param eventType Character: type of event
#' @param filename Character: path to the annotation file
#' @param showID Boolean: show the events' ID (FALSE by default)
#' 
#' @importFrom utils write.table
#' 
#' @return Invisible TRUE if everything's okay
#' @export
writeAnnotation <- function(jointEvents, eventType,
                            filename = paste0("data/annotation_",
                                              eventType, ".txt"),
                            showID = FALSE) {
    res <- jointEvents[[eventType]]
    # Show the columns Chromosome, Strand and coordinates of interest
    by <- c("Chromosome", "Strand", getCoordinates(eventType))
    ord <- 0
    
    # Show the events' ID if desired
    if (showID) {
        cols <- grep("Event.ID", names(res), value = TRUE)
        by <- c(cols, by)
        ord <- length(cols)
    }
    res <- subset(res, select = by)
    
    ## TODO(NunoA): clean this mess
    # Order by chromosome and coordinates
    orderBy <- lapply(c(1 + ord, (3 + ord):ncol(res)),
                      function(x) return(res[[x]]))
    res <- res[do.call(order, orderBy), ]
    
    res <- unique(res)
    write.table(res, file = filename, quote = FALSE, row.names = FALSE,
                sep = "\t")
    return(invisible(TRUE))
}

#' Read the annotation of an event type from a file
#' 
#' @inheritParams writeAnnotation
#' 
#' @importFrom utils read.table
#' 
#' @return Data frame with the annotation
#' @export
readAnnotation <- function(eventType, filename = paste0("data/annotation_",
                                                        eventType, ".txt")) {
    res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
    return(res)
}

#' Compare the number of events from the different programs in a Venn diagram
#' 
#' @param join List of lists of data frame
#' @param eventType Character: type of event
#' 
#' @importFrom gplots venn
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
    venn(p)
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
    first <- ifelse(plus, junc5, junc3)
    last <- ifelse(plus, junc3, junc5)
    res <- sprintf("chr%s:%s:%s,chr%s:%s:%s",
                   chr, first, strand, chr, last, strand)
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
    chr <- annotation$Chromosome
    strand <- annotation$Strand
    
    if (eventType == "SE") {
        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand,
                                  annotation$C1.end, annotation$A1.start)
        incBstr <- junctionString(chr, strand,
                                  annotation$A1.end, annotation$C2.start)
        exclstr <- junctionString(chr, strand, 
                                  annotation$C1.end, annotation$C2.start)
        
        # Get specific junctions
        fmatch <- fmatch
        coords <- rownames(junctionQuant)
        incA <- junctionQuant[fmatch(incAstr, coords), ]
        incB <- junctionQuant[fmatch(incBstr, coords), ]
        excl <- junctionQuant[fmatch(exclstr, coords), ]
        
        # Calculate inclusion levels
        inc <- (incA + incB) / 2
        psi <- inc/(excl + inc)
        rownames(psi) <- sprintf("%s_%s_%s_%s_%s_%s", chr, strand,
                                 annotation$C1.end, annotation$A1.start,
                                 annotation$A1.end, annotation$C2.start)
        
        # Clear rows with nothing but NAs
        naRows <- rowSums(!is.na(psi)) == 0
        return(psi[!naRows, ])
    }
}