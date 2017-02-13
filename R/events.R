#' @include events_mats.R
#' @include events_miso.R
#' @include events_vastTools.R
#' @include events_suppa.R
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
#' 
#' @examples
#' psichomics:::createJunctionsTemplate(nrow = 8)
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

#' Returns the coordinates of interest for a given event type
#' 
#' @param type Character: alternative splicing event type
#' @param sorting Boolean: get coordinates used for sorting and comparison 
#' between different programs? FALSE by default
#' 
#' @return Coordinates of interest according to the alternative splicing event
#' type
getSplicingEventCoordinates <- function(type, sorting=FALSE) {
    coords <- switch(type,
                     "SE"   = c("C1.end", "A1.start", "A1.end", "C2.start"),
                     "A3SS" = c("C1.end", "C2.start", "A1.start"),
                     "A5SS" = c("C1.end", "C2.start", "A1.end"),
                     "AFE"  = c("C1.start", "C1.end", "A1.start", "A1.end"),
                     "ALE"  = c("A1.start", "A1.end", "C2.start", "C2.end"),
                     "RI"   = c("C1.start", "C1.end", "C2.start", "C2.end"),
                     "MXE"  = c("C1.end", "A1.start", "A1.end",
                                "A2.start", "A2.end", "C2.start"),
                     "TandemUTR" = c("C1.start", "C1.end", "A1.end"))
    
    if (sorting) {
        coords <- switch(type,
                         "A3SS" = c("C2.start", "A1.start"),
                         "A5SS" = c("C1.end", "A1.end"),
                         "AFE"  = c("A1.start", "A1.end", "C1.start", "C1.end"),
                         "ALE"  = c("A1.start", "A1.end", "C2.start", "C2.end"),
                         "MXE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
                         "TandemUTR" = c("A1.end", "C1.end"))
    }
    return(coords)
}

#' Convert a column to numeric if possible and ignore given columns composed
#' of lists
#'
#' @param table Data matrix: table
#' @param by Character: column names of interest
#' @param toNumeric Boolean: which columns to convert to numeric (FALSE by
#' default)
#'
#' @return Processed data matrix
#' @examples
#' event <- read.table(text = "ABC123 + 250 300 350
#'                             DEF456 - 900 800 700")
#' names(event) <- c("Event ID", "Strand", "C1.end", "A1.end", "A1.start")
#'
#' # Let's change one column to character
#' event[ , "C1.end"] <- as.character(event[ , "C1.end"])
#' is.character(event[ , "C1.end"])
#'
#' event <- psichomics:::getNumerics(event, by = c("Strand", "C1.end", "A1.end",
#'                                   "A1.start"),
#'                                   toNumeric = c(FALSE, TRUE, TRUE, TRUE))
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

#' Full outer join all given events based on select columns
#' 
#' @param events Data frame or matrix: alternative splicing events
#' @param types Character: alternative splicing types
#' 
#' @return List of events joined by alternative splicing event type
joinEventsPerType <- function(events, types) {
    if (missing(types)) types <- names(events)
    joint <- lapply(types, function(type, events) {
        cat(type, fill=TRUE)
        # Create vector with comparable columns
        id <- c("Strand", "Chromosome", "Event.type")
        by <- c(id, getSplicingEventCoordinates(type))
        toNumeric <- !by %in% id
        
        # Convert given columns to numeric if possible
        tables <- lapply(events[[type]], getNumerics, by, toNumeric)
        
        # Make the names of non-comparable columns distinct
        cols <- lapply(names(tables), function(k) {
            ns <- names(tables[[k]])
            inBy <- ns %in% by
            ifelse(inBy, ns, paste(k, ns, sep="."))
        })
        
        # Full join all the tables
        res <- Reduce(function(x, y) dplyr::full_join(x, y, by), tables)
        names(res) <- unique(unlist(cols))
        
        # Remove equal rows
        return(unique(res))
    }, events)
    names(joint) <- types
    return(joint)
}

#' Prepare annotation from alternative splicing events
#'
#' In case more than one data frame with alternative splicing events is given,
#' the events are cross-referenced according to the chromosome, strand and 
#' relevant coordinates per event type (see details).
#'
#' @param ... Data frame(s) of alternative splicing events to include in the
#' annotation
#'
#' @details Events from two or more data frames are cross-referenced based on
#' each event's chromosome, strand and specific coordinates relevant for each
#' event type:
#' \itemize{
#'      \item Skipped exon: constitutive exon 1 end, alternative exon (start
#'      and end) and constitutive exon 2 start
#'      \item Mutually exclusive exon: constitutive exon 1 end, alternative exon
#'      1 and 2 (start and end) and constitutive exon 2 start
#'      \item Alternative 5' splice site: constitutive exon 1 end, alternative 
#'      exon 1 end and constitutive exon 2 start
#'      \item Alternative first exon: same as alternative 5' splice site
#'      \item Alternative 3' splice site: constitutive exon 1 end, alternative
#'      exon 1 start and constitutive exon 2 start
#'      \item Alternative last exon: same as alternative 3' splice site
#' }
#'
#' @note When cross-referencing events, gene information is discarded.
#'
#' @importFrom plyr rbind.fill dlply
#'
#' @return List of data frames with the annotation from different data frames
#' joined by event type
#' @export
#' @examples 
#' # Load sample files (SUPPA annotation)
#' folder <- "extdata/eventsAnnotSample/suppa_output/suppaEvents"
#' suppaOutput <- system.file(folder, package="psichomics")
#' 
#' # Parse and prepare SUPPA annotation
#' suppa <- parseSuppaAnnotation(suppaOutput)
#' annot <- prepareAnnotationFromEvents(suppa)
#' 
#' # Load sample files (rMATS annotation)
#' folder <- "extdata/eventsAnnotSample/mats_output/ASEvents/"
#' matsOutput <- system.file(folder, package="psichomics")
#' 
#' # Parse rMATS annotation and prepare combined annotation from rMATS and SUPPA
#' mats <- parseMatsAnnotation(matsOutput)
#' annot <- prepareAnnotationFromEvents(suppa, mats)
prepareAnnotationFromEvents <- function(...) {
    events <- list(...)
    if (!all(vapply(events, is, "ASevents", FUN.VALUE=logical(1))))
        warning("All variables should be an object of class ASevents")
    
    # Remove the "chr" prefix from the chromosome field
    for (each in seq_along(events)) {
        events[[each]]$Chromosome <- gsub("chr", "", events[[each]]$Chromosome)
    }
    
    # Organise splicing events by event type and then by program in a list of 
    # list of dataframes
    events <- rbind.fill(events)
    events <- dlply(events, "Event.type")
    events <- lapply(events, dlply, "Program")
    
    cat("Sorting coordinates...", fill=TRUE)
    events <- sortCoordinates(events)
    
    cat("Joining events per event type...", fill=TRUE)
    join <- joinEventsPerType(events)
    
    # If available, add 1st constitutive exon's end and 2nd constituve exon's 
    # start from SUPPA or rMATS to AFE and ALE events, respectively, as other 
    # programs may not state these coordinates
    suppaAFE <- join$AFE$SUPPA.C2.start
    matsAFE  <- join$AFE$MATS.C2.start
    AFE.C2.start <- as.numeric( ifelse(
        sapply(suppaAFE, is.null),
        ifelse(sapply(matsAFE, is.null), NA, unlist(matsAFE)),
        unlist(suppaAFE)))
    
    suppaALE <- join$ALE$SUPPA.C1.end
    matsALE  <- join$ALE$MATS.C1.end
    ALE.C1.end <- as.numeric( ifelse(
        sapply(suppaALE, is.null),
        ifelse(sapply(matsALE, is.null), NA, unlist(matsALE)),
        unlist(suppaALE)))
    
    # Organise columns
    annot <- lapply(names(join), function(i) {
        # Include genes if there is only one column
        geneCols <- grep(".Gene", names(join[[i]]), fixed=TRUE)
        if (length(geneCols) == 1) {
            gene <- "Gene"
            names(join[[i]])[geneCols] <- gene
            join[[i]][[gene]] <- as.list(join[[i]][[gene]])
        } else {
            gene <- NULL
        }
        
        cols <- c("Chromosome", "Strand", gene, getSplicingEventCoordinates(i),
                  grep("Event.ID", names(join[[i]]), value = TRUE))
        return(join[[i]][, cols])
    })
    names(annot) <- names(join)
    annot$AFE["C2.start"] <- AFE.C2.start
    annot$ALE["C1.end"]   <- ALE.C1.end
    events <- annot
    
    cat("Cleaning the annotation...", fill=TRUE)
    types <- c(SE="Skipped exon", MXE="Mutually exclusive exon",
               A3SS="Alternative 3' splice site", 
               A5SS="Alternative 5' splice site",
               AFE="Alternative first exon", ALE="Alternative last exon",
               RI="Retained intron", TandemUTR="Tandem UTR")
    
    for (type in names(types)) {
        if (!is.null(events[[type]]))
            events[[type]] <- cbind("Event type"=types[[type]], events[[type]])
    }
    events <- rbind.fill(events)
    
    coords <- c("C1.start"="Constitutive exon 1 start",
                "C1.end"="Constitutive exon 1 end",
                "A1.start"="Alternative exon 1 start",
                "A1.end"="Alternative exon 1 end",
                "A2.start"="Alternative exon 2 start",
                "A2.end"="Alternative exon 2 end",
                "C2.start"="Constitutive exon 2 start",
                "C2.end"="Constitutive exon 2 end")
    m <- match(names(coords), names(events))
    names(events)[m[!is.na(m)]] <- coords[!is.na(m)]
    if (is.null(events$Gene)) events$Gene <- NA
    eventId <- grep("Event.ID", names(events), value = TRUE)
    
    # Order rows by event type, chromosome and the first exons coordinates and
    # order columns
    ord <- order(events$`Event type`, events$Chromosome,
                 events$`Constitutive exon 1 start`,
                 events$`Constitutive exon 1 end`,
                 events$`Alternative exon 1 start`,
                 events$`Alternative exon 1 end`)
    events <- events[ord, c("Event type", "Chromosome", "Strand", "Gene",
                            coords[!is.na(m)], eventId)]
    
    # Organise by event type and remove columns with NAs only
    events <- split(events[-1], events$`Event type`)
    for (type in names(events)) {
        naCols <- apply(events[[type]], 2, function(col) all(is.na(col)))
        events[[type]] <- events[[type]][!naCols]
    }
    return(events)
}

#' Compare the number of events from the different programs in a Venn diagram
#' 
#' @param join List of lists of data frame
#' @param eventType Character: type of event
#' 
#' @return Venn diagrams for a given event type
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
#'
#' @param chr Character: chromosome
#' @param strand Character: strand
#' @param junc5 Integer: 5' end junction
#' @param junc3 Integer: 3' end junction
#' @param showStrand Boolean: include strand?
#' 
#' @return Formatted character string
junctionString <- function(chr, strand, junc5, junc3, showStrand) {
    plus <- strand == "+"
    first <- ifelse(plus, junc5, junc3)
    last <- ifelse(plus, junc3, junc5)
    if (showStrand)
        res <- sprintf("chr%s:%s:%s:%s", chr, first, last, strand)
    else
        res <- sprintf("chr%s:%s:%s", chr, first, last)
    return(res)
}

colsAsNumbers <- function(type, annotation) {
    # Create vector with comparable columns
    id <- c("Strand", "Chromosome", "Event.type")
    by <- c(id, getSplicingEventCoordinates(type))
    toNumeric <- !by %in% id
    
    # Convert given columns to numeric if possible
    tables <- lapply(annotation[[type]], getNumerics, by, toNumeric)
    return(tables)
}

#' Sort coordinates for some event types
#' 
#' Some programs sort the coordinates of specific event types differently. To
#' make them all comparable across programs, the coordinates are ordered by
#' increasing (plus strand) or descresing order (minus strand)
#' 
#' @param events List of data frames with alternative splicing events for a 
#' given program
#' 
#' @return List of data frames with alternative splicing events for a given
#' program
sortCoordinates <- function(events) {
    types <- names(events)
    for (type in types) {
        coord <- getSplicingEventCoordinates(type, sorting=TRUE)
        events[[type]] <- colsAsNumbers(type, events)
        if (!is.null(coord)) {
            for (program in names(events[[type]])) {
                print(paste(type, program))
                table <- events[[type]][[program]]
                plus <- table[["Strand"]] == "+"
                plusOrd <- apply(table[plus, coord], 1, sort)
                minusOrd <- apply(table[!plus, coord], 1, sort, decreasing=TRUE)
                if (length(plusOrd) > 0)
                    events[[type]][[program]][plus, coord] <- t(plusOrd)
                if (length(minusOrd) > 0)
                    events[[type]][[program]][!plus, coord] <- t(minusOrd)
            }
        }
    }
    return(events)
}

#' Calculate inclusion levels using alternative splicing event annotation and
#' junction quantification for many samples
#' 
#' @param eventType Character: type of the alternative event to calculate
#' @param junctionQuant Data.frame: junction quantification with samples as
#' columns and junctions as rows
#' @param annotation Data.frame: alternative splicing annotation related to
#' event type
#' @param minReads Integer: minimum of total reads required to consider the
#' quantification as valid (10 by default)
#' 
#' @importFrom fastmatch fmatch
#' @return Matrix with inclusion levels
calculateInclusionLevels <- function(eventType, junctionQuant, annotation,
                                     minReads = 10) {
    # Immediately return NULL if ALE and AFE events are missing coordinates
    if (eventType == "AFE" && is.null(annotation$`Constitutive exon 2 start`)) {
        return(NULL)
    } else if (eventType == "ALE" && 
               is.null(annotation$`Constitutive exon 1 end`)) {
        return(NULL)
    }
    
    if (is.null(annotation$Gene))
        geneCol <- NULL
    else
        geneCol <- "Gene"
    
    coords <- rownames(junctionQuant)
    showStrand <- any(grepl("\\+|\\-", coords))
    
    if (eventType == "SE") {
        # Remove duplicates based on columns used to create identifiers
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               "Constitutive exon 1 end",
                               "Alternative exon 1 start",
                               "Alternative exon 1 end",
                               "Constitutive exon 2 start", geneCol)
        chr <- annotation$Chromosome
        strand <- annotation$Strand
        
        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand,
                                  annotation$`Constitutive exon 1 end`,
                                  annotation$`Alternative exon 1 start`,
                                  showStrand)
        incBstr <- junctionString(chr, strand,
                                  annotation$`Alternative exon 1 end`,
                                  annotation$`Constitutive exon 2 start`,
                                  showStrand)
        exclstr <- junctionString(chr, strand, 
                                  annotation$`Constitutive exon 1 end`, 
                                  annotation$`Constitutive exon 2 start`,
                                  showStrand)
        
        # Get specific junction quantification
        coords <- rownames(junctionQuant)
        incA <- junctionQuant[fmatch(incAstr, coords), ]
        incB <- junctionQuant[fmatch(incBstr, coords), ]
        excl <- junctionQuant[fmatch(exclstr, coords), ]
        rm(incAstr, incBstr, exclstr)
        
        # Calculate inclusion levels
        inc <- (incA + incB) / 2
        rm(incA, incB)
        
        tot <- excl + inc
        rm(excl)
        if (nrow(tot) == 0)
            return(NULL)
        
        # Prepare presentation of multigenes
        multigene <- lapply(annotation$Gene, length) > 1
        gene <- annotation$Gene
        gene[multigene] <- lapply(gene[multigene], paste, collapse="/")
        
        eventNames <- paste(sep="_", eventType, chr, strand, 
                            annotation$`Constitutive exon 1 end`, 
                            annotation$`Alternative exon 1 start`, 
                            annotation$`Alternative exon 1 end`,
                            annotation$`Constitutive exon 2 start`, gene)
    } else if (eventType == "MXE") {
        # Remove duplicates based on columns used to create identifiers
        annotation <- uniqueBy(
            annotation, "Chromosome", "Strand", "Constitutive exon 1 end", 
            "Alternative exon 1 start", "Alternative exon 1 end",
            "Alternative exon 2 start", "Alternative exon 2 end",
            "Constitutive exon 2 start", geneCol)
        chr <- annotation$Chromosome
        strand <- annotation$Strand
        
        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand,
                                  annotation$`Constitutive exon 1 end`,
                                  annotation$`Alternative exon 1 start`,
                                  showStrand)
        incBstr <- junctionString(chr, strand,
                                  annotation$`Alternative exon 1 end`,
                                  annotation$`Constitutive exon 2 start`,
                                  showStrand)
        excAstr <- junctionString(chr, strand,
                                  annotation$`Constitutive exon 1 end`,
                                  annotation$`Alternative exon 2 start`,
                                  showStrand)
        excBstr <- junctionString(chr, strand,
                                  annotation$`Alternative exon 2 end`,
                                  annotation$`Constitutive exon 2 start`,
                                  showStrand)
        
        # Get specific junction quantification
        coords <- rownames(junctionQuant)
        incA <- junctionQuant[fmatch(incAstr, coords), ]
        incB <- junctionQuant[fmatch(incBstr, coords), ]
        excA <- junctionQuant[fmatch(excAstr, coords), ]
        excB <- junctionQuant[fmatch(excBstr, coords), ]
        
        # Calculate inclusion levels
        inc <- (incA + incB)
        exc <- (excA + excB)
        tot <- inc + exc
        if (nrow(tot) == 0)
            return(NULL)
        
        # Prepare presentation of multigenes
        multigene <- lapply(annotation$Gene, length) > 1
        gene <- annotation$Gene
        gene[multigene] <- lapply(gene[multigene], paste, collapse="/")
        
        eventNames <- paste(sep="_", eventType, chr, strand, 
                            annotation$`Constitutive exon 1 end`,
                            annotation$`Alternative exon 1 start`,
                            annotation$`Alternative exon 1 end`, 
                            annotation$`Alternative exon 2 start`, 
                            annotation$`Alternative exon 2 end`,
                            annotation$`Constitutive exon 2 start`, gene)
    } else if (eventType %in% c("A5SS", "AFE")) {
        # Remove duplicates based on columns used to create identifiers
        annotation <- annotation[!is.na(
            annotation$`Constitutive exon 2 start`), ]
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               "Constitutive exon 1 end",
                               "Alternative exon 1 end",
                               "Constitutive exon 2 start", geneCol)
        chr <- annotation$Chromosome
        strand <- annotation$Strand
        
        # Create searchable strings for junctions
        incStr <- junctionString(chr, strand,
                                 annotation$`Alternative exon 1 end`, 
                                 annotation$`Constitutive exon 2 start`,
                                 showStrand)
        excStr <- junctionString(chr, strand,
                                 annotation$`Constitutive exon 1 end`,
                                 annotation$`Constitutive exon 2 start`,
                                 showStrand)
        
        # Get specific junction quantification
        coords <- rownames(junctionQuant)
        inc <- junctionQuant[fmatch(incStr, coords), ]
        exc <- junctionQuant[fmatch(excStr, coords), ]
        tot <- inc + exc
        if (nrow(tot) == 0)
            return(NULL)
        
        # Prepare presentation of multigenes
        multigene <- lapply(annotation$Gene, length) > 1
        gene <- annotation$Gene
        gene[multigene] <- lapply(gene[multigene], paste, collapse="/")
        
        eventNames <- paste(sep="_", eventType, chr, strand, 
                            annotation$`Constitutive exon 1 end`, 
                            annotation$`Alternative exon 1 end`, 
                            annotation$`Constitutive exon 2 start`, gene)
    } else if (eventType %in% c("A3SS", "ALE")) {
        # Remove duplicates based on columns used to create identifiers
        annotation <- annotation[!is.na(annotation$`Constitutive exon 1 end`), ]
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               "Constitutive exon 1 end",
                               "Alternative exon 1 start",
                               "Constitutive exon 2 start", geneCol)
        chr <- annotation$Chromosome
        strand <- annotation$Strand
        
        # Create searchable strings for junctions
        incStr <- junctionString(chr, strand,
                                 annotation$`Constitutive exon 1 end`,
                                 annotation$`Alternative exon 1 start`,
                                 showStrand)
        excStr <- junctionString(chr, strand,
                                 annotation$`Constitutive exon 1 end`, 
                                 annotation$`Constitutive exon 2 start`,
                                 showStrand)
        
        # Get specific junction quantification
        coords <- rownames(junctionQuant)
        inc <- junctionQuant[fmatch(incStr, coords), ]
        exc <- junctionQuant[fmatch(excStr, coords), ]
        tot <- inc + exc
        if (nrow(tot) == 0)
            return(NULL)
        
        # Prepare presentation of multigenes
        multigene <- lapply(annotation$Gene, length) > 1
        gene <- annotation$Gene
        gene[multigene] <- lapply(gene[multigene], paste, collapse="/")
        
        eventNames <- paste(sep="_", eventType, chr, strand,
                            annotation$`Constitutive exon 1 end`,
                            annotation$`Alternative exon 1 start`, 
                            annotation$`Constitutive exon 2 start`, gene)
    }
    
    # Calculate inclusion levels
    psi <- inc/tot
    
    # Ignore PSI values when total reads are below the threshold
    psi[tot < minReads | is.na(tot)] <- NA
    colnames(psi) <- colnames(inc)
    rownames(psi) <- eventNames
    rm(inc)
    
    # Clear rows with nothing but missing values
    naRows <- rowSums(!is.na(psi)) == 0
    return(psi[!naRows, ])
}