#' @include events_mats.R
#' @include events_miso.R
#' @include events_vastTools.R
#' @include events_suppa.R
NULL

#' Creates a template of alternative splicing junctions
#'
#' @param nrow Integer: row number
#' @param program Character: program used to get the junctions
#' @param event.type Character: event type
#' @param chromosome Character: chromosome
#' @param strand Character: positive-sense (\code{+}) or negative-sense
#' (\code{-}) strand
#' @param id Character: event identifiers
#'
#' @return A data frame with the junctions coordinate names pre-filled with
#' \code{NA}
#' @keywords internal
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
#' between different programs?
#'
#' @return Coordinates of interest according to the alternative splicing event
#' type
#' @keywords internal
getSplicingEventCoordinates <- function(type, sorting=FALSE) {
    coords <- switch(type,
                     "SE"   = c("C1.end", "A1.start", "A1.end", "C2.start"),
                     "A3SS" = c("C1.end", "A2.start", "A1.start"),
                     "A5SS" = c("A2.end", "C2.start", "A1.end"),
                     "AFE"  = c("A2.start", "A2.end", "A1.start", "A1.end"),
                     "ALE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
                     "RI"   = c("C1.start", "C1.end", "C2.start", "C2.end"),
                     "MXE"  = c("C1.end", "A1.start", "A1.end",
                                "A2.start", "A2.end", "C2.start"),
                     "TandemUTR" = c("A2.start", "A2.end", "A1.end"))

    if (sorting) {
        coords <- switch(type,
                         "A3SS" = c("A2.start", "A1.start"),
                         "A5SS" = c("A2.end", "A1.end"),
                         "AFE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
                         "ALE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
                         "RI"   = c("C1.end", "C2.start"),
                         "MXE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
                         "TandemUTR" = c("A1.end", "A2.end"))
    }
    return(coords)
}

#' Convert a column to numeric if possible and ignore given columns composed
#' of lists
#'
#' @param table Data matrix: table
#' @param by Character: column names of interest
#' @param toNumeric Boolean: which columns to convert to numeric
#'
#' @return Processed data matrix
#' @keywords internal
#'
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
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return List of events joined by alternative splicing event type
#' @keywords internal
joinEventsPerType <- function(events, types=NULL) {
    if (is.null(types)) types <- names(events)

    pb <- txtProgressBar(max=length(types), style=3)
    joinEvents <- function(k, type, events, pb) {
        type <- types[k]
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
        res <- unique(res)
        setTxtProgressBar(pb, k)
        return(res)
    }
    joint <- lapply(seq(types), joinEvents, types, events, pb)
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
#' @family functions to prepare alternative splicing annotations
#' @return List of data frames with the annotation from different data frames
#' joined by event type
#' @export
#'
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
        events[[each]]$Chromosome <- gsub("^chr", "", events[[each]]$Chromosome)
    }

    # Organise splicing events by event type and then by program in a list of
    # list of dataframes
    events <- rbind.fill(events)
    events <- dlply(events, "Event.type")
    events <- lapply(events, dlply, "Program")

    display("Sorting coordinates...")
    events <- sortCoordinates(events)

    display("Joining events per event type...")
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
    organiseCols <- function(i) {
        data <- join[[i]]
        # Include gene column
        geneCol <- grep("Gene$", names(data))
        if (length(geneCol) == 1) {
            gene <- "Gene"
            names(data)[geneCol] <- gene
            data[[gene]] <- as.list(data[[gene]])
        } else if (length(geneCol) == 0) {
            gene <- NULL
        } else {
            gene <- NULL
            # De-prioritise ENSG gene names
            isENSG <- function(col, data) {
                genes <- data[col]
                # Admit that genes are ENSG if over a given threshold
                thershold <- 0.5
                ensg <- grepl("^ENS(G|T)", unique(genes[!is.na(genes)]))
                return(sum(ensg) >= thershold * length(ensg))
            }
            areENSGcols <- sapply(geneCol, isENSG, data)
            geneCol <- c(geneCol[!areENSGcols], geneCol[areENSGcols])

            vec <- unlist(data[geneCol])
            names(vec) <- rep(seq(nrow(data)), length(geneCol))
            vec <- vec[!is.na(vec)]
            vec <- vec[unique(names(vec))]
            idx <- as.numeric(names(vec))

            gene <- "Gene"
            data[[gene]] <- ""
            data[[gene]][idx] <- vec
        }

        cols <- c("Chromosome", "Strand", gene, getSplicingEventCoordinates(i),
                  grep("Event.ID", names(data), value = TRUE))
        return(data[, cols])
    }
    annot <- lapply(names(join), organiseCols)
    names(annot) <- names(join)
    annot$AFE["C2.start"] <- AFE.C2.start
    annot$ALE["C1.end"]   <- ALE.C1.end
    events <- annot

    display("Cleaning the annotation...")
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
#' @keywords internal
vennEvents <- function(join, eventType) {
    join <- join[[eventType]]

    programs <- join[grep("Program", names(join))]
    nas <- !is.na(programs)
    nas <- ifelse(nas, row(nas), NA)
    p <- lapply(seq(ncol(nas)), function(col) nas[!is.na(nas[ , col]), col])
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
#' @keywords internal
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
#' increasing (plus strand) or decreasing order (minus strand)
#'
#' @param events List of data frames with alternative splicing events for a
#' given program
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return List of data frames with alternative splicing events for a given
#' program
#' @keywords internal
sortCoordinates <- function(events) {
    types <- names(events)

    progress    <- 0
    maxProgress <- sum(sapply(events, length))
    pb <- txtProgressBar(max=maxProgress, style=3)
    for (type in types) {
        coord          <- getSplicingEventCoordinates(type, sorting=TRUE)
        events[[type]] <- colsAsNumbers(type, events)
        programs       <- names(events[[type]])
        if (!is.null(coord)) {
            for (program in programs) {
                table <- events[[type]][[program]]
                plus  <- table[["Strand"]] == "+"
                plusOrd <- apply(table[plus, coord], 1, sort)
                if (length(plusOrd) > 0) {
                    events[[type]][[program]][plus, coord] <- t(plusOrd)
                }
                minusOrd <- apply(table[!plus, coord], 1, sort, decreasing=TRUE)
                if (length(minusOrd) > 0) {
                    events[[type]][[program]][!plus, coord] <- t(minusOrd)
                }
                progress <- progress + 1
                setTxtProgressBar(pb, progress)
            }
        } else {
            progress <- progress + length(programs)
            setTxtProgressBar(pb, progress)
        }
    }
    return(events)
}

#' Prepare presentation of multiple genes for the same splicing event
#'
#' @param gene Character: gene
#' @param collapse Character: character string to separate in case of more than
#'   one gene
#'
#' @return Same object with items collapsed
#' @keywords internal
prepareGenePresentation <- function(gene, collapse="/") {
    multigene <- lapply(gene, length) > 1
    gene[multigene] <- lapply(gene[multigene], paste, collapse=collapse)
    return(gene)
}

prepareEventInfoTooltip <- function(event, data=NULL, pretty=TRUE, ...) {
    parsed <- parseSplicingEvent(event, ..., pretty=pretty, data=data)
    if (is.null(parsed)) return(NULL)
    gene    <- as.character(prepareGenePresentation(parsed$gene,
                                                    collapse=" or "))
    subtype <- parsed$subtype
    strand  <- parsed$strand
    if (!is.null(strand)) {
        strand <- ifelse(strand == "+", "forward", "reverse")
        strand <- sprintf(" (%s strand)", strand)
    } else {
        strand <- ""
    }
    coord <- sprintf("chr %s: %s to %s%s", parsed$chrom,
                     parsed$start, parsed$end, strand)
    res   <- list(gene=gene, subtype=subtype, coord=coord)
    return(res)
}

listPairs <- function(vec1, vec2, vec3=NULL, sorted=FALSE) {
    if (!is.null(vec3)) {
        FUN <- range
    } else if (sorted) {
        FUN <- sort
    } else {
        FUN <- return
    }
    processElems <- function(i, FUN) FUN(c(vec1[[i]], vec2[[i]], vec3[[i]]))
    res <- lapply(seq(length(vec1)), processElems, FUN)
    return(res)
}

#' Calculate inclusion levels using alternative splicing event annotation and
#' junction quantification for many samples
#'
#' @param eventType Character: type of the alternative event to calculate
#' @param junctionQuant Matrix: junction quantification with samples as columns
#' and junctions as rows
#' @param annotation Data.frame: alternative splicing annotation related to
#' event type
#' @param minReads Integer: minimum of total reads required to consider the
#' quantification as valid
#'
#' @importFrom fastmatch fmatch
#'
#' @return Matrix with inclusion levels
#' @keywords internal
calculateInclusionLevels <- function(eventType, junctionQuant, annotation,
                                     minReads = 10,
                                     onlyReturnASeventNames = FALSE) {
    # Immediately return NULL if ALE and AFE events are missing coordinates
    if (eventType %in% c("AFE", "AFE_exon") &&
        is.null(annotation$`Constitutive exon 2 start`)) {
        return(NULL)
    } else if (eventType %in% c("ALE", "ALE_exon") &&
               is.null(annotation$`Constitutive exon 1 end`)) {
        return(NULL)
    }

    if (is.null(annotation$Gene)) {
        geneCol <- NULL
    } else {
        geneCol <- "Gene"
    }
    gene       <- NULL
    coords     <- rownames(junctionQuant)
    showStrand <- any(grepl("\\+|\\-", coords))

    if (eventType == "SE") {
        # Remove duplicates based on columns used to create identifiers
        con1end    <- "Constitutive exon 1 end"
        alt1start  <- "Alternative exon 1 start"
        alt1end    <- "Alternative exon 1 end"
        con2start  <- "Constitutive exon 2 start"
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               con1end, alt1start, alt1end, con2start, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        con1end    <- annotation[[con1end]]
        alt1start  <- annotation[[alt1start]]
        alt1end    <- annotation[[alt1end]]
        con2start  <- annotation[[con2start]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            con1end, alt1start, alt1end, con2start, gene)
        # Prepare event data
        coords_con1 <- I(as.list(con1end))
        coords_alt1 <- I(listPairs(alt1start, alt1end))
        coords_alt2 <- NA
        coords_con2 <- I(as.list(con2start))
        coords_pos  <- I(listPairs(con1end, con2start, sorted=TRUE))

        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand, con1end, alt1start, showStrand)
        incBstr <- junctionString(chr, strand, alt1end, con2start, showStrand)
        excAstr <- junctionString(chr, strand, con1end, con2start, showStrand)

        # Get specific junction quantification
        incA <- fmatch(incAstr, coords)
        incB <- fmatch(incBstr, coords)
        excA <- fmatch(excAstr, coords)
        excB <- 0

        nas  <- is.na(incA) | is.na(incB) | is.na(excA)
        incA <- incA[!nas]
        incB <- incB[!nas]
        excA <- excA[!nas]
        if (length(incA) == 0) return(NULL)
    } else if (eventType == "MXE") {
        # Remove duplicates based on columns used to create identifiers
        con1end    <- "Constitutive exon 1 end"
        alt1start  <- "Alternative exon 1 start"
        alt1end    <- "Alternative exon 1 end"
        alt2start  <- "Alternative exon 2 start"
        alt2end    <- "Alternative exon 2 end"
        con2start  <- "Constitutive exon 2 start"
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               con1end, alt1start, alt1end,
                               alt2start, alt2end, con2start, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        con1end    <- annotation[[con1end]]
        alt1start  <- annotation[[alt1start]]
        alt1end    <- annotation[[alt1end]]
        alt2start  <- annotation[[alt2start]]
        alt2end    <- annotation[[alt2end]]
        con2start  <- annotation[[con2start]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            con1end, alt1start, alt1end,
                            alt2start, alt2end, con2start, gene)
        # Prepare event data
        coords_con1 <- I(as.list(con1end))
        coords_alt1 <- I(listPairs(alt1start, alt1end))
        coords_alt2 <- I(listPairs(alt2start, alt2end))
        coords_con2 <- I(as.list(con2start))
        coords_pos  <- I(listPairs(con1end, con2start, sorted=TRUE))

        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand, con1end, alt1start, showStrand)
        incBstr <- junctionString(chr, strand, alt1end, con2start, showStrand)
        excAstr <- junctionString(chr, strand, con1end, alt2start, showStrand)
        excBstr <- junctionString(chr, strand, alt2end, con2start, showStrand)

        # Get specific junction quantification
        incA <- fmatch(incAstr, coords)
        incB <- fmatch(incBstr, coords)
        excA <- fmatch(excAstr, coords)
        excB <- fmatch(excBstr, coords)

        nas  <- is.na(incA) | is.na(incB) | is.na(excA) | is.na(excB)
        incA <- incA[!nas]
        incB <- incB[!nas]
        excA <- excA[!nas]
        excB <- excB[!nas]
        if (length(incA) == 0) return(NULL)
    } else if (eventType %in% c("A5SS", "AFE")) {
        alt1end   <- "Alternative exon 1 end"
        alt2end   <- "Alternative exon 2 end"
        con2start <- "Constitutive exon 2 start"

        # Backwards compatible with previous annotations
        if (!alt2end %in% names(annotation))
            alt2end <- "Constitutive exon 1 end"

        # Remove duplicates based on columns used to create identifiers
        annotation <- annotation[!is.na(annotation[[con2start]]), ]
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               alt2end, alt1end, con2start, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        alt1end    <- annotation[[alt1end]]
        alt2end    <- annotation[[alt2end]]
        con2start  <- annotation[[con2start]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            alt2end, alt1end, con2start, gene)
        # Prepare event data
        coords_con1 <- NA
        coords_alt1 <- I(as.list(alt1end))
        coords_alt2 <- I(as.list(alt2end))
        coords_con2 <- I(as.list(con2start))
        coords_pos  <- I(listPairs(alt1end, alt2end, con2start, sorted=TRUE))

        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand, alt1end, con2start, showStrand)
        excAstr <- junctionString(chr, strand, alt2end, con2start, showStrand)

        # Get specific junction quantification
        incA <- fmatch(incAstr, coords)
        incB <- 0
        excA <- fmatch(excAstr, coords)
        excB <- 0

        nas  <- is.na(incA) | is.na(excA)
        incA <- incA[!nas]
        excA <- excA[!nas]
        if (length(incA) == 0) return(NULL)
    } else if (eventType %in% c("A3SS", "ALE")) {
        con1end   <- "Constitutive exon 1 end"
        alt1start <- "Alternative exon 1 start"
        alt2start <- "Alternative exon 2 start"

        # Backwards compatible with previous annotations
        if (!alt2start %in% names(annotation))
            alt2start <- "Constitutive exon 2 start"

        # Remove duplicates based on columns used to create identifiers
        annotation <- annotation[!is.na(annotation[[con1end]]), ]
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               con1end, alt1start, alt2start, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        con1end    <- annotation[[con1end]]
        alt1start  <- annotation[[alt1start]]
        alt2start  <- annotation[[alt2start]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            con1end, alt1start, alt2start, gene)
        # Prepare event data
        coords_con1 <- I(as.list(con1end))
        coords_alt1 <- I(as.list(alt1start))
        coords_alt2 <- I(as.list(alt2start))
        coords_con2 <- NA
        coords_pos  <- I(listPairs(con1end, alt1start, alt2start, sorted=TRUE))

        # Create searchable strings for junctions
        incAstr <- junctionString(chr, strand, con1end, alt1start, showStrand)
        excAstr <- junctionString(chr, strand, con1end, alt2start, showStrand)

        # Get specific junction quantification
        incA <- fmatch(incAstr, coords)
        incB <- 0
        excA <- fmatch(excAstr, coords)
        excB <- 0

        nas  <- is.na(incA) | is.na(excA)
        incA <- incA[!nas]
        excA <- excA[!nas]
        if (length(incA) == 0) return(NULL)
    } else if (eventType == "AFE_exon") {
        alt1end      <- "Alternative exon 1 end"
        alt2end      <- "Alternative exon 2 end"

        # Backwards compatible with previous annotations
        if (!alt2end %in% names(annotation))
            alt2end <- "Constitutive exon 1 end"

        # Remove duplicates based on columns used to create identifiers
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               alt2end, alt1end, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        alt1end <- annotation[[alt1end]]
        alt2end <- annotation[[alt2end]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            alt2end, alt1end, gene)
        # Prepare event data
        coords_con1 <- NA
        coords_alt1 <- I(as.list(alt1end))
        coords_alt2 <- I(as.list(alt2end))
        coords_con2 <- NA
        coords_pos  <- I(listPairs(alt1end, alt2end, sorted=TRUE))

        # Create searchable strings for junctions
        incStr <- junctionString(chr, strand, alt1end, "", showStrand)
        excStr <- junctionString(chr, strand, alt2end, "", showStrand)
        coords <- gsub(":[0-9]*:([\\+\\-])$", "::\\1", coords)

        # Get specific junction quantification
        inc <- lapply(incStr, grep, coords, fixed=TRUE)
        exc <- lapply(excStr, grep, coords, fixed=TRUE)

        nas <- is.na(inc) | is.na(exc)
        inc <- inc[!nas]
        exc <- exc[!nas]
        if (length(inc) == 0) return(NULL)
    } else if (eventType == "ALE_exon") {
        alt1start    <- "Alternative exon 1 start"
        alt2start    <- "Alternative exon 2 start"

        # Backwards compatible with previous annotations
        if (!alt2start %in% names(annotation))
            alt2start <- "Constitutive exon 2 start"

        # Remove duplicates based on columns used to create identifiers
        annotation <- uniqueBy(annotation, "Chromosome", "Strand",
                               alt1start, alt2start, geneCol)
        chr    <- annotation$Chromosome
        strand <- annotation$Strand
        if (!is.null(geneCol)) {
            gene <- prepareGenePresentation(annotation[[geneCol]])
        }

        alt1start <- annotation[[alt1start]]
        alt2start <- annotation[[alt2start]]
        eventNames <- paste(sep="_", eventType, chr, strand,
                            alt1start, alt2start, gene)
        # Prepare event data
        coords_con1 <- NA
        coords_alt1 <- I(as.list(alt1start))
        coords_alt2 <- I(as.list(alt2start))
        coords_con2 <- NA
        coords_pos  <- I(listPairs(alt1start, alt2start, sorted=TRUE))

        # Create searchable strings for junctions
        incStr <- junctionString(chr, strand, "", alt1start, showStrand)
        excStr <- junctionString(chr, strand, "", alt2start, showStrand)
        coords <- gsub(":.*?:", "::", coords)

        # Get specific junction quantification
        inc <- lapply(incStr, grep, coords, fixed=TRUE)
        exc <- lapply(excStr, grep, coords, fixed=TRUE)

        nas <- is.na(inc) | is.na(exc)
        inc <- inc[!nas]
        exc <- exc[!nas]
        if (length(inc) == 0) return(NULL)
    }

    # Warn about events with incomplete information to be calculated
    valid <- sum(!nas)
    total <- length(nas)
    perc  <- round(valid / total * 100)
    message(sprintf(
        paste("Using %s of %s events (%s%%) whose junctions are present in",
              "junction quantification data..."),
        valid, total, perc))

    if (eventType %in% c("ALE_exon", "AFE_exon")) {
        psi <- psiFastCalc2(junctionQuant, inc=inc, exc=exc, minReads=minReads)
    } else {
        psi <- psiFastCalc(junctionQuant, incA=incA, incB=incB, excA=excA,
                           excB=excB, minReads=minReads)
    }

    if (nrow(psi) == 0) return(NULL)
    rownames(psi) <- eventNames[!nas]

    # Clear events with nothing but missing values
    validEvents <- rowSums(!is.na(psi)) > 0
    invalid     <- sum(!validEvents)
    total       <- length(validEvents)
    perc        <- round(invalid / total * 100)
    if (invalid > 0) {
        message(sprintf(
            paste("Discarding %s of %s events (%s%%) containing less than %s",
                  "read counts across all samples..."),
            invalid, total, perc, minReads))
        psi <- psi[validEvents, , drop=FALSE]
    }

    # Finalise AS event information
    gene <- NULL
    if (!is.null(geneCol)) gene <- annotation[[geneCol]]
    if (!is.null(gene)) {
        gene <- I(gene)
    } else {
        gene <- NA
    }
    eventData <- data.frame(type=eventType, subtype=eventType,
                            chrom=chr, strand=strand, gene=gene,
                            start=sapply(coords_pos, min),
                            end=sapply(coords_pos, max),
                            pos=coords_pos,
                            constitutive1=coords_con1,
                            alternative1=coords_alt1,
                            alternative2=coords_alt2,
                            constitutive2=coords_con2)
    rownames(eventData)    <- eventNames
    attr(psi, "eventData") <- eventData[rownames(psi), , drop=FALSE]
    return(psi)
}
