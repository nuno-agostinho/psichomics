#' Get supported splicing event types
#'
#' @param psi Data frame or matrix: alternative splicing quantification data
#' @param acronymsAsNames Boolean: return acronyms as names?
#'
#' @family functions for PSI quantification
#' @return Named character vector with splicing event types
#' @export
#'
#' @examples
#' getSplicingEventTypes()
getSplicingEventTypes <- function(psi=NULL, acronymsAsNames=FALSE) {
    types <- c(
        "Skipped exon"="SE",
        "Mutually exclusive exon"="MXE",
        "Alternative 5' splice site"="A5SS",
        "Alternative 3' splice site"="A3SS",
        "Alternative first exon"="AFE",
        "Alternative last exon"="ALE",
        "Alternative first exon (exon-centred - less reliable)"="AFE_exon",
        "Alternative last exon (exon-centred - less reliable)"="ALE_exon")

    if (!is.null(psi)) {
        types     <- c(types, "Retained intron"="RI")
        dataTypes <- unique(getSplicingEventData(psi)$type)
        if (is.null(dataTypes)) return(NULL)
        # Return both known and unknown event types
        types     <- c(types[types %in% dataTypes],
                       dataTypes[!dataTypes %in% types])
    }
    if (acronymsAsNames) {
        tmp        <- names(types)
        names(tmp) <- types
        types      <- tmp
    }
    return(types)
}

getCoordPosByEventType <- function(type) {
    con1 <- NULL
    alt1 <- NULL
    alt2 <- NULL
    con2 <- NULL
    if (type == "SE") {
        con1 <- 1
        alt1 <- c(2, 3)
        con2 <- 4
    } else if (type == "MXE") {
        con1 <- 1
        alt1 <- c(2, 3)
        alt2 <- c(4, 5)
        con2 <- 6
    } else if (type %in% c("A3SS", "ALE")) {
        con1 <- 1
        alt1 <- 2
        alt2 <- 3
    } else if (type %in% c("A5SS", "AFE")) {
        alt1 <- 2
        alt2 <- 1
        con2 <- 3
    }
    return(list(con1=con1, alt1=alt1, alt2=alt2, con2=con2))
}

#' Check if string identifies splicing events
#'
#' @param char Character vector
#' @param data Object containing event data
#' @param num Integer: number of elements to check
#'
#' @return \code{TRUE} if first elements are splicing events; \code{FALSE},
#'   otherwise
#' @keywords internal
areSplicingEvents <- function(char, data=NULL, num=6) {
    probe <- head(char, num)
    hasUnderscores <- all(sapply(probe, function (i)
        sum(charToRaw(i) == charToRaw("_")) > 3))

    taggedAsPSI <- isTRUE(attr(data, "dataType") == "Inclusion levels")

    hasEventData <- !is.null(findEventData(char, data))
    return(hasUnderscores || taggedAsPSI || hasEventData)
}

#' @importFrom R.utils decapitalize
prettifyEventType <- function(eventType, source=NULL, decap=FALSE) {
    if (is.null(source)) {
        types <- getSplicingEventTypes(acronymsAsNames=TRUE)[eventType]
        types[is.na(types)] <- eventType[is.na(types)]
        names(types)[is.na(names(types))] <- types[is.na(names(types))]
        types <- setNames(sprintf("%s (%s)", types, names(types)), names(types))
    } else if (all(source == "vast-tools")) {
        types <- c("S"="Skipped exon",
                   "C1"="Skipped exon",
                   "C2"="Skipped exon",
                   "C3"="Skipped exon",
                   "ANN"="Skipped exon",
                   "MIC"="Skipped microexon",
                   "IR-C"="Retained intron",
                   "IR-S"="Retained intron",
                   "Alt3"="Alternative 3' splice site",
                   "Alt5"="Alternative 5' splice site")
        types <- setNames(rep(types, 2),
                          c(names(types), paste0("A_", names(types))))
        types <- setNames(sprintf("%s (%s)", types, names(types)), names(types))
        types <- types[eventType]
    }
    if (decap) types <- decapitalize(types)
    return(types)
}

#' Look for event data in input
#'
#' Check if event data can be found in \code{data} and then \code{event}. Event
#' data has to be an object of class \code{eventData}
#'
#' @param event Character: AS event that may contain event data in its
#' attribute \code{eventData}
#' @param data Data frame or matrix: either event data or data containing event
#' data in its attributes \code{rowData} or \code{eventData}
#'
#' @return Event data (or \code{NULL} if not found)
#' @keywords internal
findEventData <- function(event=NULL, data=NULL) {
    returnIfEventData <- function(data) {
        res <- NULL
        if (!is.null(data) && is(data, "eventData")) res <- data
        return(res)
    }

    eventData <- returnIfEventData(attr(data, "rowData"))
    if (is.null(eventData)) {
        eventData <- returnIfEventData(attr(data, "eventData"))
    }
    if (is.null(eventData)) {
        eventData <- returnIfEventData(data)
    }
    if (is.null(eventData)) {
        eventData <- returnIfEventData(attr(event, "eventData"))
    }
    return(eventData)
}

#' Get splicing event information for given alternative splicing quantification
#' data
#'
#' @param psi Matrix or data frame: alternative splicing quantification data
#'
#' @return Matrix or data frame containing splicing event information for
#' alternative splicing events in \code{psi} (if available)
#' @export
getSplicingEventData <- function(psi) findEventData(data=psi)

collapseElemsToSingleString <- function(ll) {
    multi <- sapply(ll, length) > 1
    ll[multi] <- lapply(ll[multi], paste, collapse=" ")
    return(unlist(ll))
}

#' Parse alternative splicing event identifier
#'
#' @param event Character: event identifier
#' @param char Boolean: return character vector instead of list with parsed
#' values?
#' @param pretty Boolean: return a prettier name of the event identifier?
#' @param extra Character: extra information to add (such as species and
#' assembly version); only used if \code{pretty = TRUE} and \code{char = TRUE}
#' @param coords Boolean: display extra coordinates regarding the alternative
#' and constitutive regions of alternative splicing events? Only used if
#' \code{char = FALSE}
#' @param data Matrix or data frame: alternative splicing information
#'
#' @return Data.frame containing type of event, chromosome, strand, gene and
#' position of alternative splicing events or character with that same
#' information (depending on what is available)
#' @export
#'
#' @importFrom shiny isRunning
#'
#' @examples
#' events <- c(
#'   "A3SS_15_+_63353138_63353912_63353397_TPM1",
#'   "A3SS_11_-_61118463_61117115_61117894_CYB561A3",
#'   "A5SS_21_+_48055675_48056459_48056808_PRMT2",
#'   "A5SS_1_-_1274742_1274667_1274033_DVL1",
#'   "AFE_9_+_131902430_131901928_131904724_PPP2R4",
#'   "AFE_5_-_134686513_134688636_134681747_H2AFY",
#'   "ALE_12_+_56554104_56554410_56555171_MYL6",
#'   "ALE_8_-_38314874_38287466_38285953_FGFR1",
#'   "SE_9_+_6486925_6492303_6492401_6493826_UHRF2",
#'   "SE_19_-_5218431_5216778_5216731_5215606_PTPRS",
#'   "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1",
#'   "MXE_17_-_74090495_74087316_74087224_74086478_74086410_74085401_EXOC7")
#' parseSplicingEvent(events)
parseSplicingEvent <- function(event, char=FALSE, pretty=FALSE, extra=NULL,
                               coords=FALSE, data=NULL) {
    eventData <- findEventData(event, data)
    if (is.null(eventData) || !all(event %in% rownames(eventData))) {
        if (all(grepl("(.*_){3,}.*", event))) { # If event has >= 3 underscores
            original <- event
            event <- parseEventFromStr(event=event, char=char, pretty=pretty,
                                       extra=extra, coords=coords)
            if (char) {
                names(event) <- original
            } else {
                rownames(event) <- original
            }
        } else {
            msg <- paste(
                "Cannot parse events. Try running parseSplicingEvent() and",
                "passing alternative splicing quantification to the 'data'",
                "argument.")
            # Show warning if not running Shiny app
            if (!isRunning()) warning(msg)
            if (!char) event <- NULL
        }
        return(event)
    }

    info      <- eventData[event, , drop=FALSE]
    id        <- info$id
    source    <- unique(info$source)
    gene      <- info$gene
    eventType <- info$type
    subtype   <- info$subtype
    if (is.null(subtype)) subtype <- eventType
    chrom     <- info$chrom
    start     <- info$start
    end       <- info$end
    pos       <- info$pos
    strand    <- info$strand

    exonTypes  <- c("constitutive", "alternative")
    exonCoords <- paste0(exonTypes, rep(seq(2), each=2))
    if (char) {
        gene <- prepareGenePresentation(gene)
        if (pretty) {
            sep    <- " "
            extra  <- ifelse(!is.null(extra), paste(",", extra), "")
            parsed <- sprintf(
                "%s %s (chr%s:%s-%s, %s strand%s)",
                gene, prettifyEventType(subtype, source, decap=TRUE),
                chrom, start, end, strand, extra)
        } else {
            sep <- " "
            if (all(exonCoords %in% colnames(info))) {
                isAFEorA5SS <- eventType %in% c("AFE", "A5SS")
                col1        <- ifelse(isAFEorA5SS,
                                      info$alternative2, info$constitutive1)
                col3        <- ifelse(isAFEorA5SS,
                                      info$constitutive1, info$alternative2)
                fullCoords  <- paste(
                    collapseElemsToSingleString(col1),
                    collapseElemsToSingleString(info$alternative1),
                    collapseElemsToSingleString(col3),
                    collapseElemsToSingleString(info$constitutive2))
                fullCoords  <- gsub("NA", "", fullCoords, fixed=TRUE)
                fullCoords  <- trimWhitespace(fullCoords)
                parsed      <- paste(subtype, chrom, strand, fullCoords, gene,
                                     sep=sep)
            } else {
                parsed <- paste(subtype, chrom, strand, gene, sep=sep)
            }
        }
        if (!is.null(id)) parsed <- paste(id, parsed, sep=sep)
        return(parsed)
    } else {
        if (pretty) {
            info$type <- prettifyEventType(eventType, source)
            if (!is.null(info$subtype)) {
                info$subtype <- prettifyEventType(subtype, source)
                if (all(is.na(info$type))) {
                    info$type <- gsub(" \\(.*\\)", "", info$subtype)
                }
            } else {
                info$subtype <- info$type
            }
        }
        if (!coords) info <- info[!colnames(info) %in% exonCoords]
        return(info)
    }
}

parseEventFromStr <- function(event, char=FALSE, pretty=FALSE, extra=NULL,
                              coords=FALSE) {
    if (is.null(event)) return(NULL)
    # Pre-treat special case of exon-centred AFE and ALE
    event <- gsub("AFE_exon", "AFE", event, fixed=TRUE)
    event <- gsub("ALE_exon", "ALE", event, fixed=TRUE)

    # Protect gene symbols made up of underscores
    event <- gsub("(.*)_(Arg|Und|var1|B|[0-9]+)$", "\\1::\\2", event)
    recoverGeneNamesWithUnderscore <- function(gene) {
        lapply(gene, function(g) gsub("::", "_", g, fixed=TRUE))
    }

    if (char) {
        if (pretty) {
            event     <- strsplit(event, "_", fixed=TRUE)
            eventType <- prettifyEventType(sapply(event, "[[", 1), decap=TRUE)
            chrom     <- sapply(event, "[[", 2)
            strand    <- ifelse(sapply(event, "[[", 3) == "+",
                                "positive", "negative")
            gene      <- sapply(event, function(i) i[[length(i)]])
            gene      <- recoverGeneNamesWithUnderscore(gene)

            start     <- sapply(event, "[[", 4)
            end       <- sapply(event, function(i) i[[length(i) - 1]])

            if (is.null(extra)) {
                event <- sprintf("%s %s (chr%s: %s-%s, %s strand)", gene,
                                 eventType, chrom, start, end, strand)
            } else {
                event <- sprintf("%s %s (chr%s: %s-%s, %s strand, %s)", gene,
                                 eventType, chrom, start, end, strand, extra)
            }
        } else {
            event <- gsub("_", " ", event, fixed=TRUE)
        }
        return(event)
    }

    event <- strsplit(event, "_")
    parsed <- data.frame(matrix(nrow=length(event)))

    len <- vapply(event, length, numeric(1))
    lenMinus1 <- len - 1

    parsed$type   <- vapply(event, "[[", 1, FUN.VALUE=character(1))
    parsed$chrom  <- vapply(event, "[[", 2, FUN.VALUE=character(1))
    parsed$strand <- vapply(event, "[[", 3, FUN.VALUE=character(1))
    parsed$gene   <- strsplit(vapply(seq_along(event),
                                     function(i) event[[i]][[len[[i]]]],
                                     FUN.VALUE=character(1)), "/")
    parsed$gene   <- recoverGeneNamesWithUnderscore(parsed$gene)

    # Parse position according to type of alternative splicing event
    parsed$pos <- lapply(seq_along(event),
                         function(i) as.numeric(event[[i]][4:lenMinus1[[i]]]))
    if (coords) {
        parsed$constitutive1 <- NA
        parsed$alternative1  <- NA
        parsed$alternative2  <- NA
        parsed$constitutive2 <- NA
        for (row in seq(nrow(parsed))) {
            type     <- parsed[row, "type"]
            coordPos <- getCoordPosByEventType(type)
            con1     <- coordPos$con1
            alt1     <- coordPos$alt1
            alt2     <- coordPos$alt2
            con2     <- coordPos$con2

            if (!is.null(con1))
                parsed$constitutive1[[row]] <- list(parsed[[row, "pos"]][con1])
            if (!is.null(alt1))
                parsed$alternative1[[row]]  <- list(parsed[[row, "pos"]][alt1])
            if (!is.null(alt2))
                parsed$alternative2[[row]]  <- list(parsed[[row, "pos"]][alt2])
            if (!is.null(con2))
                parsed$constitutive2[[row]] <- list(parsed[[row, "pos"]][con2])
        }

        parsed$alternative1  <- unlist(parsed$alternative1,  recursive=FALSE)
        parsed$constitutive1 <- unlist(parsed$constitutive1, recursive=FALSE)
        parsed$constitutive2 <- unlist(parsed$constitutive2, recursive=FALSE)
        parsed$alternative2  <- unlist(parsed$alternative2,  recursive=FALSE)
    }

    parsed$pos <- suppressWarnings( # Simply ignore non-numeric items
        lapply(parsed$pos, function(i) range(as.numeric(i), na.rm=TRUE)))
    parsed$start <- sapply(parsed$pos, "[[", 1)
    parsed$end   <- sapply(parsed$pos, "[[", 2)

    if (pretty) parsed$type <- prettifyEventType(parsed$type)
    parsed$subtype <- parsed$type

    parsed[,1] <- NULL
    return(parsed)
}

#' Match splicing events with respective genes
#'
#' @param ASevents Character: alternative splicing events to be matched
#' @inheritParams parseSplicingEvent
#'
#' @return Named character vector containing the splicing events and their
#' respective gene as their name
#' @keywords internal
matchSplicingEventsWithGenes <- function(ASevents, data=NULL) {
    ASeventParsed <- parseSplicingEvent(ASevents, data=data)$gene
    if (is.null(ASeventParsed)) return(NULL)

    ASeventGenes  <- rep(ASevents, sapply(ASeventParsed, length))
    names(ASeventGenes) <- unlist(ASeventParsed)
    return(ASeventGenes)
}

#' Get alternative splicing events from genes or vice-versa
#'
#' @details
#' A list of alternative splicing events is required to run
#' \code{getSplicingEventFromGenes}
#'
#' @param genes Character: gene symbols (or TCGA-styled gene symbols)
#' @param ASevents Character: alternative splicing events
#' @inheritParams matchSplicingEventsWithGenes
#'
#' @return Named character containing alternative splicing events or genes and
#' their respective genes or alternative splicing events as names (depending on
#' the function in use)
#'
#' @export
#'
#' @examples
#' ASevents <- c("SE_1_+_201763003_201763300_201763374_201763594_NAV1",
#'               "SE_1_+_183515472_183516238_183516387_183518343_SMG7",
#'               "SE_1_+_183441784_183471388_183471526_183481972_SMG7",
#'               "SE_1_+_181019422_181022709_181022813_181024361_MR1",
#'               "SE_1_+_181695298_181700311_181700367_181701520_CACNA1E")
#' genes <- c("NAV1", "SMG7", "MR1", "HELLO")
#'
#' # Get splicing events from genes
#' matchedASevents <- getSplicingEventFromGenes(genes, ASevents)
#'
#' # Names of matched events are the matching input genes
#' names(matchedASevents)
#' matchedASevents
#'
#' # Get genes from splicing events
#' matchedGenes <- getGenesFromSplicingEvents (ASevents)
#'
#' # Names of matched genes are the matching input alternative splicing events
#' names(matchedGenes)
#' matchedGenes
getSplicingEventFromGenes <- function(genes, ASevents, data=NULL) {
    if (!is(ASevents, "matched")) {
        ASeventGenes <- matchSplicingEventsWithGenes(ASevents, data=data)
    } else {
        ASeventGenes <- ASevents
    }
    genes <- gsub("\\|.*$", "", genes) # Process TCGA-styled gene symbols
    ASeventGenes <- ASeventGenes[names(ASeventGenes) %in% genes]
    return(ASeventGenes)
}

#' @rdname getSplicingEventFromGenes
#' @export
getGenesFromSplicingEvents <- function(ASevents, data=NULL) {
    genes <- parseSplicingEvent(ASevents, data=data)$gene
    match <- unlist(genes)
    names(match) <- rep(ASevents, sapply(genes, length))
    return(match)
}
