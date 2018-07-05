## Auxiliary functions used throughout the program

# Print how to start the graphical interface when attaching the package
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Start the visual interface by running the function ",
                          "psichomics()")
}

#' Round down/up the minimum/maximum value
#' @param x Numeric: values
#' @param digits Numeric: number of maximum digits
#' 
#' @return Rounded numeric value
roundMinDown <- function(x, digits=0) floor  (min(x) * 10^digits) / 10^digits

#' @rdname roundMinDown
roundMaxUp   <- function(x, digits=0) ceiling(max(x) * 10^digits) / 10^digits

#' Get psichomics file inside a given directory
#' @inheritParams base::system.file
#' @return Loaded file
insideFile <- function(...) {
    return(system.file(..., package="psichomics"))
}

#' Sidebar without a well
#' 
#' Modified version of \code{shiny::sidebarPanel} without a well
#' 
#' @importFrom shiny div tags
#' 
#' @inherit shiny::sidebarPanel
sidebar <- function(..., width=4) {
    div(class = paste0("col-sm-", width), tags$form(...))
}

#' Load local file
#' @param file Character: path to the file
#' 
#' @return Loaded file
#' @export
#' 
#' @examples 
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
readFile <- function(file) {
    readRDS(insideFile("extdata", file))
}

#' Link to run arbitrary JavaScript code
#' 
#' @param text Character: text label
#' @param code Character: JavaScript code
#' 
#' @return HTML elements
linkToRunJS <- function(text, code) {
    HTML(sprintf('<a href="#" onclick="%s; return false;">%s</a>', code, text))
}

#' Create a row for a HTML table
#' 
#' @param ... Elements to include in the row
#' @param th Boolean: is this row the table head?
#' 
#' @return HTML elements
tableRow <- function (..., th=FALSE) {
    args <- list(...)
    if (th) row <- tags$th
    else    row <- tags$td
    do.call(tags$tr, lapply(args, row))
}

#' Splicing event types available
#' 
#' @param acronymsAsNames Boolean: return acronyms as names? FALSE by default
#' 
#' @return Named character vector with splicing event types
#' @export
#' 
#' @examples 
#' getSplicingEventTypes()
getSplicingEventTypes <- function(acronymsAsNames=FALSE) {
    types <- c(
        "Skipped exon"="SE",
        "Mutually exclusive exon"="MXE",
        "Alternative 5' splice site"="A5SS",
        "Alternative 3' splice site"="A3SS",
        "Alternative first exon"="AFE",
        "Alternative last exon"="ALE",
        "Alternative first exon (exon-centred - less reliable)"="AFE_exon",
        "Alternative last exon (exon-centred - less reliable)"="ALE_exon")
    if (acronymsAsNames) {
        tmp        <- names(types)
        names(tmp) <- types
        types      <- tmp
    }
    return(types)
}

#' Check if string identifies splicing events
#' 
#' @param char Character vector
#' @param num Integer: number of elements to check
#' 
#' @return TRUE if first elements of the vector identify splicing events; FALSE,
#'  otherwise
areSplicingEvents <- function(char, num=6) {
    all(sapply(head(char, num), function (i)
        sum(charToRaw(i) == charToRaw("_")) > 3))
}

#' Parse an alternative splicing event based on a given identifier
#' 
#' @param event Character: event identifier
#' @param char Boolean: return a single character instead of list with parsed
#' values?
#' @param pretty Boolean: return a prettier name of the event identifier?
#' @param extra Character: extra information to add (such as species and
#' assembly version); only used if \code{pretty} and \code{char} are \code{TRUE}
#' @param coords Boolean: extra coordinates regarding the alternative and 
#' constitutive regions of alternative splicing events; only used if \code{char}
#' is \code{FALSE}
#' 
#' @return Parsed event
#' @export
#' @examples 
#' events <- c("SE_1_-_123_456_789_1024_TST",
#'             "MXE_3_+_473_578_686_736_834_937_HEY/YOU")
#' parseSplicingEvent(events)
parseSplicingEvent <- function(event, char=FALSE, pretty=FALSE, extra=NULL, 
                               coords=FALSE) {
    # Pre-treat special case of exon-centred AFE and ALE
    event <- gsub("AFE_exon", "AFE", event, fixed=TRUE)
    event <- gsub("ALE_exon", "ALE", event, fixed=TRUE)
    
    if (char) {
        if (pretty) {
            event     <- strsplit(event, "_", fixed=TRUE)
            
            eventType <- sapply(event, "[[", 1)
            eventType <- getSplicingEventTypes(acronymsAsNames=TRUE)[eventType]
            eventType <- tolower(eventType)
            
            chrom     <- sapply(event, "[[", 2)
            strand    <- ifelse(sapply(event, "[[", 3) == "+", "positive",
                                "negative")
            gene      <- sapply(event, function(i) i[[length(i)]])
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
    
    # Parse position according to type of alternative splicing event
    parsed$pos <- lapply(seq_along(event), 
                         function(i) event[[i]][4:lenMinus1[[i]]])
    if (coords) {
        parsed$constitutive1 <- NA
        parsed$alternative1  <- NA
        parsed$alternative2  <- NA
        parsed$constitutive2 <- NA
        for (row in seq(nrow(parsed))) {
            type <- parsed[row, "type"]
            
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
    
    if (pretty)
        parsed$type <- getSplicingEventTypes(acronymsAsNames=TRUE)[parsed$type]
    
    parsed[,1] <- NULL
    return(parsed)
}

parseEvent <- parseSplicingEvent

#' Match splicing events with respective genes
#' 
#' @param ASevents Character: alternative splicing events to be matched
#' 
#' @return Named character vector containing the splicing events and their
#' respective gene as their name
matchSplicingEventsWithGenes <- function(ASevents) {
    ASeventParsed <- parseSplicingEvent(ASevents)$gene
    ASeventGenes  <- rep(ASevents, sapply(ASeventParsed, length))
    names(ASeventGenes) <- unlist(ASeventParsed)
    return(ASeventGenes)
}

#' Retrieve alternative splicing events based on given genes
#' 
#' @param genes Character: gene symbols (or TCGA-styled gene symbols)
#' @inheritParams matchSplicingEventsWithGenes
#' 
#' @return Character containing respective alternative splicing events
#' @export
#' 
#' @examples 
#' ASevents <- c("SE_1_+_201763003_201763300_201763374_201763594_NAV1", 
#'               "SE_1_+_183515472_183516238_183516387_183518343_SMG7",
#'               "SE_1_+_183441784_183471388_183471526_183481972_SMG7",
#'               "SE_1_+_181019422_181022709_181022813_181024361_MR1",
#'               "SE_1_+_181695298_181700311_181700367_181701520_CACNA1E")
#' genes <- c("NAV1", "SMG7", "MR1", "HELLO")
#' getSplicingEventFromGenes(genes, ASevents)
getSplicingEventFromGenes <- function(genes, ASevents) {
    if (!is(ASevents, "matched"))
        ASeventGenes <- matchSplicingEventsWithGenes(ASevents)
    else
        ASeventGenes <- ASevents
    
    genes <- gsub("\\|.*$", "", genes) # Process TCGA-styled gene symbols
    ASeventGenes <- ASeventGenes[names(ASeventGenes) %in% genes]
    return(ASeventGenes)
}

#' Retrieve genes based on given alternative splicing events
#' 
#' @param ASevents Character: alternative splicing events
#' 
#' @return Named character containing alternative splicing events and their
#' respective genes as names
#' @export
getGenesFromSplicingEvents <- function(ASevents) {
    genes <- parseSplicingEvent(ASevents)$gene
    match <- unlist(genes)
    names(match) <- rep(ASevents, sapply(genes, length))
    return(match)
}

#' Trims whitespace from a word
#'
#' @param word Character to trim
#'
#' @return Character without whitespace
#'
#' @examples
#' psichomics:::trimWhitespace("    hey   there     ")
#' psichomics:::trimWhitespace(c("pineapple    ", "one two three", 
#'                               " sunken    ship   "))
trimWhitespace <- function(word) {
    # Remove leading and trailing whitespace
    word <- gsub("^\\s+|\\s+$", "", word)
    # Replace multiple spaces between words with one single space
    word <- gsub("\\s+", " ", word)
    return(word)
}

#' Filter NULL elements from vector or list
#' 
#' @param v Vector or list
#' 
#' @return Filtered vector or list with no NULL elements; if the input is a
#' vector composed of only NULL elements, it returns a NULL (note that it will
#' returns an empty list if the input is a list with only NULL elements)
rm.null <- function(v) Filter(Negate(is.null), v)

#' Escape symbols for use in regular expressions
#'
#' @param ... Characters to be pasted with no space
#' 
#' @return Escaped string
escape <- function(...) {
    # return(gsub("/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]", "\\$&", string))
    return(gsub("(\\W)", "\\\\\\1", paste0(...)))
}

#' Convert vector of values to JavaScript array
#' 
#' @param values Character vector
#' 
#' @return Character with valid JavaScript array
toJSarray <- function(values) {
    paste0("[", paste0(paste0("\'", values, "\'"), collapse=", "), "]")
}

#' Check if a number is whole
#' 
#' @param x Object to be tested
#' @param tol Numeric: tolerance used for comparison
#' 
#' @return TRUE if number is whole; otherwise, FALSE
is.whole <- function(x, tol=.Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}

#' Calculate mean for each row of a matrix
#' 
#' @param mat Matrix
#' @param na.rm Boolean: remove NAs?
#' 
#' @return Vector of means
#' @export
#' 
#' @examples
#' df <- rbind("Gene 1"=c(3, 5, 7), "Gene 2"=c(8, 2, 4), "Gene 3"=c(9:11)) 
#' rowMeans(df)
rowMeans <- function(mat, na.rm=FALSE) {
    if ( !is.null(dim(mat)) ) {
        nas <- 0
        if (na.rm) nas <- rowSums(is.na(mat))
        rowSums(mat, na.rm=na.rm) / (ncol(mat) - nas)
    } else {
        mean(mat, na.rm=na.rm)
    }
}

#' Calculate variance for each row of a matrix
#' 
#' @inheritParams rowMeans
#' 
#' @return Vector of variances
#' @export
#' 
#' @examples 
#' df <- rbind("Gene 1"=c(3, 5, 7), "Gene 2"=c(8, 2, 4), "Gene 3"=c(9:11)) 
#' rowVars(df)
rowVars <- function(mat, na.rm=FALSE) {
    if ( !is.null(dim(mat)) ) {
        means      <- rowMeans(mat, na.rm=na.rm)
        meansSqDev <- (mat - means) ** 2
        squaresSum <- rowSums(meansSqDev, na.rm=na.rm)
        
        nas <- 0
        if (na.rm) nas <- rowSums(is.na(mat))
        dem <- ncol(mat) - nas - 1
        squaresSum/dem
    } else {
        var(mat, na.rm=na.rm)
    }
}

#' Rename vector to avoid duplicated values with another vector
#'
#' Renames values by adding an index to the end of duplicates. This allows to
#' prepare unique values in two vectors before a merge, for instance.
#'
#' @param check Character: values to rename if duplicated
#' @param comp Character: values to compare with
#'
#' @importFrom utils head
#'
#' @return Character vector with renamed values if duplicated; else, it
#' returns the usual values. It does not return the comparator values.
#'
#' @examples
#' psichomics:::renameDuplicated(check = c("blue", "red"), comp = c("green",
#'                                                                  "blue"))
renameDuplicated <- function(check, comp) {
    # If there's nothing to compare with, return the values
    # if (length(comp) == 0) return(check)
    
    repeated <- check %in% comp | duplicated(check)
    uniq <- check[!repeated]
    
    for (dup in which(repeated)) {
        # Locate matches (don't forget the counter)
        all <- c(comp, head(check[1:dup], -1))
        expr <- paste0(escape(check[dup]), " \\([0-9]+\\)|", escape(check[dup]))
        locate <- grep(expr, all, value = TRUE)
        
        # Get the maximum counter and add one
        counter <- sub(".* \\(([0-9]+)\\)", "\\1", locate)
        
        # Replace strings with 0
        counter[grep("^[0-9]*$", counter, invert =TRUE)] <- 0
        check[dup] <- sprintf("%s (%i)", check[dup], 
                              max(as.numeric(counter)) + 1)
    }
    return(check)
}

#' Style button used to initiate a process
#' 
#' @param id Character: button identifier
#' @param label Character: label
#' @inheritDotParams shiny::actionButton -inputId -label
#' @param class Character: class
#' 
#' @importFrom shinyjs hidden
#' @importFrom shiny tags actionButton
#' @return HTML for a button
processButton <- function(id, label, ..., class="btn-primary") {
    spinner <- tags$i(id=paste0(id, "Loading"), class="fa fa-spinner fa-spin")
    button  <- actionButton(id, class=class, type="button", 
                            label=div(hidden(spinner), label), ...)
    return(button)
}

#' Signal the program that a process is starting
#' 
#' Style button to show processing is in progress
#' 
#' @param id Character: button identifier
#' @importFrom shinyjs show
#' @return Start time of the process
startProcess <- function(id) {
    disable(id)
    show(paste0(id, "Loading"))
    return(Sys.time())
}

#' Signal the program that a process has ended
#' 
#' Style button to show processing is not occurring. Also, close the progress
#' bar (if TRUE) and print the difference between the current time and a given
#' time (if given time is not NULL)
#' 
#' @param id Character: button identifier
#' @param time \code{POSIXct} object: start time needed to show the interval
#' time (if NULL, the time interval is not displayed)
#' @param closeProgressBar Boolean: close progress bar? TRUE by default
#' 
#' @importFrom shinyjs hide
#' @return NULL (this function is used to modify the Shiny session's state)
endProcess <- function(id, time=NULL, closeProgressBar=TRUE) {
    enable(id)
    hide(paste0(id, "Loading"))
    if (closeProgressBar) suppressWarnings(closeProgress())
    if (!is.null(time)) display(Sys.time() - time, "Process finished in")
}

#' Get patients from given samples
#'
#' @param sampleId Character: sample identifiers
#' @param patientId Character: patient identifiers to filter by (optional; if a
#' matrix or data frame is given, its rownames will be used to infer the patient
#' identifiers)
#' @param na Boolean: return NA for samples with no matching patients
#' @param sampleInfo Data frame or matrix: sample information containing the
#' sample identifiers as rownames and a column named "Subject ID" with the
#' respective subject identifiers
#'
#' @return Character: patient identifiers corresponding to the given samples
#' 
#' @export
#' @examples
#' samples <- paste0("GTEX-", c("ABC", "DEF", "GHI", "JKL", "MNO"), "-sample")
#' getPatientFromSample(samples)
#' 
#' # Filter returned samples based on available patients
#' patients <- paste0("GTEX-", c("DEF", "MNO"))
#' getPatientFromSample(samples, patients)
getPatientFromSample <- function(sampleId, patientId=NULL, na=FALSE,
                                 sampleInfo=NULL) {
    if (!is.null(patientId) && 
        (is.matrix(patientId) || is.data.frame(patientId))) {
        patientId <- rownames(patientId)
    }
    
    # Extract patient identifiers from sample ID and then retrieve their index
    extractPatientIndex <- function(pattern, samples, allPatients) {
        patient <- gsub(pattern, "\\1", samples)
        names(patient) <- samples
        
        if (!is.null(allPatients)) {
            # Filter by patients of interest
            if (na) {
                # Return NA as the corresponding patient
                patient[!patient %in% allPatients] <- NA
            } else {
                patient <- patient[patient %in% allPatients]
            }
        }
        return(patient)
    }
    
    if ( any(grepl("^TCGA", sampleId)) ) {
        # Retrieve TCGA patient index
        extractPatientIndex("(TCGA-.*?-.*?)-.*", sampleId, patientId)
    } else if ( any(grepl("^GTEX", sampleId)) ) {
        # Retrieve GTEx patient index
        extractPatientIndex("(GTEX-.*?)-.*", sampleId, patientId)
    } else if ( "Subject ID" %in% colnames(sampleInfo) ) {
        # Based on user-owned files
        patients <- as.character(sampleInfo[ , "Subject ID"])
        names(patients) <- rownames(sampleInfo)
        return(patients)
    } else {
        return(NULL)
    }
}

#' @rdname getPatientFromSample
#' @export
getSubjectFromSample <- getPatientFromSample

#' Get samples matching the given patients
#' 
#' @inheritParams getPatientFromSample
#' @param patients Character or list of characters: patient identifiers
#' @param samples Character: sample identifiers
#' @param clinical Data frame or matrix: clinical dataset
#' @param rm.NA Boolean: remove NAs? TRUE by default
#' @param match Integer: vector of patient index with the sample identifiers as
#' name to save time (optional)
#' @param showMatch Boolean: show matching patient index? FALSE by default
#' 
#' @return Names of the matching samples (if \code{showMatch} is \code{TRUE},
#' a character with the patients as values and their respective samples as names
#' is returned)
#' @export
#' 
#' @examples 
#' patients <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
#' samples <- paste0(patients, "-sample")
#' clinical <- data.frame(samples=samples)
#' rownames(clinical) <- patients
#' getMatchingSamples(patients[c(1, 4)], samples, clinical)
getMatchingSamples <- function(patients, samples, clinical=NULL, rm.NA=TRUE,
                               match=NULL, showMatch=FALSE) {
    if (is.null(match))
        match <- getPatientFromSample(samples, clinical)
    
    if (is.list(patients)) {
        samples <- lapply(patients, function(i) {
            res <- match[match %in% i]
            if (!showMatch) res <- unique(names(res))
            if (rm.NA) res <- res[!is.na(res)]
            return(res)
        })
    } else {
        samples <- match[match %in% patients]
        if (!showMatch) samples <- unique(names(samples))
        if (rm.NA) samples <- samples[!is.na(samples)]
    }
    return(samples)
}

#' @rdname getMatchingSamples
#' @export
getSampleFromPatient <- getMatchingSamples

#' @rdname getMatchingSamples
#' @export
getSampleFromSubject <- getMatchingSamples

#' Assign one group to each element
#' 
#' @param groups List of integers: groups of elements
#' @param elem Character: all elements available (NULL by default)
#' @param outerGroupName Character: name to give to outer group (NA by default; 
#' set to NULL to only show elements matched to their respective groups)
#' 
#' @return Character vector where each element corresponds to the group of the
#' respective element
#' @export
#' 
#' @examples 
#' groups <- list(1:3, 4:7, 8:10)
#' names(groups) <- paste("Stage", 1:3)
#' groupPerElem(groups)
groupPerElem <- function(groups, elem=NULL, outerGroupName=NA) {
    if (length(groups) == 0) {
        singleGroup <- "Single group"
        if (length(elem) > 0)
            return(rep(singleGroup, length(elem)))
        else
            return(singleGroup)
    }
    
    all <- unlist(groups)
    names(all) <- rep(names(groups), sapply(groups, length))
    
    finalGroups <- NULL
    if (!is.null(elem) && !is.null(outerGroupName)) {
        finalGroups <- rep(outerGroupName, length(elem))
        names(finalGroups) <- elem
    }
    
    colour <- NULL
    assignedColours <- attr(groups, "Colour")
    for (each in unique(all)) {
        each <- as.character(each) # Force to use numeric identifiers as names
        groupName <- names(all[all == each])
        groupNameStr <- paste(groupName, collapse=", ")
        if (!is.null(assignedColours) && !groupNameStr %in% names(colour)) {
            cols   <- assignedColours[groupName]
            colour <- c(colour, setNames(Reduce(blendColours, cols),
                                         groupNameStr))
        }
        finalGroups[each] <- groupNameStr
    }
    attr(finalGroups, "Colour") <- colour
    return(finalGroups)
}

#' @rdname groupPerElem
#' 
#' @param patients Integer: total number of patients
#' @param includeOuterGroup Boolean: join the patients that have no groups?
#' 
#' @export
groupPerPatient <- function(groups, patients=NULL, includeOuterGroup=FALSE, 
                            outerGroupName="(Outer data)") {
    .Deprecated("groupPerElem")
    if (!includeOuterGroup) outerGroupName <- NULL
    groupPerElem(groups, patients, outerGroupName)
}

#' @rdname groupPerElem
#' @param samples Character: all available samples
#' @export
groupPerSample <- function(groups, samples, includeOuterGroup=FALSE, 
                           outerGroupName="(Outer data)") {
    .Deprecated("groupPerElem")
    if (!includeOuterGroup) outerGroupName <- NULL
    groupPerElem(groups, samples, outerGroupName)
}

#' @rdname missingDataModal
#' @param modal Character: modal identifier
loadRequiredData <- function( modal=NULL ) {
    modal <- ifelse(is.null(modal), "null", modal)
    return(sprintf("showDataPanel('#%s');", modal))
}

#' Style and show a modal
#' 
#' You can also use \code{errorModal} and \code{warningModal} to use a template 
#' modal already stylised to show errors and warnings, respectively.
#' 
#' @param session Current Shiny session
#' @param title Character: modal title
#' @param ... Extra arguments to pass to \code{shiny::modalDialog}
#' @param style Character: style of the modal (NULL, "warning", "error" or 
#' "info"; NULL by default)
#' @param iconName Character: FontAwesome icon name to appear with the title
#' @param footer HTML elements to use in footer
#' @param echo Boolean: print to console? FALSE by default
#' @param size Character: size of the modal - "medium" (default), "small" or 
#' "large"
#' @param dismissButton Boolean: show dismiss button in footer? TRUE by default
#' 
#' @importFrom shiny renderUI div icon showModal modalButton modalDialog
#' @importFrom shinyBS toggleModal
#' @seealso \code{\link{showAlert}}
#' @return NULL (this function is used to modify the Shiny session's state)
styleModal <- function(session, title, ..., style=NULL,
                       iconName="exclamation-circle", footer=NULL, echo=FALSE, 
                       size="medium", dismissButton=TRUE) {
    
    size <- switch(size, "small"="s", "large"="l", "medium"="m")
    if (dismissButton) footer <- tagList(modalButton("Dismiss"), footer)
    
    modal <- modalDialog(..., title=div(icon(iconName), title), size=size,
                         footer=footer, easyClose=FALSE)
    if (!is.null(style)) {
        style <- match.arg(style, c("info", "warning", "error"))
        modal[[3]][[1]][[3]][[1]][[3]][[1]] <-
            tagAppendAttributes(modal[[3]][[1]][[3]][[1]][[3]][[1]],
                                class=style)
    }
    showModal(modal, session)
    if (echo) display(content)
}

#' @rdname styleModal
errorModal <- function(session, title, ..., size="small", footer=NULL) {
    styleModal(session, title, ..., footer=footer, style="error", size=size,
               echo=FALSE, iconName="times-circle")
}

#' @rdname styleModal
warningModal <- function(session, title, ..., size="small", footer=NULL) {
    styleModal(session, title, ..., footer=footer, style="warning", size=size,
               echo=FALSE, iconName="exclamation-circle")
}

#' @rdname styleModal
infoModal <- function(session, title, ..., size="small", footer=NULL) {
    styleModal(session, title, ..., footer=footer, style="info", size=size,
               echo=FALSE, iconName="info-circle")
}

#' Show or remove an alert
#' 
#' You can also use \code{errorAlert} and \code{warningAlert} to use template 
#' alerts already stylised to show errors and warnings respectively.
#' 
#' @param session Shiny session
#' @param ... Arguments to render as elements of alert
#' @param title Character: title of the alert (optional)
#' @param style Character: style of the alert ("alert-danger", "alert-warning" 
#' or NULL)
#' @param dismissible Boolean: is the alert dismissible? TRUE by default
#' @param alertId Character: alert identifier
#' 
#' @seealso \code{\link{showModal}}
#' @importFrom shiny span h3 renderUI div tagList
#' @return NULL (this function is used to modify the Shiny session's state)
showAlert <- function(session, ..., title=NULL, style=NULL, dismissible=TRUE, 
                      alertId="alert") {
    ns <- session$ns
    
    if (dismissible) {
        dismissible <- "alert-dismissible"
        dismiss <- tags$button(type="button", class="close",
                               "data-dismiss"="alert", "aria-label"="Close",
                               span("aria-hidden"="true", "\u00D7"))
    } else {
        dismissible <- NULL
        dismiss <- NULL
    }
    
    if (!is.null(title)) title <- h4(title)
    
    output <- session$output
    output[[alertId]] <- renderUI({
        tagList(
            div(title, id="myAlert", class="alert", class=style, role="alert",
                class="animated bounceInUp", class=dismissible, dismiss, ...)
        )
    })
}

#' @rdname showAlert
errorAlert <- function(session, ..., title=NULL, dismissible=TRUE,
                       alertId="alert") {
    showAlert(session, ..., style="alert-danger", title=title, 
              dismissible=dismissible, alertId=alertId)
}

#' @rdname showAlert
warningAlert <- function(session, ..., title=NULL, dismissible=TRUE,
                         alertId="alert") {
    showAlert(session, ..., style="alert-warning", title=title,
              dismissible=dismissible, alertId=alertId)
}

#' @rdname showAlert
#' 
#' @param output Shiny output
removeAlert <- function(output, alertId="alert") {
    output[[alertId]] <- renderUI(NULL)
}

#' Alert in the style of a dialogue box with a button
#' 
#' @param id Character: identifier (NULL by default)
#' @param description Character: description
#' @param buttonId Character: button identifier (NULL by default)
#' @param buttonLabel Character: button label (NULL by default)
#' @param buttonIcon Character: button icon (NULL by default)
#' @param ... Extra parameters when creating the alert
#' @param type Character: type of alert (error or warning)
#' @param bigger Boolean: wrap the \code{description} in a \code{h4} tag?
#'
#' @importFrom shiny icon div actionButton
#'
#' @return HTML elements
inlineDialog <- function(description, ..., buttonLabel=NULL, buttonIcon=NULL, 
                         buttonId=NULL, id=NULL, type=c("error", "warning"),
                         bigger=FALSE) {
    type <- match.arg(type)
    if (identical(type, "error")) type <- "danger"
    typeIcon <- switch(type, danger="exclamation-circle",
                       warning="exclamation-triangle")
    
    if (!is.null(buttonLabel)) {
        if (!is.null(buttonIcon))
            icon <- icon(buttonIcon)
        else
            icon <- NULL
        
        typeClass <- sprintf("btn-%s btn-block", type)
        button <- tagList(br(), br(), actionButton(buttonId, icon=icon, 
                                                   buttonLabel,
                                                   class=typeClass))
    } else {
        button <- NULL
    }
    
    if (bigger) {
        description <- h4(style="margin-top: 5px !important;",
                          icon(typeIcon), description)
    } else {
        description <- tagList(icon(typeIcon), description)
    }
    
    typeClass <- sprintf("alert alert-%s", type)
    div(id=id, class=typeClass, role="alert", style="margin-bottom: 0px;",
        description, button, ...)
}

#' @rdname inlineDialog
errorDialog <- function(description, ...)
    inlineDialog(description, ..., type="error")

#' @rdname inlineDialog
warningDialog <- function(description, ...)
    inlineDialog(description, ..., type="warning")

#' Get the Downloads folder of the user
#' @return Path to Downloads folder
#' @export
#' 
#' @examples 
#' getDownloadsFolder()
getDownloadsFolder <- function() {
    if (Sys.info()['sysname'] == "Windows")
        folder <- dirname("~")
    else
        folder <- path.expand("~")
    folder <- file.path(folder, "Downloads")
    folder <- paste0(folder, .Platform$file.sep)
    return(folder)
}

#' Return the type of a given sample
#' 
#' @param sample Character: ID of the sample
#' @param filename Character: path to RDS file containing corresponding type
#' 
#' @return Types of the TCGA samples
#' @export
#' @examples 
#' parseSampleGroups(c("TCGA-01A-Tumour", "TCGA-10B-Normal"))
parseSampleGroups <- function(sample, 
                              filename = system.file("extdata",  
                                                     "TCGAsampleType.RDS",
                                                     package="psichomics")) {
    typeList <- readRDS(filename)
    type <- gsub(".*?-([0-9]{2}).-.*", "\\1", sample, perl = TRUE)
    return(typeList[type])
}

#' Create a progress object
#' 
#' @param message Character: progress message
#' @param divisions Integer: number of divisions in the progress bar
#' @param global Shiny's global variable
#' 
#' @importFrom shiny isRunning Progress
#' @importFrom utils txtProgressBar
#' 
#' @return NULL (this function is used to modify the Shiny session's state or
#' internal hidden variables)
startProgress <- function(message, divisions,
                          global=if (isRunning()) sharedData else getHidden()) {
    display(message)
    if (isRunning()) {
        global$progress <- Progress$new()
        global$progress$set(message = message, value = 0)
    } else {
        global$progress <- txtProgressBar(style=3)
    }
    global$progress.divisions <- divisions
    return(invisible(global))
}

#' Update a progress object
#' 
#' @details If \code{divisions} is not NULL, a progress bar is started with the 
#' given divisions. If \code{value} is NULL, the progress bar will be 
#' incremented by one; otherwise, the progress bar will be incremented by the
#' integer given in value.
#' 
#' @inheritParams startProgress
#' @param value Integer: current progress value
#' @param max Integer: maximum progress value
#' @param detail Character: detailed message
#' @param console Boolean: print message to console? (TRUE by default)
#' 
#' @importFrom shiny isRunning Progress
#' @importFrom utils setTxtProgressBar
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
updateProgress <- function(message="Loading", value=NULL, max=NULL, detail=NULL,
                           divisions=NULL, 
                           global=if (isRunning()) sharedData else getHidden(),
                           console=TRUE) {
    if (!interactive()) return(NULL)
    if (!is.null(divisions)) {
        if (!isRunning()) # CLI version
            setHidden(startProgress(message, divisions, new.env()))
        else # GUI version
            startProgress(message, divisions, global)
        return(NULL)
    }

    divisions <- global$progress.divisions
    if (is.null(value)) {
        if (!isRunning()) { # CLI version
            currentValue <- global$progress$getVal()
            max   <- 1
        } else {
            currentValue <- global$progress$getValue()
            max          <- global$progress$getMax()
        }
        value <- currentValue + (max - currentValue)
    }
    amount <- ifelse(is.null(max), value/divisions, 1/max/divisions)

    # Print message to console
    if (console) {
        if (!is.null(detail) && !identical(detail, ""))
            display(paste(message, detail, sep=": "))
        else
            display(message)
    }

    # Increment progress
    if (!isRunning()) { # CLI version
        if (!is.null(global)) {
            value <- min(global$progress$getVal() + amount, 1)
            setTxtProgressBar(global$progress, value)
            setHidden(global)
        }
    } else { # GUI version
        if (is.null(detail)) detail <- ""
        global$progress$inc(amount=amount, message=message, detail=detail)
    }
    return(invisible(TRUE))
}

#' Close the progress even if there's an error
#' 
#' @param message Character: message to show in progress bar
#' @param global Global Shiny variable where all data is stored
#' 
#' @importFrom shiny isRunning Progress
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
closeProgress <- function(message=NULL, 
                          global=if (isRunning()) sharedData else getHidden()) {
    # Close the progress even if there's an error
    if (!is.null(message)) display(message)
    
    if (isRunning())
        global$progress$close()
    else
        close(global$progress)
}

#' Display characters in the command-line
#' 
#' @param char Character: message
#' @param timeStr Character: message when a \code{difftime} object is passed to
#' the \code{char} argument
#' 
#' @importFrom shiny isRunning
#' 
#' @return NULL (display message in command-line)
display <- function(char, timeStr="Time difference of") {
    if (!isRunning()) cat("", fill=TRUE)
    if (is(char, "difftime")) {
        message(timeStr, " ", format(unclass(char), digits=3), " ", 
                attr(char, "units"))
    } else {
        cat(char, fill=TRUE)
    }
}

#' Get number of significant digits
#' @param n Numeric: number to round
#' @importFrom shiny isolate
#' @return Formatted number with a given number of significant digits
signifDigits <- function(n) {
    return(isolate(formatC(n, getSignificant(), format="g")))
}

#' Round by the given number of digits
#' @param n Numeric: number to round
#' @return Formatted number with a given numeric precision
roundDigits <- function(n) {
    return(isolate(formatC(n, getPrecision(), format="f")))
}

#' Modified version of \code{shinyBS::bsModal}
#' 
#' \code{bsModal} is used within the UI to create a modal window. This allows to
#' modify the modal footer.
#' 
#' @inheritParams shinyBS::bsModal
#' @param footer UI set: List of elements to include in the footer
#' @param style Character: message style can be \code{warning}, \code{error}, 
#' \code{info} or \code{NULL}
#' @param size Character: Modal size (\code{small}, \code{default} or 
#' \code{large})
#' 
#' @importFrom shiny tagAppendAttributes
#' @importFrom shinyBS bsModal
#' 
#' @return HTML elements
bsModal2 <- function (id, title, trigger, ..., size=NULL, footer=NULL, 
                      style = NULL)  {
    if (is.null(size))
        modal <- bsModal(id, title, trigger, ...)
    else
        modal <- bsModal(id, title, trigger, ..., size=size)
    
    if (!is.null(style)) {
        style <- match.arg(style, c("info", "warning", "error"))
        modal[[3]][[1]][[3]][[1]][[3]][[1]] <-
            tagAppendAttributes(modal[[3]][[1]][[3]][[1]][[3]][[1]],
                                class=style)
    }
    
    modal[[3]][[1]][[3]][[1]][[3]][[3]] <-
        tagAppendChild(modal[[3]][[1]][[3]][[1]][[3]][[3]], footer)
    return(modal)
}

#' Enable or disable a tab from the \code{navbar}
#' @importFrom shinyjs disable addClass
#' @param tab Character: tab
#' @return NULL (this function is used to modify the Shiny session's state)
disableTab <- function(tab) {
    # Style item as disabled
    addClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
             class = "disabled")
    # Disable link itself
    disable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' @rdname disableTab
#' @importFrom shinyjs removeClass enable
enableTab <- function(tab) {
    # Style item as enabled
    removeClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
                class = "disabled")
    # Enable link itself
    enable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Create script for auto-completion of text input
#' 
#' Uses the JavaScript library \code{jquery.textcomplete}
#' 
#' @param id Character: input ID
#' @param words Character: words to suggest
#' @param novalue Character: string when there's no matching values
#' @param char Character to succeed accepted word
#'
#' @return HTML string with the JavaScript script prepared to run
#' 
#' @examples 
#' words <- c("tumor_stage", "age", "gender")
#' psichomics:::textSuggestions("textareaid", words)
textSuggestions <- function(id, words, novalue="No matching value", char=" ") {
    varId <- paste0(gsub("-", "_", id), "_words")
    var <- paste0(varId, ' = ["', paste(words, collapse = '", "'), '"];')
    
    js <- paste0('$("#', escape(id), '").textcomplete([{
                 match: /([a-zA-Z0-9_\\.]{1,})$/,
                 search: function(term, callback) {
                 var words = ', varId, ', sorted = [];
                 for (i = 0; i < words.length; i++) {
                 sorted[i] = fuzzy(words[i], term);
                 }
                 sorted.sort(fuzzy.matchComparator);
                 sorted = sorted.map(function(i) { return i.term; });
                 callback(sorted);
                 },
                 index: 1,
                 cache: true,
                 replace: function(word) {
                 return word + "', char ,'";
                 }}], { noResultsMessage: "', novalue, '"});')
    js <- HTML("<script>", var, js, "</script>")
    return(js)
}


#' Blend two HEX colours
#' 
#' @param colour1 Character: HEX colour
#' @param colour2 Character: HEX colour
#' @param colour1Percentage Character: percentage of colour 1 mixed in blended 
#' colour (default is 0.5)
#'
#' @source Code modified from \url{https://stackoverflow.com/questions/5560248}
#'
#' @return Character representing an HEX colour
#' @examples 
#' psichomics:::blendColours("#3f83a3", "#f48000")
blendColours <- function (colour1, colour2, colour1Percentage=0.5) {
    colour1 <- gsub("#", "", colour1)
    colour1 <- as.hexmode(colour1)
    R1 <- bitwShiftR(colour1, 16)
    G1 <- bitwAnd(bitwShiftR(colour1, 8), 0x00FF)
    B1 <- bitwAnd(colour1, 0x0000FF)
    
    colour2 <- gsub("#", "", colour2)
    colour2 <- as.hexmode(colour2)
    R2 <- bitwShiftR(colour2, 16)
    G2 <- bitwAnd(bitwShiftR(colour2, 8), 0x00FF)
    B2 <- bitwAnd(colour2, 0x0000FF)
    
    # Round to biggest integer if ending in .5
    mround <- function(x) trunc(x + 0.5)
    
    red   <- 0x1000000 + (mround((R2 - R1) * colour1Percentage) + R1) * 0x10000;
    green <- (mround((G2 - G1) * colour1Percentage) + G1) * 0x100;
    blue  <- mround((B2 - B1) * colour1Percentage) + B1;
    blended <- substr(as.hexmode(red + green + blue), 2, 16)
    return(paste0("#", blended));
}

#' Plot survival curves using Highcharts
#' 
#' @param object A survfit object as returned from the \code{survfit} function
#' @inheritDotParams highcharter::hc_add_series -hc -data
#' @param fun Name of function or function used to transform the survival curve:
#' \code{log} will put y axis on log scale, \code{event} plots cumulative events
#' (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard function (f(y) =
#' -log(y)), and \code{cloglog} creates a complimentary log-log survival plot
#' (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param markTimes Label curves marked at each censoring time? TRUE by default
#' @param symbol Symbol to use as marker (plus sign by default)
#' @param markerColor Colour of the marker ("black" by default); use NULL to use
#' the respective colour of each series
#' @param ranges Plot interval ranges? FALSE by default
#' @param rangesOpacity Opacity of the interval ranges (0.3 by default)
#' 
#' @method hchart survfit
#' @importFrom highcharter %>% hc_add_series highchart hc_tooltip hc_yAxis
#' hc_plotOptions fa_icon_mark JS
#' @importFrom stats setNames
#' 
#' @return \code{highcharter} object to plot survival curves
#' 
#' @examples
#' 
#' # Plot Kaplan-Meier curves
#' require("survival")
#' require("highcharter")
#' leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
#' hchart(leukemia.surv)
#' 
#' # Plot the cumulative hazard function
#' lsurv2 <- survfit(Surv(time, status) ~ x, aml, type='fleming') 
#' hchart(lsurv2, fun="cumhaz")
#' 
#' # Plot the fit of a Cox proportional hazards regression model
#' fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian)
#' ovarian.surv <- survfit(fit, newdata=data.frame(age=60))
#' hchart(ovarian.surv, ranges = TRUE)
hchart.survfit <- function(object, ..., fun = NULL, markTimes = TRUE,
                           symbol = "plus", markerColor = "black",
                           ranges = FALSE, rangesOpacity = 0.3) {
    groups <- NULL
    # Check if there are groups
    if (is.null(object$strata))
        strata <- c("Series 1" = length(object$time))
    else
        strata <- object$strata
    
    # Modify data according to functions (adapted from survival:::plot.survfit)
    if (is.character(fun)) {
        tfun <- switch(fun,
                       log = function(x) x,
                       event = function(x) 1 - x,
                       cumhaz = function(x) -log(x),
                       cloglog = function(x) log(-log(x)),
                       pct = function(x) x * 100,
                       logpct = function(x) 100 * x,
                       identity = function(x) x,
                       function(x) x)
    } else if (is.function(fun)) {
        tfun <- fun
    } else {
        tfun <- function(x) x
    }
    
    firsty <- tfun(1)
    object$surv <- tfun(object$surv)
    if (ranges && !is.null(object$upper)) {
        object$upper <- tfun(object$upper)
        object$lower <- tfun(object$lower)
    }
    
    # Prepare data
    data <- data.frame(x=object$time, y=object$surv,
                       up=object$upper, low=object$lower,
                       group=rep(names(strata), strata), 
                       stringsAsFactors = FALSE)
    # Data markers
    marker <- list(list(fillColor=markerColor, symbol=symbol, enabled=TRUE))
    if(markTimes)
        mark <- object$n.censor == 1
    else
        mark <- FALSE
    
    # Adjust Y axis range
    yValues <- object$surv
    ymin <- ifelse(min(yValues) >= 0, 0, min(yValues))
    ymax <- ifelse(max(yValues) <= 1, 1, max(yValues))
    
    hc <- highchart() %>%
        hc_tooltip(pointFormat="{point.y}") %>%
        hc_yAxis(min=ymin, max=ymax) %>%
        hc_plotOptions(line = list(marker = list(enabled = FALSE)))
    
    count <- 0
    
    # Process groups by columns (CoxPH-like) or in a single column
    if(!is.null(ncol(object$surv))) {
        groups <- seq(ncol(object$surv))
    } else {
        groups <- names(strata)
    }
    
    summ <- summary(object)$table
    for (name in groups) {
        if (!is.null(ncol(object$surv))) {
            df <- df[c("x", paste(c("y", "low", "up"), col, sep="."))]
            names(df) <- c("x", "y", "low", "up")
            submark <- mark
        } else {
            this <- data$group == name
            df <- data[this, ]
            submark <- mark[this]
        }
        
        # Add first value if there is no value for time at 0 in the data
        first <- NULL
        if (!0 %in% df$x) first <- list(list(x=0, y=firsty))
        
        # Mark events
        ls <- lapply(seq(nrow(df)), function(i) as.list(df[i, , drop=FALSE]))
        if (markTimes) ls[submark] <- lapply(ls[submark], c, marker=marker)
        
        if (is.matrix(summ))
            curveSumm <- summ[name, ]
        else
            curveSumm <- summ
        
        # Prepare group colours
        colour <- attr(object, "Colour")
        if (!is.null(colour))
            colour <- unname(colour[name])
        else
            colour <- JS("Highcharts.getOptions().colors[", count, "]")
        
        hc <- do.call(hc_add_series, c(list(
            hc, data=c(first, ls), step="left", name=name, zIndex=1,
            color=colour, ...), curveSumm))
        
        if (ranges && !is.null(object$upper)) {
            # Add interval range
            range <- lapply(ls, function(i) 
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=rangesOpacity,
                lineWidth=0, color=colour, ...)
        }
        count <- count + 1
    }
    
    return(hc)
}

#' Render a data table with sparkline HTML elements
#' 
#' @details This slightly modified version of \code{\link{renderDataTable}}
#' calls a JavaScript function to convert the sparkline HTML elements to
#' interactive Highcharts
#' 
#' @importFrom DT renderDataTable JS
#' 
#' @inheritDotParams shiny::renderDataTable -options -escape -env
#' @param options List of options to pass to \code{\link{renderDataTable}}
#' @return NULL (this function is used to modify the Shiny session's state)
renderDataTableSparklines <- function(..., options=NULL) {
    # Escape is set to FALSE to render the Sparkline HTML elements
    renderDataTable(..., escape=FALSE, env=parent.frame(n=1), options=c(
        list(drawCallback=JS("drawSparklines")), options))
}

#' Check unique rows of a data frame based on a set of its columns
#' 
#' @param data Data frame or matrix
#' @param ... Name of columns
#' 
#' @return Data frame with unique values based on set of columns
uniqueBy <- function(data, ...) {
    sub <- subset(data, select=c(...))
    uniq <- !duplicated(sub)
    return(data[uniq, ])
}

#' Add an exporting feature to a \code{highcharts} object
#' 
#' @param hc A \code{highcharts} object
#' @param fill Character: colour fill
#' @param text Character: button text
#' 
#' @importFrom highcharter hc_exporting JS
#' 
#' @return A \code{highcharts} object with an export button
export_highcharts <- function(hc, fill="transparent", text="Export") {
    export <- list(
        list(text="PNG image",
             onclick=JS("function () { 
                            this.exportChart({ type: 'image/png' }); }")),
        list(text="JPEG image",
             onclick=JS("function () { 
                            this.exportChart({ type: 'image/jpeg' }); }")),
        list(text="SVG vector image",
             onclick=JS("function () { 
                            this.exportChart({ type: 'image/svg+xml' }); }")),
        list(text="PDF document",
             onclick=JS("function () { 
                            this.exportChart({ type: 'application/pdf' }); }")),
        list(separator=TRUE),
        list(text="CSV document",
             onclick=JS("function () { this.downloadCSV(); }")),
        list(text="XLS document",
             onclick=JS("function () { this.downloadXLS(); }"))
    )
    
    hc_exporting(hc, enabled=TRUE,
                 formAttributes=list(target="_blank"),
                 buttons=list(contextButton=list(
                     text=text, theme=list(fill=fill),
                     menuItems=export)))
}

#' Create scatter plot
#' 
#' Create a scatter plot using \code{highcharter}
#' 
#' @param hc \code{Highchart} object
#' @param x Numeric: X axis
#' @param y Numeric: Y axis
#' @param z Numeric: Z axis to set the bubble size (optional)
#' @param label Character: data label for each point (optional)
#' @param showInLegend Boolean: show the data in the legend box? FALSE by
#' default
#' @inheritDotParams highcharter::hc_add_series -hc -data
#' 
#' @importFrom highcharter hc_add_series list_parse
#' 
#' @return \code{highcharter} object containing information for a scatter plot
hc_scatter <- function (hc, x, y, z=NULL, label=NULL, showInLegend=FALSE, ...) {
    df <- data.frame(x, y)
    if (!is.null(z)) df <- cbind(df, z=z)
    if (!is.null(label)) df <- cbind(df, label=label)
    
    args <- list(...)
    for (i in seq_along(args)) {
        if (!is.list(args[[i]]) && length(x) == length(args[[i]]) && 
            names(args[i]) != "name") {
            attr <- list(args[i])
            names(attr) <- names(args)[i]
            df <- cbind(df, attr)
            args[[i]] <- character(0)
        }
    }
    
    args <- Filter(length, args)
    ds <- list_parse(df)
    type <- ifelse(!is.null(z), "bubble", "scatter")
    
    if (!is.null(label))
        dlopts <- list(enabled=TRUE, format="{point.label}")
    else
        dlopts <- list(enabled=FALSE)
    
    do.call("hc_add_series", c(list(hc, data=ds, type=type, 
                                    showInLegend=showInLegend,
                                    dataLabels=dlopts), args))
}

#' Create an icon based on set operations
#' 
#' Based on the \code{\link[shiny]{icon}} function
#' 
#' @param name Character: icon name
#' @param class Character: additional classes to customise the icon element
#' @param ... Extra arguments for the icon HTML element
#' 
#' @importFrom shiny icon
#' @importFrom htmltools htmlDependency htmlDependencies htmlDependencies<-
#' 
#' @return Icon element
setOperationIcon <- function (name, class=NULL, ...) {
    if (length(list(...)) == 0) {
        style <- paste("font-size: 20px;", "line-height: 0;",
                       "vertical-align: bottom;", "display: inline-block;")
    } else {
        style <- NULL
    }
    
    prefix <- "set"
    iconClass <- ""
    if (!is.null(name)) 
        iconClass <- paste0(prefix, " ", prefix, "-", name)
    if (!is.null(class)) 
        iconClass <- paste(iconClass, class)
    iconTag <- tags$i(class=iconClass, style=style, ...)
    htmlDependencies(iconTag) <- htmlDependency(
        "set-operations", "1.0",
        c(href="set-operations"), stylesheet = "css/set-operations.css")
    return(iconTag)
}


# File browser dialog -----------------------------------------------------


#' Interactive folder selection using a native dialogue
#'
#' @param default Character: path to initial folder
#' @param caption Character: caption on the selection dialogue
#' @param multiple Boolean: allow to select multiple files?
#' @param directory Boolean: allow to select directories instead of files?
#' @param system Character: system name
#'
#' @details
#' For macOS, it uses an Apple Script to display a folder selection dialogue. 
#' With \code{default = NA}, the initial folder selection is determined by 
#' default behaviour of the "choose folder" Apple Script command.  Otherwise, 
#' paths are expanded with \link{path.expand}.
#'
#' In Windows, it uses either `utils::choose.files` or `utils::choose.dir`.
#'
#' @source Original code by wleepang:
#' \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @return A length one character vector, character NA if 'Cancel' was selected.
fileBrowser <- function(default=NULL, caption=NULL, multiple=FALSE,
                        directory=FALSE, system=Sys.info()['sysname']) {
    if (system == 'Darwin') {
        directory <- ifelse(directory, "folder", "file")
        multiple  <- ifelse(multiple, "with multiple selections allowed", "")
        
        if (!is.null(caption) && nzchar(caption)) {
            prompt <- sprintf("with prompt \\\"%s\\\"", caption)
        } else {
            prompt <- ""
        }
        
        # Default location
        if (!is.null(default) && nzchar(default)) {
            default <- sprintf("default location \\\"%s\\\"",
                               path.expand(default))
        } else {
            default <- ""
        }
        
        app  <- "path to frontmost application as text"
        args <- '-e "tell app (%s) to POSIX path of (choose %s %s %s %s)"'
        args <- sprintf(args, app, directory, multiple, prompt, default)
        
        suppressWarnings({
            path <- system2("osascript", args=args, stderr=TRUE)
        })
        
        # Return NA if user cancels the action
        if (!is.null(attr(path, "status")) && attr(path, "status")) {
            # user canceled
            return(NA)
        }
    } else if (system == 'Linux') {
        directory <- ifelse(directory, "--directory", "")
        multiple  <- ifelse(multiple, "--multiple", "")
        
        if (!is.null(caption) && nzchar(caption)) {
            prompt <- sprintf("--title='%s'", caption)
        } else {
            prompt <- ""
        }
        
        args <- " --file-selection %s %s %s"
        args <- sprintf(args, directory, multiple, prompt)
        
        suppressWarnings({
            path <- system2("zenity", args=args, stderr=TRUE)
        })
        
        # Return NA if user cancels the action
        if (!is.null(attr(path, "status")) && attr(path, "status")) {
            return(NA) # Cancelled by user
        }
        # Error: Gtk-Message: GtkDialog mapped without a transient parent
        if(length(path) == 2) path <- path[2]
    } else if (system == "Windows") {
        if (is.null(default)) default <- ""
        if (is.null(caption)) caption <- ""
        
        if (directory)
            path <- utils::choose.dir(default, caption)
        else
            path <- utils::choose.files(default, caption, multiple)
    }
    
    if (identical(path, "")) path <- NULL
    return(path)
}

#' File browser input
#'
#' Input to interactively select a file or directory on the server
#'
#' @param id Character: input identifier
#' @param label Character: input label (NULL to show no label)
#' @param value Character: initial value (paths are expanded via
#' \code{\link{path.expand}})
#' @param placeholder Character: placeholder when no file or folder is selected
#' @param info Boolean: add information icon for tooltips and popovers
#' @param infoFUN Function to use to provide information (e.g.
#' \code{shinyBS::bsTooltip} and \code{shinyBS::bsPopover})
#' @param infoPlacement Character: placement of the information (top, bottom,
#' right or left)
#' @param infoTitle Character: text to show as title of information
#' @param infoContent Character: text to show as content of information
#'
#' @details
#' To show the dialog for file input, the \code{\link{prepareFileBrowser}}
#' function needs to be included in the server logic.
#' 
#' This widget relies on \code{\link{fileBrowser}} to present an interactive
#' dialogue to users for selecting a directory on the local filesystem. 
#' Therefore, this widget is intended for shiny apps that are run locally - i.e. 
#' on the same system that files/directories are to be accessed - and not from 
#' hosted applications (e.g. from \url{https://www.shinyapps.io}).
#'
#' @importFrom shinyBS bsPopover bsTooltip
#' 
#' @source Original code by wleepang:
#' \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @return
#' A file browser input control that can be added to a UI definition.
#'
#' @seealso
#' \code{\link{updateFileBrowserInput}} and \code{\link{prepareFileBrowser}}
fileBrowserInput <- function(id, label, value=NULL, placeholder=NULL,
                             info=FALSE, infoFUN=NULL, infoPlacement="right",
                             infoTitle="", infoContent="") {
    if (!is.null(value) && !is.na(value)) value <- path.expand(value)
    if (is.null(placeholder)) placeholder <- ""
    
    check <- function (x, y) {
        # Based on shiny:::`%AND%`
        if (!is.null(x) && !is.na(x)) 
            if (!is.null(y) && !is.na(y)) 
                return(y)
        return(NULL)
    }
    
    if (info) {
        infoId   <- paste0(id, "-info")
        infoElem <- div(class="input-group-addon", id=infoId,
                        icon("question-circle"))
    } else {
        infoId   <- id
        infoElem <- NULL
    }
    
    infoTitle   <- gsub("\n", "", as.character(infoTitle),   fixed=TRUE)
    infoContent <- gsub("\n", "", as.character(infoContent), fixed=TRUE)
    if (identical(infoFUN, bsPopover)) {
        showInfo <- infoFUN(
            infoId, placement=infoPlacement, options=list(container="body"), 
            title=infoTitle, content=infoContent)
    } else if (identical(infoFUN, bsTooltip)) {
        showInfo <- infoFUN(
            infoId, placement=infoPlacement, options=list(container="body"), 
            title=infoTitle)
    } else {
        showInfo <- NULL
    }
    
    tagList(
        div(class='form-group fileBrowser-input-container',
            check(label, tags$label(label)),
            div(class='input-group shiny-input-container', style='width:100%;',
                div(class="input-group-btn",
                    div(class="btn btn-default fileBrowser-input",
                        id=sprintf("%sButton", id), 'Browse...')),
                tags$input(
                    id=id, value=value, type='text', readonly='readonly',
                    placeholder=placeholder,
                    class='form-control fileBrowser-input-chosen-dir'),
                infoElem)),
        showInfo)
}

#' Change the value of a fileBrowserInput on the client
#'
#' @param session Shiny session
#' @param id Character: input identifier
#' @param value Character: file or directory path
#' @param ... Additional arguments passed to \code{\link{fileBrowser}}. Only
#' used if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object. For \code{fileBrowserInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event. A directory
#' selection dialogue is not displayed.
#'
#' @source Original code by wleepang:
#' \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @return NULL (this function is used to modify the Shiny session's state)
updateFileBrowserInput <- function(session, id, ..., value=NULL) {
    if (is.null(value)) value <- fileBrowser(...)
    
    button <- sprintf("%sButton", id)
    session$sendInputMessage(button, list(path=value))
}

#' Prepare file browser dialogue and update the input's value accordingly to
#' selected file or directory
#'
#' @param session Shiny session
#' @param input Shiny input
#' @param id Character: input identifier
#' @inheritDotParams fileBrowser
#'
#' @return NULL (this function is used to modify the Shiny session's state)
prepareFileBrowser <- function(session, input, id, ...) {
    buttonId <- sprintf("%sButton", id)
    observeEvent(input[[buttonId]], {
        if (input[[buttonId]] > 0) { # Prevent execution on initial launch
            updateFileBrowserInput(session, id, ...)
        }
    })
}


# Interactive ggplot ------------------------------------------------------

#' Create HTML table from data frame or matrix
#' 
#' @param data Data frame or matrix
#' @param rownames Boolean: print row names?
#' @param colnames Boolean: print column names?
#' @param class Character: table class
#' @param style Character: table style
#' @param thead Boolean: add a \code{thead} tag to the first row?
#' 
#' @importFrom xtable xtable print.xtable
#' @importFrom shiny HTML
#' 
#' @return HTML elements
table2html <- function(data, rownames=TRUE, colnames=TRUE, class=NULL, 
                       style=NULL, thead=FALSE) {
    table <- xtable(data)
    table <- print(table, type="html", print.results=FALSE,
                   include.rownames=rownames, include.colnames=colnames)
    html <- HTML(table)
    
    if (thead)
        html <- gsub("(<tr>[[:space:]]*<th>.*</th>[[:space:]]*</tr>)",
                     "<thead>\\1</thead>", html)
    
    if (!is.null(class))
        class <- sprintf('class="%s"', paste(class, collapse=" "))
    
    if (!is.null(style))
        style <- sprintf('style="%s"', paste(style, collapse=" "))
    
    rep <- paste(class, style)
    if (length(rep) > 0) 
        html <- gsub("border=1", rep, html, fixed=TRUE)
    
    return(html)
}

#' Interface for interactive ggplot
#' 
#' @param id Character: identifier
#' 
#' @importFrom shiny tagList plotOutput brushOpts hoverOpts uiOutput 
#' actionButton
#' @importFrom shinyjs hidden
#' 
#' @return HTML elements
ggplotUI <- function(id) {
    idd <- function(str) paste(id, str, sep="-")
    plotId    <- idd("plot")
    tooltipId <- idd("tooltip")
    brushId   <- idd("brush")
    hoverId   <- idd("hover")
    resetId   <- idd("resetZoom")
    tagList(
        # Mimic Highcharts button to reset zoom level
        hidden(actionButton(
            resetId, "Reset zoom",
            style="font-size: 13px;",
            style="background-color: #f7f7f7;",
            style="border-color: #cccccc;",
            style="padding-bottom: 5px;", style="padding-top: 5px;",
            style="padding-left: 9px;", style="padding-right: 9px;",
            style="position: absolute;", style="z-index: 1;",
            style="top: 20px;", style="right: 35px;")),
        plotOutput(plotId,
                   brush=brushOpts(brushId, resetOnNew=TRUE),
                   hover=hoverOpts(hoverId, delay=50, delayType="throttle")),
        uiOutput(tooltipId))
}

#' Create the interface for the tooltip of a plot
#' 
#' @param df Data frame
#' @param hover Mouse hover information for a given plot as retrieved from
#' \code{\link[shiny]{hoverOpts}}
#' @param x Character: name of the variable used for the X axis
#' @param y Character: name of the variable used for the Y axis
#' 
#' @importFrom shiny tags nearPoints wellPanel
#' 
#' @return HTML elements
ggplotTooltip <- function(df, hover, x, y) {
    point <- nearPoints(df, hover, threshold=10, maxpoints=1, addDist=TRUE,
                        xvar=x, yvar=y)
    if (nrow(point) == 0) return(NULL)
    
    # Calculate point position inside the image as percent of total 
    # dimensions from left (horizontal) and from top (vertical)
    xDomain   <- hover$domain$right - hover$domain$left
    right_pct <- (hover$domain$right - hover$x) / xDomain
    yDomain   <- hover$domain$top - hover$domain$bottom
    top_pct   <- (hover$domain$top - hover$y) / yDomain
    
    # Calculate distance from left and bottom in pixels
    xRange   <- hover$range$right - hover$range$left
    right_px <- right_pct * xRange + 25
    yRange   <- hover$range$bottom - hover$range$top
    top_px   <- hover$range$top + top_pct * yRange + 2
    
    trItem <- function(key, value) tags$tr(tags$td(tags$b(key)), tags$td(value))
    
    thisPoint <- rownames(point)
    if ( areSplicingEvents(thisPoint) ) {
        event  <- parseSplicingEvent(thisPoint, pretty=TRUE)
        strand <- ifelse(event$strand == "+", "forward", "reverse")
        gene   <- paste(event$gene[[1]], collapse=" or ")
        type   <- trItem("Event type", event$type)
        coord  <- trItem(
            "Coordinates", 
            sprintf("chr %s: %s to %s (%s strand)", event$chrom,
                    event$pos[[1]][[1]], event$pos[[1]][[2]], strand))
    } else {
        gene  <- thisPoint
        type  <- NULL
        coord <- NULL
    }
    
    # Tooltip
    wellPanel(
        class="well-sm",
        style=paste0("position: absolute; z-index: 100;",
                     "background-color: rgba(245, 245, 245, 0.85); ",
                     "right:", right_px, "px; top:", top_px, "px;"),
        tags$table(class="table table-condensed", style="margin-bottom: 0;",
                   tags$thead( trItem("Gene", gene)),
                   tags$tbody( type, coord,
                               trItem(x, roundDigits(point[[x]])),
                               trItem(y, roundDigits(point[[y]])))))
}

#' Logic set to create an interactive ggplot
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param id Character: identifier
#' @param plot Character: plot expression (NULL renders no plot)
#' @inheritParams ggplotTooltip
#' 
#' @importFrom shiny renderPlot renderUI
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
ggplotServer <- function(input, output, id, plot=NULL, df=NULL, x=NULL, 
                         y=NULL) {
    idd <- function(str) paste(id, str, sep="-")
    output[[idd("plot")]] <- renderPlot(plot)
    
    if (is.null(plot)) {
        output[[idd("tooltip")]] <- renderUI(NULL)
    } else {
        output[[idd("tooltip")]] <- renderUI(
            ggplotTooltip(df, input[[idd("hover")]], x, y))
    }
}

#' @rdname ggplotServer
#' 
#' @note Insert \code{ggplotAuxSet} outside any observer (so it is only run 
#' once)
ggplotAuxServer <- function(input, output, id) {
    idd <- function(str) paste(id, str, sep="-")
    
    # Save zoom coordinates according to brushed area of the plot
    observe({
        brush <- input[[idd("brush")]]
        if (!is.null(brush)) {
            setZoom(id, brush)
            setSelectedPoints(id, NULL)
        }
    })
    
    # Toggle visibility of reset zoom button
    observe({
        zoom <- getZoom(id)
        if (is.null(zoom))
            hide(idd("resetZoom"))
        else
            show(idd("resetZoom"))
    })
    
    # Reset zoom when clicking the respective button
    observeEvent(input[[idd("resetZoom")]], {
        setZoom(id, NULL)
        setSelectedPoints(id, NULL)
    })
}