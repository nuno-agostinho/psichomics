## Auxiliary functions used throughout the program

# Print how to start the graphical interface when attaching the package
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Start the visual interface by running the function ",
                          "psichomics()")
}

#' Get psichomics file inside a given directory
#' @param ... character vectors, specifying subdirectory and file(s) within some
#' package. The default, none, returns the root of the package. Wildcards are
#' not supported.
#' @return Loaded file
insideFile <- function(...) {
    return(system.file(..., package="psichomics"))
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

#' Parse an alternative splicing event based on a given identifier
#' 
#' @param event Character: event identifier
#' 
#' @return Parsed event
#' @export
#' @examples 
#' events <- c("SE_1_-_123_456_789_1024_TST",
#'             "MX_3_+_473_578_686_736_834_937_HEY/YOU")
#' parseSplicingEvent(events)
parseSplicingEvent <- function(event) {
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
    parsed$pos    <- lapply(seq_along(event), 
                            function(i) event[[i]][4:lenMinus1[[i]]])
    parsed$pos    <- lapply(parsed$pos, 
                            function(i) range(as.numeric(i)))
    
    parsed[,1] <- NULL
    return(parsed)
}

parseEvent <- parseSplicingEvent

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

#' Check if a number is whole
#' 
#' @param x Object to be tested
#' @param tol Numeric: tolerance used for comparison
#' 
#' @return TRUE if number is whole; otherwise, FALSE
is.whole <- function(x, tol=.Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
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
#' returns the usual values. It doesn't return the comparator values.
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
#' @param ... Extra parameters to pass to \code{actionButton}
#' @param class Character: class
#' 
#' @importFrom shinyjs hidden
#' @importFrom shiny actionButton
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
#' @param time POSIXct: start time needed to show the interval time (if NULL, 
#' the time interval is not displayed)
#' @param closeProgressBar Boolean: close progress bar? TRUE by default
#' 
#' @importFrom shinyjs hide
#' @return NULL (this function is used to modify the Shiny session's state)
endProcess <- function(id, time=NULL, closeProgressBar=TRUE) {
    enable(id)
    hide(paste0(id, "Loading"))
    if (closeProgressBar) closeProgress()
    if (!is.null(time)) 
        message("Process finished in ", format(Sys.time() - time))
}

#' Match given sample identifiers and return the respective row in clinical data
#'
#' @param sampleId Character: sample identifiers
#' @param clinical Matrix or data.frame: clinical data 
#' @param prefix Character: prefix to search for in clinical data
#' @param lower Boolean: convert samples to lower case? TRUE by default
#' @param rmNoMatches Boolean: remove non-matching identifiers
#'
#' @return Integer vector of the row number in clinical data corresponding to 
#' the given IDs (named with the ID)
#' @export
#' @examples 
#' samples <- c("ABC", "DEF", "GHI", "JKL", "MNO")
#' clinical <- data.frame(patient=paste0("patient-", samples),
#'                        samples=tolower(samples))
#' getPatientFromSample(samples, clinical, prefix="")
getPatientFromSample <- function(sampleId, clinical, prefix="^tcga", 
                                 lower=TRUE, rmNoMatches=TRUE) {
    # All IDs are lower case in the clinical data
    if (lower) sampleId <- tolower(sampleId)
    
    # Get all possible identifiers starting in "tcga" from the clinical data
    idsIndex <- vapply(1:ncol(clinical), 
                       function(i) length(grep(prefix, clinical[[i]])) != 0,
                       logical(1))
    clinicalIds <- clinical[, idsIndex, drop=FALSE]
    
    # Get the clinical data row corresponding to the given IDs
    clinicalRows <- rep(NA, length(sampleId))
    names(clinicalRows) <- sampleId
    for (col in 1:ncol(clinicalIds)) {
        # Check if any ID matches the current column of clinical IDs
        match <- sapply(sampleId, grep, clinicalIds[ , col], fixed = TRUE)
        
        # All matched IDs will save their respective rows
        clinicalRows[lapply(match, length) != 0] <- unlist(match)
    }
    # Remove non-matching identifiers
    if (rmNoMatches) clinicalRows <- clinicalRows[!is.na(clinicalRows)]
    return(clinicalRows)
}

#' Search samples in the clinical dataset and return the ones matching the given
#' index
#' 
#' @inheritParams getPatientFromSample
#' @param index Numeric or list of numeric: patient row indexes
#' @param samples Character: samples
#' @param clinical Data frame or matrix: clinical dataset
#' @param upper Boolean: convert identifiers to upper case? TRUE by default
#' @param rm.NA Boolean: remove NAs? TRUE by default
#' 
#' @return Names of the matching rows
#' @export
#' 
#' @examples 
#' samples <- c("ABC", "DEF", "GHI", "JKL", "MNO")
#' clinical <- data.frame(patient=paste0("patient-", samples), 
#'                        samples=tolower(samples))
#' getMatchingSamples(c(1, 4), samples, clinical, prefix="")
getMatchingSamples <- function(index, samples, clinical, upper=TRUE, 
                               rm.NA=TRUE, prefix="^tcga") {
    patient <- rep(NA, nrow(clinical))
    p <- getPatientFromSample(samples, clinical, prefix=prefix)
    patient[p] <- names(p)
    
    if (is.list(index)) {
        match <- lapply(index, function(i) {
            res <- patient[i]
            if (upper) res <- toupper(res)
            if (rm.NA) res <- res[!is.na(res)]
            return(res)
        })
    } else {
        match <- patient[index]
        if (upper) match <- toupper(match)
        if (rm.NA) match <- match[!is.na(match)]
    }    
    return(match)
}

#' Assign one group to each patient
#' 
#' @param groups List of integers: clinical groups
#' @param patients Integer: total number of clinical patients (remaining 
#' patients will be filled with missing values)
#' @param includeOuterGroup Boolean: join the patients that have no groups?
#' @param outerGroupName Character: name to give to outer group
#' 
#' @return Character vector where each element corresponds to the group of a
#' clinical patient
#' @export
#' 
#' @examples 
#' groups <- list(1:3, 4:7, 8:10)
#' names(groups) <- paste("Stage", 1:3)
#' groupPerPatient(groups)
groupPerPatient <- function(groups, patients, includeOuterGroup=FALSE, 
                            outerGroupName="(Outer data)") {
    if (length(groups) == 0) return(rep("Single group", patients))
    
    all <- unlist(groups)
    names(all) <- rep(names(groups), sapply(groups, length))
    
    finalGroups <- rep(NA, patients)
    for (each in unique(all))
        finalGroups[each] <- paste(names(all[all == each]), collapse=", ")
    
    # Assign patients with no groups to the outer group
    if (includeOuterGroup) finalGroups[is.na(finalGroups)] <- outerGroupName
    
    return(finalGroups)
}

#' Assign one group to each sample
#' 
#' @param groups List of characters: list of samples
#' @param samples Character: all available samples
#' @param includeOuterGroup Boolean: join the patients that have no groups?
#' @param outerGroupName Character: name to give to outer group
#' 
#' @return Character vector where each element corresponds to the group of a
#' sample
#' @export
#' 
#' @examples 
#' groups <- list(letters[1:3], letters[10:12], letters[5:8])
#' names(groups) <- paste("Stage", 1:3)
#' samples <- letters
#' groupPerSample(groups, samples)
groupPerSample <- function(groups, samples, includeOuterGroup=FALSE, 
                           outerGroupName="(Outer data)") {
    if (length(groups) == 0) return(rep("Single group", length(samples)))
    
    all <- unlist(groups)
    names(all) <- rep(names(groups), sapply(groups, length))
    
    finalGroups <- rep(NA, length(samples))
    for (each in unique(all)) {
        i <- match(each, samples)
        finalGroups[i] <- paste(names(all[all == each]), collapse=", ")
    }
    
    # Assign patients with no groups to the outer group
    if (includeOuterGroup) finalGroups[is.na(finalGroups)] <- outerGroupName
    
    return(finalGroups)
}

#' Style and show a modal
#' 
#' You can also use \code{errorModal} and \code{warningModal} to use template 
#' modals already stylised to show errors and warnings respectively.
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
#' @seealso showAlert
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
    if (echo) cat(content, fill=TRUE)
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

#' Show an alert
#' 
#' You can also use \code{errorAlert} and \code{warningAlert} to use template 
#' alerts already stylised to show errors and warnings respectively.
#' 
#' @param session Shiny session
#' @param ... Arguments to render as elements of alert
#' @param title Character: title of the alert (optional)
#' @param style Character: style of the alert ("alert-danger", "alert-warning" 
#' or NULL)
#' @param dismissable Boolean: is the alert dismissable? TRUE by default
#' @param alertId Character: alert identifier
#' 
#' @seealso showModal
#' @importFrom shiny span h3 renderUI div tagList
#' @return NULL (this function is used to modify the Shiny session's state)
showAlert <- function(session, ..., title=NULL, style=NULL, dismissable=TRUE, 
                      alertId="alert") {
    ns <- session$ns
    
    if (dismissable) {
        dismissable <- "alert-dismissible"
        dismiss <- tags$button(type="button", class="close",
                               "data-dismiss"="alert", "aria-label"="Close",
                               span("aria-hidden"="true", "\u00D7"))
    } else {
        dismissable <- NULL
        dismiss <- NULL
    }
    
    if (!is.null(title)) title <- h4(title)
    
    output <- session$output
    output[[alertId]] <- renderUI({
        tagList(
            div(title, id="myAlert", class="alert", class=style, role="alert",
                class="animated bounceInUp", class=dismissable, dismiss, ...)
        )
    })
}

#' @rdname showAlert
errorAlert <- function(session, ..., title=NULL, dismissable=TRUE,
                       alertId="alert") {
    showAlert(session, ..., style="alert-danger", title=title, 
              dismissable=dismissable, alertId=alertId)
}

#' @rdname showAlert
warningAlert <- function(session, ..., title=NULL, dismissable=TRUE,
                         alertId="alert") {
    showAlert(session, ..., style="alert-warning", title=title,
              dismissable=dismissable, alertId=alertId)
}

removeAlert <- function(output, alertId="alert") {
    output[[alertId]] <- renderUI(NULL)
}

#' Sample variance by row
#' 
#' Calculate the sample variance of each row in the given matrix
#' 
#' @param x Matrix
#' @param na.rm Boolean: should the NAs be ignored? FALSE by default
#' 
#' @return Variance for each row
rowVar <- function (x, na.rm = FALSE) {
    means <- rowMeans(x, na.rm = na.rm)
    meansSqDev <- (x - means)^2
    squaresSum <- rowSums(meansSqDev, na.rm = na.rm)
    nas <- rowSums(is.na(x))
    return(squaresSum/(ncol(x) - nas - 1))
}

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
#' @importFrom shiny Progress
#' @return NULL (this function is used to modify the Shiny session's state)
startProgress <- function(message, divisions, global = sharedData) {
    cat(message, fill=TRUE)
    global$progress.divisions <- divisions
    global$progress <- Progress$new()
    global$progress$set(message = message, value = 0)
}

#' Update a progress object
#' 
#' @details If \code{divisions} isn't NULL, a progress bar is started with the 
#' given divisions. If \code{value} is NULL, the progress bar will be 
#' incremented by one; otherwise, the progress bar will be incremented by the
#' integer given in value.
#' 
#' @inheritParams startProgress
#' @param value Integer: current progress value
#' @param max Integer: maximum progress value
#' @param detail Character: detailed message
#' @param console Boolean: print message to console? (TRUE by default)
#' @return NULL (this function is used to modify the Shiny session's state)
updateProgress <- function(message="Hang in there", value=NULL, max=NULL,
                           detail=NULL, divisions=NULL, 
                           global=sharedData, console=TRUE) {
    if (!is.null(divisions)) {
        startProgress(message, divisions, global)
        return(NULL)
    }
    divisions <- global$progress.divisions
    if (is.null(value)) {
        value <- global$progress$getValue()
        max <- global$progress$getMax()
        value <- value + (max - value)
    }
    amount <- ifelse(is.null(max), value/divisions, 1/max/divisions)
    if (is.null(detail)) detail <- "" # Force detail to not clean
    global$progress$inc(amount = amount, message = message, detail = detail)
    
    if (!console)
        return(invisible(TRUE))
    
    # Print message to console
    if (!is.null(detail))
        cat(paste(message, detail, sep=": "), fill=TRUE)
    else
        cat(message, fill=TRUE)
    return(invisible(TRUE))
}

#' Close the progress even if there's an error
#' 
#' @param message Character: message to show in progress bar
#' @param global Global Shiny variable where all data is stored
#' @return NULL (this function is used to modify the Shiny session's state)
closeProgress <- function(message=NULL, global = sharedData) {
    # Close the progress even if there's an error
    if (!is.null(message)) cat(message, fill=TRUE)
    global$progress$close()
}

#' Get number of significant digits
#' @param n Numeric: number to round
#' @importFrom shiny isolate
#' @return Formatted number with a given number of significant digits
signifDigits <- function(n) {
    return(isolate(formatC(n, getSignificant(), format="g")))
}

#' Round by the given number of digits
#' @param n Numeric: number to roundhf
#' @return Formatted number with a given numeric precision
roundDigits <- function(n) {
    return(isolate(formatC(n, getPrecision(), format="f")))
}

#' Modified version of shinyBS::bsModal
#' 
#' bsModal is used within the UI to create a modal window. This allows to use
#' the footer.
#' 
#' @inheritParams shinyBS::bsModal
#' @param footer UI set: List of elements to include in the footer
#' @param style Character: message style can be "warning", "error", "info" or 
#' NULL
#' @param size Character: Modal size ("small", "default" or "large")
#' 
#' @importFrom shiny tagAppendAttributes
#' @importFrom shinyBS bsModal
#' 
#' @return HTML element to create a modified modal
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

#' Disable a tab from the navbar
#' @importFrom shinyjs disable addClass
#' @param tab Character: tab to disable
#' @return NULL (this function is used to modify the Shiny session's state)
disableTab <- function(tab) {
    # Style item as disabled
    addClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
             class = "disabled")
    # Disable link itself
    disable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Enable a tab from the navbar
#' @importFrom shinyjs removeClass enable
#' @param tab Character: tab to enable
#' @return NULL (this function is used to modify the Shiny session's state)
enableTab <- function(tab) {
    # Style item as enabled
    removeClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
                class = "disabled")
    # Enable link itself
    enable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Create script for autocompletion of text input
#' 
#' Uses the JavaScript library jquery.textcomplete
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

#' Plot survival curves using Highcharts
#' 
#' @param object A survfit object as returned from the \code{survfit} function
#' @param ... Extra parameters to pass to \code{hc_add_series} function
#' @param fun Name of function or function used to transform the survival curve:
#' \code{log} will put y axis on log scale, \code{event} plots cumulative events
#' (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard function (f(y) =
#' -log(y)), and \code{cloglog} creates a complimentary log-log survival plot
#' (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param markTimes Label curves marked at each censoring time? TRUE by default
#' @param symbol Symbol to use as marker (plus sign by default)
#' @param markerColor Color of the marker ("black" by default); use NULL to use
#' the respective color of each series
#' @param ranges Plot interval ranges? FALSE by default
#' @param rangesOpacity Opacity of the interval ranges (0.3 by default)
#' 
#' @method hchart survfit
#' @importFrom highcharter %>% hc_add_series highchart hc_tooltip hc_yAxis
#' hc_plotOptions fa_icon_mark JS
#' @importFrom stats setNames
#' 
#' @return Highcharts object to plot survival curves
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
        if (!0 %in% df$x)
            first <- list(list(x=0, y=firsty))
        else
            first <- NULL
        
        # Mark events
        ls <- lapply(seq(nrow(df)), function(i) as.list(df[i, , drop=FALSE]))
        if (markTimes)
            ls[submark] <- lapply(ls[submark], c, marker=marker)
        
        if (is.matrix(summ))
            curveSumm <- summ[name, ]
        else
            curveSumm <- summ
        
        hc <- do.call(hc_add_series, c(list(
            hc, data=c(first, ls), step="left", name=name, zIndex=1,
            color=JS("Highcharts.getOptions().colors[", count, "]"), ...),
            curveSumm))
        
        if (ranges && !is.null(object$upper)) {
            # Add interval range
            range <- lapply(ls, function(i) 
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=rangesOpacity, 
                lineWidth=0,
                color=JS("Highcharts.getOptions().colors[", count, "]"),
                ...)
        }
        count <- count + 1
    }
    
    return(hc)
}

#' Render a data table with Sparkline HTML elements
#' 
#' @details This slightly modified version of \code{\link{renderDataTable}}
#' calls a JavaScript function to convert the Sparkline HTML elements to
#' interactive Highcharts
#' 
#' @importFrom DT renderDataTable JS
#' 
#' @param ... Arguments to pass to \code{\link{renderDataTable}}
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
#' @param y Numeric: position
#' @param verticalAlign Character: vertical alignment
#' @param fill Character: colour fill
#' @param text Character: button text
#' 
#' @return A \code{highcharts} object with an export button
export_highcharts <- function(hc, y=-45, verticalAlign="bottom", 
                              fill="transparent", text="Export") {
    hc_exporting(hc, enabled=TRUE,
                 formAttributes = list(target = "_blank"),
                 buttons=list(contextButton=list(text=text, y=y,
                                                 verticalAlign=verticalAlign, 
                                                 theme=list(fill=fill))))
}

#' Modified function of highcharter::hc_add_series_scatter
#' 
#' @inheritParams highcharter::hc_add_series_scatter
#' 
#' @importFrom dplyr data_frame mutate
#' @importFrom assertthat assert_that
#' @importFrom highcharter colorize list_parse hc_add_series
#' 
#' @return Highchart object containing information for a scatter plot
hc_add_series_scatter <- function (hc, x, y, z = NULL, color = NULL, 
                                   label = NULL, showInLegend = FALSE, ...) 
{
    assert_that(length(x) == length(y), is.numeric(x), is.numeric(y))
    df <- data_frame(x, y)
    if (!is.null(z)) {
        assert_that(length(x) == length(z))
        df <- df %>% mutate(z = z)
    }
    if (!is.null(color)) {
        assert_that(length(x) == length(color))
        cols <- colorize(color)
        df <- df %>% mutate(valuecolor = color, color = cols)
    }
    if (!is.null(label)) {
        assert_that(length(x) == length(label))
        df <- df %>% mutate(label = label)
    }
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
    if (!is.null(label)) {
        dlopts <- list(enabled = TRUE, format = "{point.label}")
    }
    else {
        dlopts <- list(enabled = FALSE)
    }
    do.call("hc_add_series", c(list(hc, data = ds, type = type, 
                                    showInLegend = showInLegend, dataLabels = dlopts), args))
}