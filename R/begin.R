#' @include utils.R 
NULL

# TODO(NunoA): increase allowed size and warn the user to wait for large files
# Refuse files with size greater than the specified
MB = 5000 # File size in MB
options(shiny.maxRequestSize = MB * 1024^2)

tabsFolder <- system.file("R", package="psichomics")

#' Get number of significant digits
#' @param n Numeric: number to round
#' @importFrom shiny isolate
signifDigits <- function(n) {
    return(isolate(formatC(n, getSignificant(), format="g")))
}

#' Round by the given number of digits
#' @param n Numeric: number to round
roundDigits <- function(n) {
    return(isolate(formatC(n, getPrecision(), format="f")))
}

# Global variable with all the data of a session
sharedData <- reactiveValues()

#' Get global data
getData <- reactive(sharedData$data)

#' Get number of cores
getCores <- reactive(sharedData$cores)

#' Get number of significant digits
getSignificant <- reactive(sharedData$significant)

#' Get number of decimal places
getPrecision <- reactive(sharedData$precision)

#' Get selected alternative splicing event
getEvent <- reactive(sharedData$event)

#' Get available data categories
getCategories <- reactive(names(getData()))

#' Get selected data category
getCategory <- reactive(sharedData$category)

#' Get data of selected data category
getCategoryData <- reactive(
    if(!is.null(getCategory())) getData()[[getCategory()]])

#' Get selected dataset
getActiveDataset <- reactive(sharedData$activeDataset)

#' Get clinical data of the data category
getClinicalData <- reactive(getCategoryData()[["Clinical data"]])

#' Get junction quantification data
#' @note Needs to be called inside reactive function
#' 
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' @return List of data frames of junction quantification
getJunctionQuantification <- function(category=getCategory()) {
    data <- getData()[[getCategory()]]
    match <- sapply(data, attr, "dataType") == "Junction quantification"
    if (any(match)) return(data[match])
}

#' Get inclusion leves of the selected data category
getInclusionLevels <- reactive(getCategoryData()[["Inclusion levels"]])

#' Get data from shared data
#' @param ... Arguments to identify a variable
#' @param sep Character to separate identifiers
getGlobal <- function(..., sep="_") sharedData[[paste(..., sep=sep)]]

#' Get the table of differential analyses of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @return Data frame of differential analyses
getDifferentialAnalyses <- function(category = getCategory())
    getGlobal(category, "differentialAnalyses")

#' Get the species of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @return Character value with the species
getSpecies <- function(category = getCategory())
    getGlobal(category, "species")

#' Get the assembly version of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @return Character value with the assembly version
getAssemblyVersion <- function(category = getCategory())
    getGlobal(category, "assemblyVersion")

#' Get groups from a given data type
#' @note Needs to be called inside reactive function
#' 
#' @param dataset Character: data set (e.g. "Clinical data")
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @return Matrix with groups of a given dataset
getGroupsFrom <- function(dataset, category = getCategory())
    getGlobal(category, dataset, "groups")

#' Get clinical matches from a given data type
#' @note Needs to be called inside reactive function
#' 
#' @param dataset Character: data set (e.g. "Junction quantification")
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @return Integer with clinical matches to a given dataset
getClinicalMatchFrom <- function(dataset, category = getCategory())
    getGlobal(category, dataset, "clinicalMatch")

#' Set element as globally accessible
#' @details Set element inside the global variable
#' @note Needs to be called inside reactive function
#' 
#' @param ... Arguments to identify a variable
#' @param value Any value to attribute to element
#' @param sep Character to separate identifier
setGlobal <- function(..., value, sep="_") {
    sharedData[[paste(..., sep=sep)]] <- value
}

#' Set data of the global data
#' @note Needs to be called inside reactive function
#' @param data Data frame or matrix to set as data
setData <- function(data) setGlobal("data", value=data)

#' Set number of cores
#' @param cores Character: number of cores
#' @note Needs to be called inside reactive function
setCores <- function(cores) setGlobal("cores", value=cores)

#' Set number of significant digits
#' @param significant Character: number of significant digits
#' @note Needs to be called inside reactive function
setSignificant <- function(significant) setGlobal("significant", value=significant)

#' Set number of decimal places
#' @param precision Character: number of decimal places
#' @note Needs to be called inside reactive function
setPrecision <- function(precision) setGlobal("precision", value=precision)

#' Set event
#' @param event Character: event
#' @note Needs to be called inside reactive function
setEvent <- function(event) setGlobal("event", value=event)

#' Set data category
#' @param category Character: data category
#' @note Needs to be called inside reactive function
setCategory <- function(category) setGlobal("category", value=category)

#' Set active dataset
#' @param dataset Character: dataset
#' @note Needs to be called inside reactive function
setActiveDataset <- function(dataset) setGlobal("activeDataset", value=dataset)

#' Set inclusion levels for a given data category
#' @note Needs to be called inside reactive function
#' 
#' @param value Data frame or matrix: inclusion levels
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setInclusionLevels <- function(value, category = getCategory())
    sharedData$data[[category]][["Inclusion levels"]] <- value

#' Set groups from a given data type
#' @note Needs to be called inside reactive function
#' 
#' @param dataset Character: data set (e.g. "Clinical data")
#' @param groups Matrix: groups of dataset
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setGroupsFrom <- function(dataset, groups, category = getCategory())
    setGlobal(category, dataset, "groups", value=groups)

#' Set the table of differential analyses of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param table Character: differential analyses table
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setDifferentialAnalyses <- function(table, category = getCategory())
    setGlobal(category, "differentialAnalyses", value=table)

#' Set the species of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param value Character: species
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setSpecies <- function(value, category = getCategory())
    setGlobal(category, "species", value=value)

#' Set the assembly version of a data category
#' @note Needs to be called inside reactive function
#' 
#' @param value Character: assembly version
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setAssemblyVersion <- function(value, category = getCategory())
    setGlobal(category, "assemblyVersion", value=value)

#' Set clinical matches from a given data type
#' @note Needs to be called inside reactive function
#' 
#' @param dataset Character: data set (e.g. "Clinical data")
#' @param matches Vector of integers: clinical matches of dataset
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
setClinicalMatchFrom <- function(dataset, matches, category = getCategory())
    setGlobal(category, dataset, "clinicalMatch", value=matches)

#' Trims whitespace from a word
#'
#' @param word Character to trim
#'
#' @return Character without whitespace
#' @export
#'
#' @examples
#' trimWhitespace("    hey   there     ")
#' trimWhitespace(c("pineapple    ", "one two three", " sunken    ship   "))
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
#' @export
rm.null <- function(v) Filter(Negate(is.null), v)

#' Escape symbols for use in regular expressions
#'
#' @param ... Characters to be pasted with no space
#' 
#' @return Escaped string
#' @export
escape <- function(...) {
    # return(gsub("/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]", "\\$&", string))
    return(gsub("(\\W)", "\\\\\\1", paste0(...)))
}

#' Rename vector to avoid duplicated values with comparison
#'
#' Renames values by adding an index to the end of duplicates.
#'
#' @param check Character: values to rename if duplicated
#' @param comp Character: values to compare with
#'
#' @return Character vector with renamed values if duplicated; else, it
#' returns the usual values. Doesn't return the comparator values.
#' @export
#'
#' @examples
#' renameDuplicated(check = c("blue", "red"), comp = c("green", "blue"))
renameDuplicated <- function(check, comp) {
    # If there's nothing to compare with, return the values
    if (length(comp) == 0) return(check)
    
    repeated <- check %in% comp
    uniq <- check[!repeated]
    
    for (dup in check[repeated]) {
        # Locate matches (don't forget the counter)
        all <- c(comp, uniq)
        expr <- paste0(escape(dup), " \\([0-9]+\\)|", escape(dup))
        locate <- grep(expr, all, value = TRUE)
        
        # Get the maximum counter and add one
        counter <- sub(".* \\(([0-9]+)\\)", "\\1", locate)
        
        # Replace strings with 0
        counter[grep("^[0-9]*$", counter, invert =TRUE)] <- 0
        dup <- sprintf("%s (%i)", dup, max(as.numeric(counter)) + 1)
        
        # Append value to the unique group
        uniq <- c(uniq, dup)
    }
    return(uniq)
}

#' Match given IDs with the clinical data
#'
#' @param ids Character: IDs of interest
#' @param clinical Matrix or data.frame: clinical data 
#'
#' @return Integer vector of the row number in clinical data corresponding to 
#' the given IDs (named with the ID)
#' @export
matchIdWithClinical <- function(ids, clinical) {
    # All IDs are lower case in the clinical data
    ids <- tolower(ids)
    
    # Get all possible identifiers starting in "tcga" from the clinical data
    idsIndex <- sapply(1:ncol(clinical), 
                       function(i) length(grep("^tcga", clinical[[i]])) != 0)
    clinicalIds <- clinical[, idsIndex]
    
    # Get the clinical data row corresponding to the given IDs
    clinicalRows <- rep(NA, length(ids))
    names(clinicalRows) <- ids
    for (col in 1:ncol(clinicalIds)) {
        # Check if any ID matches the current column of clinical IDs
        match <- sapply(ids, grep, clinicalIds[ , col], fixed = TRUE)
        
        # All matched IDs will save their respective rows
        clinicalRows[lapply(match, length) != 0] <- unlist(match)
    }
    # Remove non-matching IDs
    clinicalRows <- clinicalRows[!is.na(clinicalRows)]
    return(clinicalRows)
}

#' Get row names matching selected groups
#' @param selected Character: name of selected groups
#' @param groups Matrix: named groups with row indexes
#' 
#' @return Names of the matching rows
getMatchingRowNames <- function(selected, groups, matches) {
    # Get selected groups
    rows <- groups[selected, "Rows"]
    
    # Get names of the matching rows with the data
    ns <- names(matches[matches %in% unlist(rows)])
    ns <- unique(ns)
    return(ns)
}

#' Assign one group for each clinical patient
#' 
#' @param groups Matrix: clinical groups
#' @param patients Integer: total number of clinical patients
#' @param includeOuterGroup Boolean: join the patients that have no groups?
#' @param outerGroupName Character: name to give to outer group
#' @param allDataName Character: name to give in case there are no groups
#' 
#' @return Character vector where each element corresponds to the group of a
#' clinical patient
#' @export
groupPerPatient <- function(groups, patients, includeOuterGroup=FALSE, 
                            outerGroupName="(Outer data)",
                            allDataName="All data") {
    rows <- groups[, "Rows", drop=FALSE]
    if (length(rows) == 0) return(rep(allDataName, patients))
    
    all <- unlist(rows)
    names(all) <- rep(rownames(rows), sapply(rows, length))

    finalGroups <- rep(NA, patients)
    for (each in unique(all))
        finalGroups[each] <- paste(names(all[all == each]), collapse=", ")
    
    # Assign patients with no groups to the outer group
    if (includeOuterGroup) finalGroups[is.na(finalGroups)] <- outerGroupName
    
    return(finalGroups)
}

#' Show a modal
#' 
#' You can also use \code{errorModal} and \code{warningModal} to use template 
#' modals already stylised to show errors and warnings respectively.
#' 
#' @inheritParams bsModal2
#' @param session Current Shiny session
#' @param iconName Character: FontAwesome icon name to appear with the title
#' @param printMessage Boolean: print to console? FALSE by default
#' @param modalId Character: identifier
#' 
#' @importFrom shiny renderUI div icon
#' @importFrom shinyBS toggleModal
#' @seealso showAlert
#' @export
showModal <- function(session, title, ..., style = NULL,
                      iconName = "exclamation-circle", footer = NULL,
                      printMessage = FALSE, size = NULL,
                      modalId = "modal") {
    ns <- session$ns
    session$output[[modalId]] <- renderUI({
        bsModal2(ns("showModal"), style=style, div(icon(iconName), title),
                 trigger=NULL, size=size, ..., footer=footer)})
    toggleModal(session, "showModal", toggle = "open")
    if (printMessage) print(content)
}

#' @rdname showModal
#' @export
errorModal <- function(session, title, ..., size = "small", footer = NULL) {
    showModal(session, title, ..., footer=footer, style = "error", size = size,
              printMessage = FALSE, iconName = "times-circle")
}

#' @rdname showModal
#' @export
warningModal <- function(session, title, ..., size = "small", footer = NULL) {
    showModal(session, title, ..., footer=footer, style="warning", size = size,
              printMessage = FALSE, iconName = "exclamation-circle")
}

#' @rdname showModal
#' @export
infoModal <- function(session, title, ..., size = "small", footer = NULL) {
    showModal(session, title, ..., footer=footer, style = "info", size = size,
              printMessage = FALSE, iconName = "info-circle")
}

#' Show a modal
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
#' @export
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
    
    if (!is.null(title)) title <- h3(title)
    
    session$output[[alertId]] <- renderUI({
        tagList(
            div(title, id="myAlert", class="alert", class=style, class="fade",
                class=dismissable, role="alert", dismiss, ...),
            tags$script("window.setTimeout(
                        function() { $('#myAlert').addClass('now'); }, 100)")
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

#' Return the type of a given sample
#' 
#' @param sample Character: ID of the sample
#' @param filename Character: path to RDS file containing corresponding type
getSampleTypes <- function(sample, 
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
#' 
#' @export
startProgress <- function(message, divisions, global = sharedData) {
    print(message)
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
#' 
#' @export
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
    global$progress$inc(amount = amount, message = message, detail = detail)
    
    if (!console)
        return(invisible(TRUE))
    
    # Print message to console
    if (!is.null(detail))
        print(paste(message, detail, sep=": "))
    else
        print(message)
    return(invisible(TRUE))
}

#' Close the progress even if there's an error
#' 
#' @param message Character: message to show in progress bar
#' @param global Global Shiny variable where all data is stored
#' 
#' @export
closeProgress <- function(message=NULL, global = sharedData) {
    # Close the progress even if there's an error
    if (!is.null(message)) print(message)
    global$progress$close()
}