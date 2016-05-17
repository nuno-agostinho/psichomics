#' @import shiny shinyBS shinyjs ggplot2 highcharter survival miscTools
#' @include utils.R 
NULL

# Global variable with all the data of a session
sharedData <- reactiveValues()

# Get data from sharedData
getData <- reactive(sharedData$data)
getCategory <- reactive(sharedData$category)
getCategories <- reactive(names(getData()))
getCategoryData <- reactive(
    if(!is.null(getCategory())) getData()[[getCategory()]])
getClinicalData <- reactive(getCategoryData()[["Clinical data"]])
getJunctionQuantification <- reactive(
    getCategoryData()[["Junction quantification"]])
getInclusionLevels <- reactive(
    getCategoryData()[["Inclusion levels"]])
getGroupsFrom <- function(dataType, category = getCategory())
    sharedData[[paste(category, dataType, "groups", sep = "_")]] 
getClinicalMatchFrom <- function(dataType, category = getCategory())
    sharedData[[paste(category, dataType, "clinicalMatch", sep = "_")]] 

# Set data from sharedData (needs to be inside reactive functions)
setElement <- function(item, value) sharedData[[item]] <- value
setData <- function(value) setElement("data", value)
setCategory <- function(value) setElement("category", value)
setInclusionLevels <- function(value, category = getCategory())
    sharedData$data[[category]][["Inclusion levels"]] <- value
setGroupsFrom <- function(dataType, value, category = getCategory())
    setElement(paste(category, dataType, "groups", sep = "_"), value)
setClinicalMatchFrom <- function(dataType, value, category = getCategory())
    setElement(paste(category, dataType, "clinicalMatch", sep = "_"), value)

#' Create an identifier for a given object
#' 
#' @param ... Arguments to identify an object
#' 
#' @details To make an object's identifier unique, use the names of the module,
#' submodule, subsubmodule, ... as arguments before passing the name. Check the
#' example.
#' 
#' @return Character with underscores instead of spaces
#' 
#' @examples
#' psichomics:::objectId("Exploratory analysis", "PCA", "plot")
objectId <- function(...) {
    return(gsub(" ", "_", paste(...)))
}

# TODO(NunoA): increase allowed size and warn the user to wait for large files
# Refuse files with size greater than the specified
MB = 5000 # File size in MB
options(shiny.maxRequestSize = MB * 1024^2)

# TODO(NunoA): remove this (it's only for documentation purposes)
# options(shiny.trace=TRUE)

tabsFolder <- "R/"

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

#' Checks if a given script defines the given objects
#'
#' Loads the script into a new environment and checks if all the given objects
#' are present.
#'
#' @param script Character: file path to the script
#' @param check Character: objects to check
#' @param parentEnv Environment: enclosing environment inherited
#'
#' @return Environment with the loaded script if all the given objects are
#' present; otherwise, returns NULL
#' @export
checkObjects <- function(script, check, parentEnv = NULL) {
    if (is.null(parentEnv))
        env <- new.env()
    else
        env <- new.env(parent = parentEnv)
    sys.source(script, env)
    # Check for the given variables
    varsDefined <- vapply(check, exists, logical(1), envir = env)
    if (all(varsDefined)) return(env)
}

#' Sources scripts containing the given variables from a given folder
#'
#' @param folder Character: folder where the scripts are located
#' @param ... Extra parameters to be passed to list.files
#' @inheritParams checkObjects
#'
#' @return List of environments with sourced scripts
#' @export
sourceScripts <- function(folder, check, parentEnv = NULL, ...){
    files <- list.files(folder, full.names = TRUE, ...)
    # Get every script that defines the desired variables
    envs <- lapply(files, checkObjects, check, parentEnv)
    envs <- Filter(Negate(is.null), envs)
    return(envs)
}

#' Call a given function from valid scripts
#' 
#' Scripts from a given folder are checked to see if they have the given
#' objects. If they do, a given function will be called.
#' 
#' @note It's a good idea to check if the function is included in the script.
#'
#' @inheritParams sourceScripts
#' @param func Character: name of function to call
#' @param ... Arguments to pass to the given function call
#'
#' @return Variable from valid script
#' @export
callScriptsFunction <- function(func, ..., check = func, folder = "R/") {
    # Get scripts given the variables of interest
    scripts <- sourceScripts(folder, check)
    # Get a given variable from those scripts
    f <- lapply(scripts, "[[", func)
    # Remove nulls (needed?)
    f <- Filter(Negate(is.null), f)
    # Calls the function of each script with the given parameters
    loaded <- lapply(f, do.call, list(...))
    return(loaded)
}

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
    return(clinicalRows)
}

#' Start graphical interface of PSICHOMICS
#' 
#' @param ... Parameters to pass to the function runApp
#' 
#' @export
psichomics <- function(..., reload = FALSE) {
    if (reload) devtools::load_all()
    runApp(...)
}

#' Assign one group for each clinical patient
#' 
#' @param groups Matrix: clinical groups
#' @param patientsNumber Integer: total number of clinical patients
#' @param includeOuterGroup Boolean: join the patients that have no groups?
#' @param outerGroupName Character: name to give to outer group
#' @param allDataName Character: name to give in case there are no groups
#' 
#' @return Character vector where each element corresponds to the group of a
#' clinical patient
#' @export
groupPerPatient <- function(groups, patientsNumber, includeOuterGroup=FALSE, 
                            outerGroupName="(Outer data)",
                            allDataName="All data") {
    ## TODO(NunoA): join groups if a patient belongs to more than one group?
    finalGroups <- rep(NA, patientsNumber)
    rows <- groups[, "Rows", drop=FALSE]
    
    if (length(rows) == 0) return(rep(allDataName, patientsNumber))
    
    for (i in seq_along(rows))
        finalGroups[as.numeric(rows[[i]])] <- rownames(rows)[i]
    
    # Assign patients with no groups to the outer group
    if (includeOuterGroup) finalGroups[is.na(finalGroups)] <- outerGroupName
    
    return(finalGroups)
}

#' Simply show a modal
#' 
#' You can also use \code{errorModal} and \code{warningModal} to use template 
#' modals already stylised to show errors and warnings respectively.
#' 
#' @param session Current Shiny session
#' @param title Character: modal title
#' @param content Character: 
#' @param style Character: style of the modal header; NULL (default), danger, 
#' info or warning
#' @param iconName Character: FontAwesome icon name to appear with the title
#' @param footer List of interface elements: Custom modal footer
#' @param printMessage Boolean: print to console? TRUE by default
#' @param size Character: modal size can be "large", "small" (default) or NULL
#' (medium)
#' 
#' @export
showModal <- function(session, title, ..., style = NULL,
                      iconName = "exclamation-circle", footer = NULL,
                      printMessage = FALSE, size = NULL) {
    session$output[["globalModal"]] <- renderUI(
        bsModal2(id(title), style = style, div(icon(iconName), title),
                 trigger = NULL, size = size, ..., footer = NULL))
    toggleModal(session, id(title), toggle = "open")
    if (printMessage) print(content)
}

#' @rdname showModal
#' @export
errorModal <- function(session, title, ..., footer = NULL) {
    showModal(session, title, ..., footer, style = "danger", size = "small",
              printMessage = FALSE, iconName = "times-circle")
}

#' @rdname showModal
#' @export
warningModal <- function(session, title, ..., footer = NULL) {
    showModal(session, title, ..., footer, style = "warning", size = "small",
              printMessage = FALSE, iconName = "exclamation-circle")
}