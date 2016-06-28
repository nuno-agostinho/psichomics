#' @import shiny utils
#' @include utils.R 
NULL

#' Get number of significant digits
#' @param n Numeric: number to round
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

#' Get junction quantification data of the data category
getJunctionQuantification <- reactive(
    getCategoryData()[["Junction quantification"]])

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

#' Create an identifier for a given object
#' 
#' To make an object's identifier unique, use the names of the module,
#' submodule, subsubmodule, ... as arguments before passing the name. Check the
#' example.
#' 
#' @param ... Arguments to identify an object
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

tabsFolder <- system.file("R", package="psichomics")

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
#' @param pattern Character: regular expression to filter file names
#' @param ... Extra parameters to be passed to list.files
#' @inheritParams checkObjects
#'
#' @return List of environments with sourced scripts
#' @export
sourceScripts <- function(folder, check, parentEnv=NULL, pattern=NULL, ...){
    files <- list.files(folder, pattern=pattern, full.names = TRUE, ...)
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
    # Get scripts given the variables of interest to show at primary level
    scripts <- sourceScripts(folder, check)
    primary <- vapply(scripts, function(i) isTRUE(i[["primary"]]), logical(1))
    # Get a given variable from those scripts
    f <- lapply(scripts[primary], "[[", func)
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
    ## TODO(NunoA): join groups if a patient belongs to more than one group?
    rows <- groups[, "Rows", drop=FALSE]
    if (length(rows) == 0) return(rep(allDataName, patients))
    
    finalGroups <- rep(NA, patients)
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
#' @inheritParams bsModal2
#' @param session Current Shiny session
#' @param iconName Character: FontAwesome icon name to appear with the title
#' @param printMessage Boolean: print to console? FALSE by default
#' 
#' @importFrom shinyBS toggleModal
#' @export
showModal <- function(session, title, ..., style = NULL,
                      iconName = "exclamation-circle", footer = NULL,
                      printMessage = FALSE, size = NULL,
                      modalId = "modal") {
    ns <- session$ns
    session$output[[modalId]] <- renderUI({
        bsModal2(ns(id(title)), style=style, div(icon(iconName), title),
                 trigger=NULL, size=size, ..., footer=footer)})
    toggleModal(session, id(title), toggle = "open")
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