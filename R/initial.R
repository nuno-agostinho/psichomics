#' @import shiny shinyBS shinyjs ggplot2 highcharter
NULL
#> NULL

library(highcharter)

# Global variable with all the data of a session
sharedData <- reactiveValues()

# Get data from sharedData
getData <- reactive(sharedData$data)
getCategory <- reactive(sharedData$category)
getCategories <- reactive(names(getData()))
getCategoryData <- reactive(
    if(!is.null(getCategory())) getData()[[getCategory()]])
getJunctionQuantification <- reactive(
    getCategoryData()[["Junction quantification"]])
getInclusionLevels <- reactive(sharedData[["psi"]])
getGroupsFrom <- function(dataType, category = getCategory())
    sharedData[[paste(category, dataType, "groups", sep = "_")]] 

# Set data from sharedData (needs to be inside reactive functions)
setElement <- function(item, value) sharedData[[item]] <- value
setData <- function(value) setElement("data", value)
setCategory <- function(value) setElement("category", value)
setInclusionLevels <- function(value) setElement("psi", value)
setGroupsFrom <- function(dataType, value, category = getCategory())
    setElement(paste(category, dataType, "groups", sep = "_"), value)

#' Create an identifier for a given object
#' 
#' @param ... Arguments to identify an object
#' 
#' @details To make an object's identifier unique, use the names of the module,
#' submodule, subsubmodule, ... as arguments before passing the name. Check the
#' example.
#' 
#' @return Character with underscores instead of spaces
#' @examples
#' objectId("Exploratory analysis", "PCA", "plot")
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
    if (all(varsDefined))
        return(env)
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
#' @param func Character: function to call
#' @param ... Arguments to pass to the given function call
#'
#' @return Variable from valid script
#' @export
callScriptsFunction <- function(func, ..., check = func, folder = "R/") {
    # Get scripts given the variables of interest
    scripts <- sourceScripts(folder, check)
    # Get a given variable from those script
    f <- lapply(scripts, "[[", func)
    # Remove nulls (needed?)
    f <- Filter(Negate(is.null), f)
    # Calls the function of each script with the given parameters
    loaded <- lapply(f, do.call, list(...))
    return(loaded)
}

#' Rename new vector to avoid duplicated values with old vector
#'
#' Renames duplicated values by adding an index
#'
#' @param check Character: values to rename if duplicated
#' @param comp Character: values to compare with
#'
#' @return Character vector with renamed values if duplicated; else, it
#' returns the usual values
#' @export
#'
#' @examples
#' renameDuplicated(check = c("Nuno", "Carolina", "Mariana", "Teresa"),
#'        comp = c("Nuno", "Lina", "Marie"))
renameDuplicated <- function(check, comp) {
    # If there's nothing to compare with, return the values
    if (length(comp) == 0) return(check)
    
    repeated <- check %in% comp
    uniq <- c(comp, check[!repeated])
    
    for (dup in check[repeated]) {
        # Locate matches (don't forget the counter)
        expr <- paste0(dup, " \\([0-9]+\\)|", dup)
        locate <- grep(expr, uniq, value = TRUE)
        
        # Get the maximum counter and add one
        counter <- sub(".* \\(([0-9]+)\\)", "\\1", locate)
        
        # Replace strings with 0
        counter[grep("^[0-9]*$", counter, inver =TRUE)] <- 0
        dup <- sprintf("%s (%i)", dup, max(as.numeric(counter)) + 1)
        
        # Append value to the unique group
        uniq <- c(uniq, dup)
    }
    return(uniq)
}