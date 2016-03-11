#' @import shiny shinyBS shinyjs ggplot2
NULL
#> NULL

# Global variable with all the data inside combos!
shared.data <- reactiveValues(combos = list())

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
#'
#' @return Environment with the loaded script if all the given objects are
#' present; otherwise, returns NULL
#' @export
checkObjects <- function(script, check) {
    env <- new.env()
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
sourceScripts <- function(folder, check, ...){
    files <- list.files(folder, full.names = TRUE, ...)
    # Get every script that defines the desired variables
    envs <- lapply(files, checkObjects, check)
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