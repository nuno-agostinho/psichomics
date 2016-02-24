#' @import shiny shinyBS shinyjs ggplot2

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


#' Checks the format of a file
#'
#' @param format Environment: format of the file
#' @param head Data.frame: head of the file to check
#'
#' @return TRUE if the file is of the given format; otherwise, returns FALSE
#' @export
checkFileFormat <- function(format, head) {
    ## TODO(NunoA): account for comments
    # Transpose data
    if (!is.null(format$transpose) && format$transpose) head <- t(head)
    
    # Check if header is comparable
    if (ncol(head) < length(format$header)) return(FALSE)
    
    # Check if headers match
    headerCheck <- ifelse(!is.null(format$headerCheck), format$headerCheck, 1)
    fileHeader <- head[headerCheck, 1:length(format$header)]
    valid <- all(fileHeader == format$header)
    return(valid)
}

#' Loads a file according to its format
#' 
#' @inheritParams checkFileFormat
#' @param file Character: file to load
#' 
#' @details The resulting data frame includes the attribute "tablename" with the
#' name of the data frame
#' 
#' @return Data frame with the loaded file
#' @export
loadFile <- function(format, file) {
    ## TODO(NunoA): account for the comment character
    delim <- ifelse(!is.null(format$delim), format$delim, "\t")
    loaded <- readr::read_delim(file, delim=delim, col_names = FALSE)
    
    # Transpose data
    if (!is.null(format$transpose) && format$transpose)
        loaded <- data.frame(t(loaded),
                             stringsAsFactors = FALSE,
                             row.names = NULL)
    
    # Column names
    headerUse <- ifelse(!is.null(format$headerUse), format$headerUse, 1)
    names(loaded) <- loaded[headerUse, ]
    
    # Filter in only content of interest
    contentStart <- ifelse(!is.null(format$contentStart),
                           format$contentStart, 2)
    loaded <- loaded[contentStart:nrow(loaded), ]
    attr(loaded, "tablename") <- format$tablename
    return(loaded)
}

#' Parse file given a folder with recognised formats
#' 
#' Tries to recognise the file format and parses the content of the given file
#' accordingly.
#' 
#' @param file Character: file to parse
#' @param formatsFolder Character: folder with recognised file formats
#' 
#' @details The resulting data frame includes the attribute "tablename" with the
#' name of the data frame
#' 
#' @return Data frame with the contents of the given file if the file format is
#' recognised; otherwise, returns NULL
parseValidFile <- function(file, formatsFolder) {
    # Get all available formats information
    formats <- sourceScripts(formatsFolder,
                             c("name", "header"))
    
    # The number of rows to read will be the maximum value asked by all the file
    # formats; if no format aks for a specific number of rows, the default is 6
    headRows <- lapply(formats, "[[", "header_rows")
    headRows <- unlist(rm.null(headRows))
    headRows <- ifelse(!is.null(headRows), max(headRows), 6)
    
    ## TODO(NunoA): readr can't be used because it gives an error for
    ## col_names=FALSE if the first line isn't the same data type of the rest
    ## of the file; e.g. if the first line is character and the rest is integer
    head <- read.delim(file, header = FALSE, nrows = 6, stringsAsFactors = F)
    
    # Check if the file is recognised by at least one file format
    recognised <- lapply(formats, checkFileFormat, head)
    recognised <- unlist(recognised)
    
    if (sum(recognised) > 1) {
        ## TODO(NunoA): If more than one format is recognised, check if the
        # filename is helpful in distinguishing the file
        stop("Error: more than one file format was recognised.")
        ## TODO(NunoA): if there is still a conflict, ask the user which file
        ## format to use
    } else if (sum(recognised) == 1) {
        format <- formats[recognised][[1]]
        loaded <- loadFile(format, file)
        return(loaded)
    } 
}