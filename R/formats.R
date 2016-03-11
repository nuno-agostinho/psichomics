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
    if (!is.null(format$transpose) && format$transpose)
        head <- t(head)
    
    lenCheck <- length(format$check)
    # See whether using row or column to check format
    if (is.null(format$rowCheck) | !format$rowCheck) {
        # Check if length is not satisfatory
        if (nrow(head) < lenCheck) return(FALSE)
        # Select wanted row and check for a match
        return(all(head[1:lenCheck, format$checkIndex] == format$check))
    } else {
        # Check if length is not satisfatory
        if (ncol(head) < lenCheck) return(FALSE)
        # Select wanted column and check for a match
        return(all(head[format$checkIndex, 1:lenCheck] == format$check))
    }
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
    skip <- ifelse(!is.null(format$skip), format$skip, 0)
    
    loaded <- data.table::fread(file, sep = delim, header = FALSE,
                                stringsAsFactors = FALSE, data.table = FALSE,
                                skip = skip)
    
    # Transpose data
    if (!is.null(format$transpose) && format$transpose)
        loaded <- data.frame(t(loaded), stringsAsFactors = FALSE, 
                             row.names = NULL)
    
    # Remove duplicated rows
    if (!is.null(format$unique) && format$unique) loaded <- unique(loaded)
    
    # Add column names
    if (!is.null(format$colNames)) {
        if (skip != 0) {
            header <- data.table::fread(file, sep = delim, header = FALSE,
                                        stringsAsFactors = FALSE,
                                        data.table = FALSE, nrows = skip)
            names(loaded) <- header[format$colNames, ]
        } else {
            names(loaded) <- loaded[format$colNames, ]
        }
    }
    
    # Add row names
    if (!is.null(format$rowNames))
        rowNames <- unlist(loaded[, format$rowNames])
    else
        rowNames <- NULL
    
    # Filter out unwanted columns
    if (!is.null(format$ignoreCols)) loaded <- loaded[ , -format$ignoreCols]
    if (!is.null(format$ignoreRows)) {
        rowNames <- rowNames[-format$ignoreRows]
        loaded <- loaded[-format$ignoreRows, ]
    }
    
    # Add table name and description
    attr(loaded, "tablename") <- format$tablename
    attr(loaded, "description") <- format$description
    attr(loaded, "show") <- format$show
    
    # Add row names (it doesn't work placed before for some reason...)
    rownames(loaded) <- rowNames
    attr(loaded, "rowNames") <- !is.null(rowNames)
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
    
    # Get all information from available formats
    formats <- sourceScripts(formatsFolder, c("tablename", "check"))
    
    # The number of rows to read will be the maximum value asked by all the file
    # formats; if no format aks for a specific number of rows, the default is 6
    headRows <- lapply(formats, "[[", "header_rows")
    headRows <- unlist(rm.null(headRows))
    headRows <- ifelse(!is.null(headRows), max(headRows), 6)
    
    ## TODO(NunoA): check if fread makes this faster
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