#' Checks the format of a file
#'
#' @details The name of the file may also be required to be considered of a 
#' certain format.
#'
#' @param format Environment: format of the file
#' @param head Data.frame: head of the file to check
#' @param filename Character: name of the file
#'
#' @return TRUE if the file is of the given format; otherwise, returns FALSE
checkFileFormat <- function(format, head, filename) {
    # If file name is of importance, check if the filename matches
    if (isTRUE(format$matchName) &&
        !grepl(format$filename, filename, fixed = TRUE))
        return(FALSE)
    
    ## TODO(NunoA): account for comments
    
    # Transpose data
    if (!is.null(format$transpose) && format$transpose)
        head <- t(head)
    
    lenCheck <- length(format$check)
    # Check if using row or column when checking format
    if (is.null(format$rowCheck) || !format$rowCheck) {
        if (nrow(head) < lenCheck) return(FALSE)
        # Select wanted row and check for a match
        desired <- trimws(head[1:lenCheck, format$checkIndex])
        return(all(desired == format$check))
    } else {
        if (ncol(head) < lenCheck) return(FALSE)
        # Select wanted column and check for a match
        desired <- trimws(head[format$checkIndex, 1:lenCheck])
        return(all(desired == format$check))
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
#' @importFrom data.table fread
#' @importFrom stringr str_split_fixed
#' 
#' @return Data frame with the loaded file
loadFile <- function(format, file) {
    ## TODO(NunoA): account for the comment character
    delim <- ifelse(!is.null(format$delim), format$delim, "\t")
    skip <- ifelse(!is.null(format$skip), format$skip, 0)
    
    transpose <- !is.null(format$transpose) && format$transpose
    loaded <- fread(file, sep = delim, header = FALSE,
                    stringsAsFactors = !transpose, data.table = FALSE,
                    skip = skip)
    
    # Transpose data
    if (transpose) {
        loaded <- data.frame(t(loaded), stringsAsFactors = TRUE, 
                             row.names = NULL)
    }
    
    # Add column names from given row
    if (!is.null(format$colNames)) {
        if (skip != 0) {
            header <- fread(file, sep = delim, header = FALSE, nrows = skip,
                            stringsAsFactors = FALSE, data.table = FALSE)
            names(loaded) <- header[format$colNames, ]
        } else {
            names(loaded) <- unname(vapply(loaded[format$colNames, ],
                                           as.character, character(1)))
            
        }
    }
    
    # Add row names and remove duplicated rows
    rowNames <- NULL
    if (!is.null(format$rowNames)) { 
        rowNames <- as.character(loaded[, format$rowNames])
        if (!is.null(format$unique) && format$unique) {
            loaded <- loaded[!duplicated(rowNames), ]
            rowNames <- as.character(loaded[, format$rowNames])
        }
    } else {
        ## TODO(NunoA): Slow process... try to improve this
        if (!is.null(format$unique) && format$unique) loaded <- unique(loaded)
    }
    
    # Filter out unwanted columns
    if (!is.null(format$ignoreCols)) loaded <- loaded[ , -format$ignoreCols]
    if (!is.null(format$ignoreRows)) {
        rowNames <- rowNames[-format$ignoreRows]
        loaded <- loaded[-format$ignoreRows, ]
    }
    
    # Convert columns to numeric if data was transposed
    if (!is.null(format$transpose) && format$transpose) {
        for (col in 1:ncol(loaded)) {
            try <- tryCatch(as.numeric(as.character(loaded[ , col])),
                            warning=function(e) e)
            if (!"warning" %in% class(try))
                loaded[ , col] <- try
        }
    }
    
    # Add row names (it doesn't work placed before for some reason...)
    rownames(loaded) <- rowNames
    attr(loaded, "rowNames") <- !is.null(rowNames)
    attr(loaded, "filename") <- file
    attr(loaded, "dataType") <- format$dataType
    
    # Further process the dataset if needed
    if (!is.null(format$process)) loaded <- format$process(loaded)
    
    # Add table name and description
    attr(loaded, "tablename") <- format$tablename
    attr(loaded, "description") <- format$description
    attr(loaded, "show") <- format$show
    return(loaded)
}

#' Loads file formats
#' @return Loaded file formats available
loadFileFormats <- function() {
    # Get all functions ending with "UI"
    fun <- ls(getNamespace("psichomics"), all.names=TRUE, pattern="Format$")
    
    # Get the parameters of each function
    formats <- lapply(fun, function(format) {
        # Parse function name to get the function itself
        FUN <- eval(parse(text=format))
        # Check if module should be loaded by app
        if (loadBy("formats", FUN)) {
            # Remove last "UI" from the name and use it as ID
            name <- gsub("UI$", "", format)
            FUN()
        }
    })
    # Remove NULL elements from list
    formats <- Filter(Negate(is.null), formats)
    return(formats)
}

#' Parse file given a list of file formats
#' 
#' Tries to recognise the file format and parses the content of the given file
#' accordingly.
#' 
#' @param file Character: file to parse
#' @param formats List of file formats to check
#' 
#' @details The resulting data frame includes the attribute "tablename" with the
#' name of the data frame
#' 
#' @importFrom utils read.delim
#' 
#' @return Data frame with the contents of the given file if the file format is
#' recognised; otherwise, returns NULL
parseValidFile <- function(file, formats) {
    # The number of rows to read will be the maximum value asked by all the file
    # formats; if no format aks for a specific number of rows, the default is 6
    headRows <- lapply(formats, "[[", "header_rows")
    headRows <- unlist(rm.null(headRows))
    headRows <- ifelse(!is.null(headRows), max(headRows), 6)
    
    ## TODO(NunoA): check if fread makes this faster
    head <- read.delim(file, header=FALSE, nrows=6, stringsAsFactors=FALSE)
    
    # Check if the file is recognised by at least one file format
    recognised <- lapply(formats, checkFileFormat, head, file)
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