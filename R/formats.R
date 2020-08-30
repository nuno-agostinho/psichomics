#' Checks the format of a file
#'
#' @details The name of the file may also be required to be considered of a
#' certain format.
#'
#' @param format Environment: format of the file
#' @param head Data.frame: head of the file to check
#' @param filename Character: name of the file
#'
#' @return \code{TRUE} if the file matches the given format's attributes
#' @keywords internal
checkFileFormat <- function(format, head, filename="") {
    # Replace with `browser()` to debug any format with the attribute `debug`
    if (!is.null(format$debug)) print(format$id)

    # If file name is of importance, check if the filename matches
    if (isTRUE(format$matchName) && !identical(filename, "") &&
        !grepl(format$filename, filename, fixed = TRUE))
        return(FALSE)

    ## TODO(NunoA): account for comments

    # Transpose data
    if (!is.null(format$transpose) && format$transpose)
        head <- t(head)

    lenCheck   <- length(format$check)
    checkByCol <- is.null(format$rowCheck) || !format$rowCheck
    if (checkByCol) {
        # Check for a match in desired column
        if (nrow(head) < lenCheck || ncol(head) < format$checkIndex) {
            return(FALSE)
        }
        desired <- head[seq(lenCheck), format$checkIndex]
    } else {
        # Check for a match in desired row
        if (ncol(head) < lenCheck || nrow(head) < format$checkIndex) {
            return(FALSE)
        }
        desired <- head[format$checkIndex, seq(lenCheck)]
    }
    if (all(is.na(desired))) return(FALSE)

    res <- all(trimws(desired) == format$check)
    if (!is.null(format$extraCheck)) res <- res && format$extraCheck(head)
    return(isTRUE(res))
}

addFileAttrs <- function(loaded, file, format) {
    # Add table name, description and other attributes
    rows   <- format$rows
    cols   <- format$columns
    loaded <- addObjectAttrs(loaded,
                             "filename"=file,
                             "dataType"=format$dataType,
                             "tablename"=format$tablename,
                             "description"=format$description,
                             "show"=format$show,
                             # Identify rows and columns
                             "rows"=ifelse(!is.null(rows), rows, "rows"),
                             "columns"=ifelse(!is.null(cols), cols, "columns"),
                             replace=FALSE)
    return(loaded)
}

#' Parse file according to its format
#'
#' @inheritParams checkFileFormat
#' @param file Character: file to load
#' @param ... Extra parameters passed to \link[data.table]{fread}
#' @param verbose Boolean: detail step while parsing?
#'
#' @details The resulting data frame includes the attribute \code{tablename}
#' with the name of the data frame
#'
#' @importFrom data.table fread
#' @importFrom stringr str_split_fixed
#'
#' @return Data frame with the loaded file
#' @keywords internal
parseFile <- function(format, file, ..., verbose=FALSE) {
    ## TODO(NunoA): account for the comment character
    delim            <- ifelse(!is.null(format$delim), format$delim, "\t")
    skip             <- ifelse(!is.null(format$skip), format$skip, 0)
    transpose        <- isTRUE(format$transpose)
    stringsAsFactors <- ifelse(!is.null(format$stringsAsFactors),
                               format$stringsAsFactors, TRUE)

    # Read file
    time <- Sys.time()
    if (verbose) message("Reading file...")
    loaded <- fread(file, sep=delim, header=FALSE, data.table=FALSE, skip=skip,
                    stringsAsFactors=!transpose && stringsAsFactors, ...)
    if (verbose) message("File read in ", format(Sys.time() - time))
    if (is.null(loaded)) {
        if (verbose) message("NULL returned while reading file")
        return(NULL)
    }

    time <- Sys.time()
    # Transpose data
    if (transpose) {
        if (verbose) message("Transposing data...")
        loaded <- data.frame(t(loaded), stringsAsFactors=stringsAsFactors,
                             row.names=NULL)
    }

    if (verbose) message("Processing data...")
    # Add column names from given row
    if (!is.null(format$colNames)) {
        if (skip != 0) {
            dots <- list(...)
            dots$nrows <- NULL

            header <- do.call(fread, c(list(
                input=file, sep=delim, header=FALSE, nrows=skip,
                stringsAsFactors=FALSE, data.table=FALSE), dots))
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
            dups <- duplicated(rowNames)
            if (sum(dups) > 0) {
                loaded <- loaded[!dups, ]
                warning(sprintf("Discarded %s row%s with duplicated rownames.",
                                sum(dups), if (sum(dups) > 1) "s" else ""))
            }
            rowNames <- as.character(loaded[ , format$rowNames])
        }
    } else {
        ## TODO(NunoA): Slow process... try to improve this
        if (verbose) message("Discarding duplicated rows...")
        if (!is.null(format$unique) && format$unique) loaded <- unique(loaded)
    }

    # Filter out unwanted columns
    if (!is.null(format$ignoreCols)) {
        loaded <- loaded[ , -format$ignoreCols, drop=FALSE]
    }
    if (!is.null(format$ignoreRows)) {
        rowNames <- rowNames[-format$ignoreRows]
        loaded <- loaded[-format$ignoreRows, ]
    }

    # Convert columns to numeric if data was transposed
    if (!is.null(format$transpose) && format$transpose) {
        for (col in seq(ncol(loaded))) {
            try <- tryCatch(as.numeric(as.character(loaded[ , col])),
                            warning=function(e) e)
            if (!"warning" %in% class(try))
                loaded[ , col] <- try
        }
    }

    # Add row names (it doesn't work placed before for some reason...)
    rownames(loaded) <- rowNames
    attr(loaded, "rowNames") <- !is.null(rowNames)

    # Further process the dataset if needed
    if (!is.null(format$process)) {
        if (verbose) message("Performing format-specific processing...")
        loaded <- format$process(loaded)
        if (is.null(loaded)) {
            if (verbose) message("Format-specific processing returned NULL")
            return(NULL)
        }
    }
    if (is.list(loaded) && !is.data.frame(loaded)) {
        loaded <- lapply(loaded, addFileAttrs, file, format)
    } else {
        loaded <- addFileAttrs(loaded, file, format)
    }
    if (verbose) message("Data processed in ", format(Sys.time() - time))
    return(loaded)
}

#' Load supported file formats
#'
#' @return Supported file formats
#' @keywords internal
loadFileFormats <- function() {
    # Get all functions ending with "UI"
    fun <- ls(getNamespace("psichomics"), all.names=TRUE, pattern="Format$")

    # Get the parameters of each function
    formats <- lapply(fun, function(format) {
        # Parse function name to get the function itself
        FUN <- eval(parse(text=format))
        # Check if module should be loaded by app
        if (loadBy("formats", FUN)) {
            # Remove last "UI" from the name and use it as the identifier
            item    <- FUN()
            item$id <- gsub("UI$", "", format)
            return(item)
        }
    })
    if (length(formats) == length(fun)) names(formats) <- fun
    # Remove NULL elements from list
    formats <- Filter(Negate(is.null), formats)
    return(formats)
}

#' Load file based on its format
#'
#' Tries to recognise the file format and parses the content of the given file
#' accordingly.
#'
#' @param file Character: file to parse
#' @param formats List of file formats to check
#' @param ... Extra parameters passed to \link[data.table]{fread}
#' @param verbose Boolean: detail steps while parsing
#' @param multiple Boolean: expect more than one file?
#'
#' @details The resulting data frame includes the attribute \code{tablename}
#' with the name of the data frame
#'
#' @importFrom utils read.delim
#' @importFrom data.table fread
#' @importFrom R.utils decompressFile
#'
#' @return Data frame with the contents of the given file if the file format is
#' recognised; otherwise, returns \code{NULL}
#' @keywords internal
loadFile <- function(file, formats=loadFileFormats(), ..., verbose=FALSE,
                     multiple=FALSE) {
    if (!is.list(formats[[1]])) formats <- list(formats)

    if (verbose) cat("\n")
    # The maximum number of rows to check a file is the maximum value asked by
    # the selected file formats; the default is 6
    headRows <- lapply(formats, "[[", "header_rows")
    headRows <- unlist(rm.null(headRows))
    headRows <- ifelse(!is.null(headRows), max(headRows), 6)

    # If file is compressed, decompress first (based on data.table::fread)
    ext <- file_ext(file)
    if (ext %in% c("gz", "bz2")) {
        time <- Sys.time()
        if (verbose) message("Decompressing to temporary file...")
        FUN <- switch(ext, "gz"=gzfile, "bz2"=bzfile)
        decompressFile(file, decompFile <- tempfile(tmpdir=tempdir()), ext=NULL,
                       FUN=FUN, remove=FALSE)
        file <- decompFile
        on.exit(unlink(decompFile), add=TRUE)
        if (verbose) message("Decompressed in ", format(Sys.time() - time))
    }
    head <- fread(file, header=FALSE, nrows=headRows, stringsAsFactors=FALSE,
                  data.table=FALSE)

    # Check if the file is recognised by at least one file format
    if (verbose) message("Checking file format...")
    recognised <- lapply(formats, checkFileFormat, head, file)
    recognised <- unlist(recognised)

    if (verbose) {
        message("Recognised format: ",
                paste(names(recognised[recognised]), collapse=", "))
    }

    if (sum(recognised) > 1) {
        ## TODO(NunoA): ask the user which file format to use
        stop("Error: more than one file format was recognised.")
    } else if (sum(recognised) == 1) {
        format <- formats[recognised][[1]]
        loaded <- parseFile(format, file, ..., verbose=verbose)
        # Avoid returning more than one dataset when not expected
        if (!multiple && !is.data.frame(loaded) && is.list(loaded)) {
            loaded <- loaded[[1]]
        }
        return(loaded)
    }
}
