# Get FireBrowse parameters ----------------------------------------------------

#' Get available parameters for TCGA data
#'
#' Parameters obtained via \href{http://firebrowse.org/api-docs/}{FireBrowse}
#'
#' @family functions associated with TCGA data retrieval
#' @return Parsed response
#'
#' @importFrom R.utils capitalize
#'
#' @aliases getFirebrowseDataTypes
#' @export
#'
#' @examples
#' getTCGAdataTypes()
getTCGAdataTypes <- function() {
    choices <- list("RNA sequencing"=c(
        "junction_quantification", "exon_quantification",
        "exon_expression", "junction_expression",
        "RSEM_genes", "RSEM_genes_normalized", "RSEM_isoforms", "Preprocess"))
    names(choices[[1]]) <- capitalize(gsub("_", " ", choices[[1]], fixed=TRUE))
    return(choices)
}

#' @export
getFirebrowseDataTypes <- getTCGAdataTypes

#' @rdname getTCGAdataTypes
#' @aliases getFirebrowseDates
#' @export
#'
#' @examples
#' if (isFirebrowseUp()) getTCGAdates()
getTCGAdates <- function() {
    dates <- parseFirebrowseMetadata("Dates")$Dates
    format <- getFirebrowseDateFormat()
    dates <- as.Date(dates, format$query)
    return(dates)
}

#' @export
getFirebrowseDates <- getTCGAdates

#' @rdname getTCGAdataTypes
#'
#' @param cohort Character: filter results by cohorts (optional)
#'
#' @aliases getFirebrowseCohorts
#' @export
#'
#' @examples
#' if (isFirebrowseUp()) getTCGAcohorts()
getTCGAcohorts <- function(cohort = NULL) {
    response <- parseFirebrowseMetadata("Cohorts", cohort=cohort)
    cohorts <- response$Cohorts[[2]]
    names(cohorts) <- response$Cohorts[[1]]
    return(cohorts)
}

#' @export
getFirebrowseCohorts <- getTCGAcohorts


# Process and manipulate FireBrowse data ---------------------------------------

#' @rdname parseTCGAsampleInfo
#'
#' @param filename Character: path to RDS file containing corresponding types
#'
#' @aliases parseSampleGroups
#' @export
#'
#' @examples
#' parseTCGAsampleTypes(c("TCGA-01A-Tumour", "TCGA-10B-Normal"))
parseTCGAsampleTypes <- function(samples, filename = system.file(
    "extdata", "TCGAsampleType.RDS", package="psichomics")) {
    typeList <- readRDS(filename)
    type <- gsub(".*?-([0-9]{2}).-.*", "\\1", samples, perl = TRUE)
    return(typeList[type])
}

#' @export
parseSampleGroups <- parseTCGAsampleTypes

#' Parse sample information from TCGA sample identifiers
#'
#' @param samples Character: sample identifiers
#' @param match Integer: match between samples and subjects (\code{NULL} by
#' default; performs the match)
#'
#' @aliases parseTcgaSampleInfo
#' @family functions associated with TCGA data retrieval
#'
#' @return Metadata associated with each TCGA sample
#' @export
#'
#' @examples
#' samples <- c("TCGA-3C-AAAU-01A-11R-A41B-07", "TCGA-3C-AALI-01A-11R-A41B-07",
#'              "TCGA-3C-AALJ-01A-31R-A41B-07", "TCGA-3C-AALK-01A-11R-A41B-07",
#'              "TCGA-4H-AAAK-01A-12R-A41B-07", "TCGA-5L-AAT0-01A-12R-A41B-07")
#'
#' parseTCGAsampleInfo(samples)
parseTCGAsampleInfo <- function(samples, match=NULL) {
    parsed <- parseTCGAsampleTypes(samples)
    if ( all(is.na(parsed)) ) return(NULL)

    info <- data.frame(parsed)
    colnames(info) <- "Sample types"
    rownames(info) <- samples

    if (is.null(match)) match <- getSubjectFromSample(samples)
    info <- cbind(info, "Patient ID"=match)

    # Metadata
    attr(info, "rowNames")    <- TRUE
    attr(info, "description") <- "Metadata for TCGA samples"
    attr(info, "dataType")    <- "Sample metadata"
    attr(info, "tablename")   <- "Sample metadata"
    attr(info, "rows")        <- "samples"
    attr(info, "columns")     <- "attributes"
    return(info)
}

#' @export
parseTcgaSampleInfo <- parseTCGAsampleInfo

#' Retrieve URLs from a response to a FireBrowse data query
#'
#' @param res Response from \code{httr::GET} to a FireBrowse data query
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content
#'
#' @return Named character with URLs
#' @keywords internal
#'
#' @examples
#' res <- psichomics:::queryFirebrowseData(cohort = "ACC")
#' url <- psichomics:::parseUrlsFromFirebrowseResponse(res)
parseUrlsFromFirebrowseResponse <- function(res) {
    # Parse the query response
    parsed <- content(res, "text", encoding = "UTF8")
    parsed <- fromJSON(parsed)[[1]]
    parsed$date <- parseDateResponse(parsed$date)

    # Get cohort names
    cohort <- getTCGAcohorts()
    cohort <- cohort[parsed$cohort]

    ## TODO(NunoA): maybe this could be simplified?
    # Split URLs from response by cohort and datestamp
    url <- split(parsed$url, paste(cohort, format(parsed$date, "%Y-%m-%d")))
    url <- lapply(url, unlist)
    link <- unlist(url)
    names(link) <- rep(names(url), vapply(url, length, numeric(1)))
    return(link)
}

#' Returns the date format used by the FireBrowse API
#'
#' @return Named list with date formats from FireBrowse API
#' @keywords internal
#'
#' @examples
#' format <- psichomics:::getFirebrowseDateFormat()
#'
#' # date format to use in a query to FireBrowse API
#' format$query
#'
#' # date format to parse a date in a response from FireBrowse API
#' format$response
getFirebrowseDateFormat <- function() {
    query <- "%Y_%m_%d"
    response <- "%d %b %Y"
    return(list(query=query, response=response))
}

#' Parse the date from a response
#'
#' @param string Character: dates
#'
#' @return Parsed date
#' @keywords internal
parseDateResponse <- function(string) {
    format <- getFirebrowseDateFormat()$response
    date <- strsplit(string, " ")
    date <- lapply(date, function(i) paste(i[2:4], collapse=" "))
    date <- as.Date(unlist(date), format=format)
    return(date)
}

#' Query the FireBrowse API for metadata
#'
#' @param type Character: metadata to retrieve
#' @param ... Character: parameters to pass to query (optional)
#'
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#'
#' @return List with parsed response
#' @keywords internal
#'
#' @examples
#' psichomics:::parseFirebrowseMetadata("Dates")
#' psichomics:::parseFirebrowseMetadata("Centers")
#' psichomics:::parseFirebrowseMetadata("HeartBeat")
#'
#' # Get the abbreviation and description of all cohorts available
#' psichomics:::parseFirebrowseMetadata("Cohorts")
#' # Get the abbreviation and description of the selected cohorts
#' psichomics:::parseFirebrowseMetadata("Cohorts", cohort = c("ACC", "BRCA"))
parseFirebrowseMetadata <- function(type, ...) {
    # Remove NULL arguments
    args <- Filter(Negate(is.null), list(...))

    # Query multiple items by collapsing them with a comma
    if (length(args) > 0)
        args <- lapply(args, paste, collapse = ",")

    # Query the given link and parse the response
    link <- paste0("http://firebrowse.org/api/v1/Metadata/", type)
    response <- GET(link, query = c(format = "json", args))
    response <- fromJSON(content(response, "text", encoding = "UTF8"))
    return(response)
}

# Download and load FireBrowse data --------------------------------------------

#' Check if \href{http://firebrowse.org/api-docs/}{FireBrowse API} is running
#'
#' @importFrom httr GET warn_for_status http_error
#' @importFrom methods is
#'
#' @family functions associated with TCGA data retrieval
#' @return Invisible \code{TRUE} if the
#' \href{http://firebrowse.org/api-docs/}{FireBrowse API} is working; otherwise,
#' raises a warning with the status code and a brief explanation.
#' @export
#'
#' @examples
#' isFirebrowseUp()
isFirebrowseUp <- function() {
    link <- paste0("http://firebrowse.org/api/v1/Metadata/HeartBeat")
    heartbeat <- tryCatch(GET(link, query=list(format="json")), error=return)
    if (is(heartbeat, "error")) {
        return(FALSE)
    } else if (http_error(heartbeat)) {
        warn_for_status(heartbeat, "reach FireBrowse API")
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' Prepare TCGA sample metadata from loaded datasets
#'
#' If no TCGA datasets apply, the input is returned
#'
#' @param data List of list of data frames
#'
#' @return List of list of data frames
#' @keywords internal
loadTCGAsampleMetadata <- function(data) {
    for (i in seq(data)) {
        # Retrieve sample metadata from junction quantification
        match <- sapply(data[[i]], attr, "dataType") ==
            "Junction quantification"
        junctionQuantSamples <- NULL
        if (any(match)) {
            junctionQuant <- data[[i]][match]
            if (!is.null(junctionQuant)) {
                samples <- unique(unlist(lapply(junctionQuant, colnames)))
                if (any(grepl("^TCGA", samples))) {
                    junctionQuantSamples <- samples
                    data[[i]]$"Sample metadata" <- parseTCGAsampleInfo(samples)
                }
            }
        }

        # Retrieve sample metadata from gene expression
        match <- sapply(data[[i]], attr, "dataType") == "Gene expression"
        if (any(match)) {
            geneExpr <- data[[i]][match]
            if (!is.null(geneExpr)) {
                samples <- unique(unlist(lapply(geneExpr, colnames)))
                samples <- samples[!samples %in% junctionQuantSamples]
                if (any(grepl("^TCGA", samples))) {
                    data[[i]]$"Sample metadata" <- parseTCGAsampleInfo(samples)
                }
            }
        }
    }
    return(data)
}

#' Query the FireBrowse API for TCGA data
#'
#' @param format Character: response format as \code{JSON}, \code{CSV} or
#' \code{TSV}
#' @param date Character: dates of the data retrieval by FireBrowse (by default,
#' it uses the most recent data available)
#' @param cohort Character: abbreviation of the cohorts (by default, returns
#' data for all cohorts)
#' @param data_type Character: data types (optional)
#' @param tool Character: data produced by the selected FireBrowse tools
#' (optional)
#' @param platform Character: data generation platforms (optional)
#' @param center Character: data generation centres (optional)
#' @param level Integer: data levels (optional)
#' @param protocol Character: sample characterization protocols (optional)
#' @param page Integer: page of the results to return (optional)
#' @param page_size Integer: number of records per page of results (optional)
#' @param sort_by String: column used to sort the data (by default, sort by
#' cohort)
#'
#' @importFrom httr GET
#'
#' @return Response from the FireBrowse API (it needs to be parsed)
#' @keywords internal
#'
#' @examples
#' cohort <- getTCGAcohorts()[1]
#' psichomics:::queryFirebrowseData(cohort = names(cohort),
#'                                  data_type = "mRNASeq")
#'
#' # Querying for data from a specific date
#' dates <- getTCGAdates()
#' dates <- format(dates, psichomics:::getFirebrowseDateFormat()$query)
#'
#' psichomics:::queryFirebrowseData(date = dates[2], cohort = names(cohort))
queryFirebrowseData <- function(format = "json", date = NULL, cohort = NULL,
                                data_type = NULL, tool = NULL, platform = NULL,
                                center = NULL, level = NULL, protocol = NULL,
                                page = NULL, page_size = NULL, sort_by = NULL) {
    # Only allow these response formats
    format <- match.arg(format, c("json", "csv", "tsv"))

    # Use most recent date by default
    if (is.null(date)) date <- getTCGAdates()[[1]]
    date <- gsub("-", "_", date, fixed=TRUE)

    # Process the parameters of the query
    labels <- list("format", "date", "cohort", "data_type", "tool", "platform",
                   "center", "level", "protocol", "page", "page_size",
                   "sort_by")
    query <- lapply(labels, dynGet)
    names(query) <- labels
    query <- Filter(Negate(is.null), query)

    # Collapse items with a comma to query for multiple items
    query <- lapply(query, paste, collapse = ",")

    # Query the API
    response <- GET("http://firebrowse.org", query = query,
                    path = "api/v1/Archives/StandardData")
    return(response)
}

#' Download files to a given directory
#'
#' @param url Character: download links
#' @param folder Character: directory to store the downloaded archives
#' @param ... Extra parameters passed to the download function
#' @param download Function to use to download files
#'
#' @importFrom utils download.file
#'
#' @return Invisible TRUE if every file was successfully downloaded
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' url <- paste0("https://unsplash.it/400/300/?image=", 570:572)
#' psichomics:::downloadFiles(url, "~/Pictures")
#'
#' # Download without printing to console
#' psichomics:::downloadFiles(url, "~/Pictures", quiet = TRUE)
#' }
downloadFiles <- function(url, folder, download = download.file, ...) {
    destination <- file.path(folder, basename(url))
    for (i in seq_along(url)) {
        updateProgress("Downloading file", detail = basename(url[i]), i,
                       length(url))
        download(url[i], destination[i], ...)
    }
    display("Downloading completed")
    return(destination)
}

#' Compute the 32-byte \code{MD5} hashes of one or more files and check with
#' given \code{md5} file
#'
#' @param filesToCheck Character: files to calculate and match \code{MD5} hashes
#' @param md5file Character: file containing correct \code{MD5} hashes
#'
#' @importFrom digest digest
#' @importFrom utils read.table
#'
#' @return Logical vector showing TRUE for files with matching \code{md5sums}
#' and \code{FALSE} for files with non-matching \code{md5sums}
#' @keywords internal
checkIntegrity <- function(filesToCheck, md5file) {
    if (is.na(md5file)) return(FALSE)
    md5sums <- digest(file = filesToCheck)
    md5table <- read.table(md5file, stringsAsFactors = FALSE)[[1]]
    return(md5sums %in% md5table)
}

#' Prepares FireBrowse archives in a given directory
#'
#' Checks FireBrowse archives' integrity using the MD5 files, extracts the
#' content of the archives, moves the content to newly-created folders and
#' removes the original downloaded archives.
#'
#' @param archive Character: path to downloaded archives
#' @param md5 Character: path to MD5 files of each archive
#' @param folder Character: master directory where every archive will be
#' extracted
#' @param outdir Character: subdirectories where to move the extracted content
#'
#' @importFrom utils untar
#'
#' @return Invisible TRUE if successful
#' @keywords internal
#'
#' @examples
#' file <- paste0(
#'     "~/Downloads",
#'     "ACC/20151101/gdac.broadinstitute.org_ACC.",
#'     "Merge_Clinical.Level_1.2015110100.0.0.tar.gz")
#' md5 <- paste0(file, ".md5")
#' \dontrun{
#' prepareFirebrowseArchives(archive = file, md5 = paste0(file, ".md5"))
#' }
prepareFirebrowseArchives <- function(archive, md5, folder, outdir) {
    # Check integrety of the downloaded archives with the MD5 files
    validFiles <- simplify2array(Map(checkIntegrity, archive, md5))

    # Extract the contents of the archives to the same folder
    invisible(lapply(archive, function(arc) untar(arc, exdir=dirname(arc))))

    # Remove the original downloaded files
    invisible(file.remove(archive, md5))

    # Create folders for these data
    directory <- file.path(folder, names(outdir))
    existing <- dir.exists(directory)
    lapply(directory[!existing], dir.create)

    # Move the files to the newly-created folders
    arc <- file_path_sans_ext(archive, compression=TRUE)
    basen <- basename(arc)

    ns <- sapply(basen, function(k) {
        m <- sapply(outdir, function(p) basename(p) == k)
        names(outdir)[col(m)[m]]
    })

    tmpNames <- file.path(dirname(arc),
                          paste0("PSIchomics_TMP", seq(length(arc))))
    # Rename folder name first and revert it afterwards to avoid reaching
    # the Windows' character limit (260 for paths)
    file.rename(arc, tmpNames)

    # Rename file names
    filenames <- list.files(tmpNames, full.names = TRUE)
    basenames <- basename(filenames)
    renamed <- gsub(".*?_.*?\\.(.*)\\.Level.*", "\\1", basenames)
    renamed <- gsub("rnaseqv2__|unc_edu__Level_.__", "", renamed)
    toRename <- filenames[basenames != renamed]
    renamed <- renamed[basenames != renamed]
    file.rename(toRename, file.path(dirname(toRename), renamed))
    file.rename(tmpNames, arc)

    # Organise downloaded folders
    file.rename(arc, file.path(dirname(arc), ns, basen))
    return(invisible(TRUE))
}

#' Load FireBrowse folders
#'
#' Loads the files present in each folder as a data.frame.
#'
#' @note For faster execution, this function uses the \code{readr} library. This
#' function ignores subfolders of the given folder (which means that files
#' inside subfolders are NOT loaded).
#'
#' @include formats.R
#'
#' @param folder Character: folder(s) in which to look for FireBrowse files
#' @param exclude Character: files to exclude from the loading
#'
#' @return List with loaded data.frames
#' @keywords internal
loadFirebrowseFolders <- function(folder, exclude="") {
    # Retrieve full path of the files inside the given folders
    files <- dir(folder, full.names=TRUE)

    # Exclude subdirectories and undesired files
    files <- files[!dir.exists(files)]
    exclude <- paste(exclude, collapse = "|")
    if (exclude != "") files <- files[!grepl(exclude, files)]

    # Try to load files and remove those with 0 rows
    loaded <- list()
    formats <- loadFileFormats()
    for (each in seq_along(files)) {
        updateProgress("Processing file", detail = basename(files[each]), each,
                       length(files))
        loaded[[each]] <- loadFile(files[each], formats)
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    return(loaded)
}

#' Download and process TCGA data
#'
#' TCGA data obtained via \href{http://firebrowse.org/api-docs/}{FireBrowse}
#'
#' @param folder Character: directory to store the downloaded archives (by
#' default, saves to \code{\link{getDownloadsFolder}()})
#' @param data Character: data to load (see \code{\link{getTCGAdataTypes}()})
#' @param exclude Character: files and folders to exclude from downloading and
#' from loading into R (by default, exclude files containing \code{.aux.},
#' \code{.mage-tab.} and \code{MANIFEST.TXT})
#' @inheritDotParams queryFirebrowseData -format
#' @param download Boolean: download missing files
#'
#' @include formats.R
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom httr stop_for_status
#' @importFrom utils download.file
#'
#' @aliases loadFirebrowseData
#' @family functions associated with TCGA data retrieval
#' @family functions to load data
#' @return A list with the loaded data, unless required files are unavailable
#' and \code{download = FALSE} (if so, it returns the URL of files to download)
#' @export
#'
#' @examples
#' getTCGAcohorts()
#' getTCGAdataTypes()
#' \dontrun{
#' loadTCGAdata(cohort = "ACC", data_type = "Clinical")
#' }
loadTCGAdata <- function(folder=getDownloadsFolder(),
                         data=c("clinical", "junction_quantification",
                                "RSEM_genes"),
                         exclude=c(".aux.", ".mage-tab.", "MANIFEST.txt"),
                         ..., download=TRUE) {
    args <- list(...)

    datasets <- unlist(getTCGAdataTypes())
    # Data types to load
    args$data_type <- c(data[!data %in% datasets], "mRNASeq")
    # Datasets to ignore
    exclude <- c(exclude, datasets[!datasets %in% data])

    # Ask for maximum number of records
    args$page_size <- 2000

    # Query FireBrowse and get URLs for archives
    res <- do.call(queryFirebrowseData, args)
    stop_for_status(res)
    url <- parseUrlsFromFirebrowseResponse(res)

    # "RSEM_genes" would match both "RSEM_genes" and "RSEM_genes_normalized"
    exclude[exclude=="RSEM_genes"] <- "RSEM_genes__"

    # Do not download specific items
    exclude <- paste(escape(exclude), collapse = "|")
    url <- url[!grepl(exclude, url)]

    # Get the file names without extensions
    md5  <- file_ext(url) == "md5"
    base <- basename(url)
    base[!md5] <- file_path_sans_ext(base[!md5], compression = TRUE)

    # Check which files are missing from the given directory
    downloadedFiles <- list.files(folder, recursive=TRUE, full.names=TRUE,
                                  include.dirs=TRUE)
    downloadedMD5   <- file_ext(downloadedFiles) == "md5"
    fullPath <- function(files)
        downloadedFiles[match(files, basename(downloadedFiles))]

    # Downloaded files may be archived
    possibleExtensions <- lapply(base[!md5], paste0, c("", ".tar", ".tar.gz"))

    missing <- logical(length(base))
    missing[md5]  <- !base[md5] %in% basename(downloadedFiles[downloadedMD5])
    missing[!md5] <- vapply(possibleExtensions, function (i)
        !any(i %in% basename(downloadedFiles[!downloadedMD5])),
        FUN.VALUE = logical(1))

    if (any(missing[!md5])) {
        missingFiles <- url[missing]
        class(missingFiles) <- c("missing", class(missingFiles))

        if (download) {
            # Download missing files
            updateProgress(divisions = 1)
            display("Downloading files...")

            if (identical(getOption("download.file.method"), "libcurl")) {
                dl <- download.file(missingFiles, destfile=file.path(
                    folder, basename(missingFiles)))
            } else {
                for (k in seq_along(missingFiles)) {
                    f <- missingFiles[k]
                    dl <- download.file(
                        f, destfile=file.path(folder, basename(f)))
                }
            }

            # Check again the files in the given directory
            downloadedFiles <- list.files(folder, recursive=TRUE,
                                          full.names=TRUE, include.dirs=TRUE)
            downloadedMD5   <- file_ext(downloadedFiles) == "md5"
        } else {
            return(missingFiles)
        }
    }

    # Split folders by the cohort type and date
    categories <- names(url[!md5])
    folders <- base[!md5]
    folders <- split(folders, categories)

    # Check for folders to unarchive
    archives <- unlist(lapply(possibleExtensions, function (item)
        item[item %in% basename(downloadedFiles[!downloadedMD5])]))
    tar <- grepl(".tar", archives, fixed = TRUE)

    if (length(archives[tar]) > 0) {
        # Extract the content, check intergrity, move archives to newly-created
        # folders and remove original archives
        updateProgress("Extracting archives...", divisions=1 + length(folders))
        prepareFirebrowseArchives(fullPath(archives[tar]),
                                  fullPath(base[md5][tar]), folder, folders)
        updateProgress("Archives prepared")
    } else {
        # Set the progress bar to the number of folders to load
        updateProgress("Loading data...", divisions = length(folders))
    }

    # Get the full path of the files
    downloadedFiles <- list.files(folder, recursive=TRUE, full.names=TRUE,
                                  include.dirs=TRUE)
    folders <- fullPath(base[!md5])
    folders <- split(folders, categories)

    # Load the files (but discard empty files with no rows)
    loaded <- lapply(folders, loadFirebrowseFolders, exclude)
    loaded <- lapply(loaded, function(i) Filter(nrow, i))

    loaded <- loadTCGAsampleMetadata(loaded)
    loaded <- processDatasetNames(loaded)
    return(loaded)
}

#' @export
loadFirebrowseData <- loadTCGAdata

# Shiny-specific FireBrowse functions ------------------------------------------

#' Creates a UI set with options to add data from TCGA/FireBrowse
#' @param ns Namespace function
#'
#' @importFrom shiny tagList uiOutput selectizeInput actionButton textAreaInput
#'
#' @return A UI set that can be added to a UI definition
#' @keywords internal
addTCGAdata <- function(ns) {
    cohorts <- getTCGAcohorts()
    acronyms <- names(cohorts)
    names(acronyms) <- sprintf("%s (%s)", cohorts, names(cohorts))

    dates <- as.character(getTCGAdates())
    names(dates) <- dates
    names(dates)[1] <- paste(names(dates)[1], "(most recent)")

    dataTypes <- getTCGAdataTypes()
    dataList <- dataTypes[[1]]
    names(dataList)[dataList == "RSEM_genes"] <- "Gene expression (RSEM)"
    names(dataList)[dataList == "RSEM_genes_normalized"] <-
        "Gene expression (normalised by RSEM)"
    dataList <- dataList[dataList %in% c("junction_quantification",
                                         "RSEM_genes", "RSEM_genes_normalized")]
    dataTypes[[1]] <- dataList
    dataTypes <- c("Clinical data"="Clinical", dataTypes)

    tagList(
        uiOutput(ns("firebrowseDataModal")),
        uiOutput(ns("iframeDownload")),
        selectizeInput(ns("firebrowseCohort"), "Tumour type", acronyms,
                       width = "100%", multiple = TRUE, options = list(
                           placeholder = "Select cohort(s)",
                           plugins=list("remove_button"))),
        selectizeInput(ns("firebrowseDate"), "Date", dates, multiple = TRUE,
                       width = "100%", selected = dates[1], options = list(
                           placeholder = "Select sample date",
                           plugins=list("remove_button"))),
        selectizeInput(ns("firebrowseData"), "Data type", multiple = TRUE,
                       width = "100%", dataTypes,
                       selected=c("Clinical", "junction_quantification",
                                  "RSEM_genes"),
                       options = list(
                           placeholder = "Select data types",
                           plugins=list("remove_button"))),
        browseDownloadFolderInput(ns("dataFolder")),
        processButton(ns("getFirebrowseData"), "Load data"))
}

#' @rdname appUI
#' @param panel Function to enclose interface
#'
#' @importFrom shiny NS helpText icon a
firebrowseUI <- function(id, panel) {
    ns <- NS(id)

    panel(style="info",
          title=list(icon("plus-circle"), "TCGA data loading"),
          value="Load TCGA/FireBrowse data",
          helpText("TCGA data are downloaded using the",
                   a(href="http://firebrowse.org", target="_blank",
                     "FireBrowse"), "API."),
          div(id=ns("firebrowseLoading"), class="progress",
              div(class="progress-bar progress-bar-striped active",
                  role="progressbar", style="width: 100%", "Loading...")),
          uiOutput(ns("checkFirebrowse")))
}

#' Return an user interface depending on the status of the FireBrowse API
#'
#' If the API is working, it'll be loaded. Else, a message will appear warning
#' the user that the API is down and that will let check again if the API is
#' back online.
#'
#' @param ns Namespace function
#'
#' @importFrom shiny br icon tagList actionButton
#' @importFrom shinyjs hide
#'
#' @return HTML elements
#' @keywords internal
checkFirebrowse <- function(ns) {
    # startProgress("Checking FireBrowse API to retrieve TCGA data...", 1)
    if (isFirebrowseUp()) {
        # updateProgress("Loading FireBrowse interface...")
        ui <- addTCGAdata(ns)
    } else {
        ui <- errorDialog("FireBrowse API appears to be offline at the moment.",
                          buttonId=ns("refreshFireBrowse"),
                          buttonLabel="Check FireBrowse again",
                          buttonIcon="refresh")
    }
    # closeProgress("FireBrowse interface loaded")
    hide("firebrowseLoading")
    return(ui)
}

#' Set data from FireBrowse
#'
#' @inheritParams appServer
#' @param replace Boolean: replace loaded data?
#'
#' @importFrom shinyjs disable enable
#' @importFrom shiny div fluidRow column icon tags
#' @importFrom shinyBS bsTooltip
#'
#' @inherit psichomics return
#' @keywords internal
setFirebrowseData <- function(input, output, session, replace=TRUE) {
    ns <- session$ns
    time <- startProcess("getFirebrowseData")

    # Load data from FireBrowse
    data <- loadFirebrowseData(folder = input$dataFolder,
                               cohort = input$firebrowseCohort,
                               date = gsub("-", "_", input$firebrowseDate),
                               data = input$firebrowseData, download = FALSE)

    areDataMissing <- any(is(data, "missing"))
    if (areDataMissing) {
        updateProgress(divisions = 1)
        setURLtoDownload(data)

        infoModal(
            session, "Confirm data download",
            "Do you wish to download the selected TCGA data? When the",
            "downloads finish, click", tags$b("Load data"), "with the exact",
            "same options to automatically load the data into",
            tags$i("psichomics."), tags$br(), tags$br(), tags$div(
                class="alert", class="alert-warning", role="alert",
                tags$i("psichomics"), "will check for downloaded files in",
                tags$br(), tags$kbd(prepareWordBreak(input$dataFolder)),
                if (isRStudioServer())
                    tagList(
                        tags$br(), tags$br(), icon("exclamation-circle"),
                        "Running", tags$i("psichomics"), "in a remote server?",
                        "Please make sure to move the files from your computer",
                        "to the aforementioned server's folder.")),
            modalId="firebrowseDataModal", caller="Load TCGA data",
            footer=actionButton(ns("acceptDownload"), "Download data",
                                class="btn-primary", "data-dismiss"="modal"))
    } else if (!is.null(data)) {
        if(replace) {
            setData(data)
        } else {
            data <- processDatasetNames(c(getData(), data))
            setData(data)
        }
    }
    endProcess("getFirebrowseData", if (!areDataMissing) time)
}

#' @rdname appServer
firebrowseServer <- function(input, output, session) {
    ns <- session$ns

    prepareFileBrowser(session, input, "dataFolder", directory=TRUE)

    # If FireBrowse is unaccessible, allow user to try again
    output$checkFirebrowse <- renderUI(isolate(checkFirebrowse(ns)))
    observeEvent(input$refreshFirebrowse,
                 output$checkFirebrowse <- renderUI(checkFirebrowse(ns)))

    # # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))

    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input$getFirebrowseData, {
        if (length(isolate(input$firebrowseCohort)) == 0) {
            errorModal(session, "No tumour type selected",
                       "Please, input a tumour type.",
                       caller="Load TCGA data", modalId="firebrowseDataModal")
        } else if (length(isolate(input$firebrowseDate)) == 0) {
            errorModal(session, "No date selected",
                       "Please, input date of samples of interest.",
                       caller="Load TCGA data", modalId="firebrowseDataModal")
        } else if (length(isolate(input$firebrowseData)) == 0) {
            errorModal(session, "No data types select",
                       "Please, input data types of interest.",
                       caller="Load TCGA data", modalId="firebrowseDataModal")
        } else if (!dir.exists(input$dataFolder)) {
            errorModal(
                session, "Folder not found", "The selected folder",
                tags$kbd(prepareWordBreak(input$dataFolder)), "was not found.",
                caller="Load TCGA data", modalId="firebrowseDataModal")
        } else if (!is.null(getData())) {
            loadedDataModal(session, "firebrowseDataModal", "firebrowseReplace",
                            "firebrowseAppend")
        } else {
            setFirebrowseData(input, output, session)
        }
    })

    # Load data when the user presses to replace data
    observeEvent(input$firebrowseReplace,
                 setFirebrowseData(input, output, session, replace=TRUE))

    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$firebrowseAppend,
                 setFirebrowseData(input, output, session, replace=FALSE))

    # Download data through the browser
    observeEvent(input$acceptDownload, {
        url <- getURLtoDownload()
        if (!is.null(url)) {
            display("Triggered file downloads")
            # Download missing files through the browser
            iframe <- function(url)
                tags$iframe(width=1, height=1, frameborder=0, src=url)
            output$iframeDownload <- renderUI(lapply(url, iframe))
            setURLtoDownload(NULL)
        }
    })
}

attr(firebrowseUI, "loader") <- "data"
attr(firebrowseServer, "loader") <- "data"
