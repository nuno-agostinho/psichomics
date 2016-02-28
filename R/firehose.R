#' @import httr tools
#' @importFrom jsonlite fromJSON
NULL

#' Returns the date format used by the Firehose API
#'
#' @return Named list with Firehose API's date formats
#' @export
#'
#' @examples
#' format <- getFirehoseDateFormat()
#' 
#' # date format to use in a query to Firehose API
#' format$query
#' 
#' # date format to parse a date in a response from Firehose API
#' format$response
getFirehoseDateFormat <- function() {
    query <- "%Y_%m_%d"
    response <- "%a, %d %b %Y %H:%M:%S" 
    return(list(query=query, response=response))
}

#' Check if the Firehose API is running
#'
#' The Firehose API is running if it returns the status condition 200; if
#' this is not the status code obtained from the API, the function will raise a
#' warning with the status code and a brief explanation.
#'
#' @return Invisible TRUE if the Firehose API is working; otherwise, raises a
#' warning
#' @export
#'
#' @examples
#' isFirehoseUp()
isFirehoseUp <- function() {
    link <- paste0("http://firebrowse.org/api/v1/Metadata/HeartBeat")
    heartbeat <- GET(link, query = list(format = "json"))
    if (http_error(heartbeat))
        return(warn_for_status(heartbeat, "reach Firehose API"))
    else
        return(invisible(TRUE))
}

#' Query the Firehose API for TCGA data
#'
#' @param format Character: response format as JSON (default), CSV or TSV
#' @param date Character: dates of the data retrieval by Firehose (by default,
#' it uses the most recent data available)
#' @param cohort Character: abbreviation of the cohorts (by default, returns
#' data for all cohorts)
#' @param data_type Character: data types (optional)
#' @param tool Character: data produced by the selected Firehose tools
#' (optional)
#' @param platform Character: data generation platforms (optional)
#' @param center Character: data generation centers (optional)
#' @param level Integer: data levels (optional)
#' @param protocol Character: sample characterization protocols (optional)
#' @param page Integer: page of the results to return (optional)
#' @param page_size Integer: number of records per page of results; max is 2000
#' (optional)
#' @param sort_by String: column used to sort the data (by default, it sorts by
#' cohort)
#'
#' @return Response from the Firehose API (it needs to be parsed)
#' @export
#'
#' @examples
#' cohort <- getFirehoseCohorts()[1]
#' queryFirehoseData(cohort = cohort, data_type = "mRNASeq")
#' 
#' # Querying for data from a specific date
#' dates <- getFirehoseDates()
#' dates <- format(dates, getFirehoseDateFormat()$query)
#' 
#' queryFirehoseData(date = dates[2], cohort = cohort)
queryFirehoseData <- function(format = "json", date = NULL, cohort = NULL, 
                              data_type = NULL, tool = NULL, platform = NULL,
                              center = NULL, level = NULL, protocol = NULL,
                              page = NULL, page_size = NULL, sort_by = NULL) {
    # Only allow these response formats
    format <- match.arg(format, c("json", "csv", "tsv"))
    
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

#' Query the Firehose API for metadata and parse the response
#'
#' @param type Character: metadata to retrieve
#' @param ... Character: parameters to pass to query (optional)
#'
#' @return List with parsed JSON response
#' @export
#'
#' @examples
#' parseFirehoseMetadata("Dates")
#' parseFirehoseMetadata("Centers")
#' parseFirehoseMetadata("HeartBeat")
#' 
#' # Get the abbreviation and description of all cohorts available
#' parseFirehoseMetadata("Cohorts")
#' # Get the abbreviation and description of the selected cohorts
#' parseFirehoseMetadata("Cohorts", cohort = c("ACC", "BRCA"))
parseFirehoseMetadata <- function(type, ...) {
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

#' Query the Firehose API for the datestamps of the data available and parse the
#' response
#'
#' @return Date with datestamps of the data available
#' @export
#'
#' @examples
#' getFirehoseDates()
getFirehoseDates <- function() {
    dates <- parseFirehoseMetadata("Dates")$Dates
    format <- getFirehoseDateFormat()
    dates <- as.Date(dates, format$query)
    return(dates)
}

#' Query the Firehose API for the cohorts available
#'
#' @param cohort Character: filter by given cohorts (optional)
#'
#' @return Character with cohort abbreviations (as values) and description (as 
#' names)
#' @export
#'
#' @examples
#' getFirehoseCohorts()
getFirehoseCohorts <- function(cohort = NULL) {
    response <- parseFirehoseMetadata("Cohorts", cohort=cohort)
    cohorts <- response$Cohorts[[1]]
    names(cohorts) <- response$Cohorts[[2]]
    return(cohorts)
}

#' Download files to a given directory
#'
#' @param url Character: download links
#' @param folder Character: directory to store the downloaded archives
#' @param ... Extra parameters passed to the download function
#' @param download Function to use to download files
#' @param progress Function to show the progress (default is function(...) 
#' print(paste(...)))
#' 
#' @return Invisible TRUE if every file was successfully downloaded
#' @export
#'
#' @examples
#' url <- paste0("https://unsplash.it/400/300/?image=", 570:572)
#' downloadFiles(url, "~/Pictures")
#' 
#' # Download without printing to console
#' downloadFiles(url, "~/Pictures", quiet = TRUE)
downloadFiles <- function(url, folder, progress = function(...) print(paste(...)),
                          download = download.file, ...) {
    destination <- file.path(folder, basename(url))
    for (i in seq_along(url)) {
        progress("Downloading file", detail = basename(url[i]), i, length(url))
        download(url[i], destination[i], ...)
    }
    print("Downloading completed")
    return(destination)
}

#' Compute the 32-byte MD5 hashes of one or more files and check with given md5
#' file
#'
#' @param filesToCheck Character: files to calculate and match MD5 hashes
#' @param md5file Character: file containing correct MD5 hashes
#'
#' @return Logical vector showing TRUE for files with matching md5sums and FALSE
#' for files with non-matching md5sums
#' @export
checkIntegrity <- function(filesToCheck, md5file) {
    md5sums <- digest::digest(file = filesToCheck)
    md5table <- read.table(md5file, stringsAsFactors = FALSE)[[1]]
    return(md5sums %in% md5table)
}

#' Prepares Firehose archives in a given directory
#'
#' Checks Firehose archives' integrity using the MD5 files, extracts the content
#' of the archives and removes the original downloaded archives.
#'
#' @param downloaded Character: path to downloaded archives
#' @param folder Character: local folder where the archives should be stored
#' @param progress Function to show the progress (default is function(...) 
#' print(paste(...)))
#' 
#' @return Invisible TRUE if successful
#' @export
#'
#' @examples
#' prepareFirehoseArchives(folder = "~/Downloads", url = paste0(
#'     "http://gdac.broadinstitute.org/runs/stddata__2015_11_01/data/",
#'     "ACC/20151101/gdac.broadinstitute.org_ACC.",
#'     "Merge_Clinical.Level_1.2015110100.0.0.tar.gz", c("", ".md5")))
prepareFirehoseArchives <- function (downloaded, folder,
                                     progress = function(...) print(paste(...))) {
    # Check integrety of the downloaded archives with the MD5 files
    downloadedFolders <- downloaded[tools::file_ext(downloaded) != "md5"]
    ## TODO(NunoA): don't assume every file has the respective MD5 file
    validFiles <- simplify2array(Map(checkIntegrity, downloadedFolders,
                                     paste0(downloadedFolders, ".md5")))
    
    ## TODO(NunoA): Should we try to download the invalid archives again?
    ## What if they're constantly invalid? Only try n times before giving up?
    if (!all(validFiles))
        stop(paste("Error: at least one file is not valid according to the",
                   "MD5 hashes."))
    
    ## TODO(NunoA): Check if path.expand works in Windows
    # Extract the contents of the archives to the same folder
    invisible(lapply(downloadedFolders,
                     untar, exdir = path.expand(folder)))
    
    # Remove the original downloaded files
    invisible(file.remove(downloaded))
    return(invisible(TRUE))
}

#' Retrieve URLs from a response to a Firehose data query
#'
#' @param res Response from httr::GET to a Firehose data query
#'
#' @return Named character with URLs
#' @export
#'
#' @examples
#' res <- queryFirehoseData(cohort = "ACC")
#' url <- parseUrlsFromResponse(res)
parseUrlsFromFirehoseResponse <- function(res) {
    # Parse the query response
    parsed <- content(res, "text", encoding = "UTF8")
    parsed <- fromJSON(parsed)[[1]]
    parsed$date <- as.Date(parsed$date,
                           format = getFirehoseDateFormat()$response)
    
    ## TODO(NunoA): maybe this could be simplified?
    # Split URLs from response by cohort and datestamp
    url <- split(parsed$url,
                 paste(parsed$cohort, format(parsed$date, "%Y-%m-%d")))
    url <- lapply(url, unlist)
    link <- unlist(url)
    names(link) <- rep(names(url), vapply(url, length, numeric(1)))
    return(link)
}

#' Load Firehose folders
#'
#' Loads the files present in each folder as a data.frame.
#' 
#' @note For faster execution, this function uses the \code{readr} library. This
#' function ignores subfolders of the given folder (which means that files 
#' inside subfolders are NOT loaded).
#'
#' @param folder Character: folder(s) in which to look for Firehose files
#' @param exclude Character: files to exclude from the loading
#' @param progress Function to show the progress (default is function(...) 
#' print(paste(...)))
#' 
#' @return List with loaded data.frames
#' @export
#'
#' @examples
#' # Load files from "~/Downloads"
#' loadFirehoseFolders("~/Downloads")
#' 
#' # Load files from folders inside "~/Downloads"
#' folders <- list.dirs("~/Downloads")
#' loadFirehoseFolders(folders)
#' 
#' # Exclude certain files from being loaded
#' loadFirehoseFolders(folders, exclude = c("pink.txt", "panther.txt"))
loadFirehoseFolders <- function (folder, exclude="",
                                 progress = function(...) print(paste(...))) {
    # Retrieve full path of the files inside the given folders
    files <- dir(folder, full.names=TRUE)
    
    # Exclude subdirectories and undesired files
    files <- files[!dir.exists(files)]
    exclude <- paste(exclude, collapse = "|")
    if (exclude != "") files <- files[!grepl(exclude, files)]
    
    # Try to load files and remove those with 0 rows
    loaded <- list()
    for (each in seq_along(files)) {
        progress("Processing file", detail = basename(files[each]),
                 each, length(files))
        loaded[[each]] <- parseValidFile(files[each], "R/formats")
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    return(loaded)
}

#' Downloads and processes data from the Firehose API and loads it into R
#' 
#' @param folder Character: directory to store the downloaded archives (by
#' default, it saves in the user's "Downloads" folder)
#' @param exclude Character: files and folders to exclude from downloading and
#' from loading into R (by default, it excludes ".aux.", ".mage-tab." and
#' "MANIFEST.TXT" files)
#' @param ... Extra parameters to be passed to \code{\link{queryFirehoseData}}
#' @param progress Function to show the progress (default is function(...) 
#' print(paste(...)))
#' @param download Function to download the files (default is download.file)
#' 
#' @export
#' @examples 
#' loadFirehoseData()
loadFirehoseData <- function(folder = "~/Downloads",
                             exclude = c(".aux.", ".mage-tab.", "MANIFEST.txt"),
                             ..., progress = function(...) print(paste(...)),
                             download = download.file) {
    # Check if folder exists
    if (!dir.exists(folder)) stop("Directory doesn't exist!")
    
    ## TODO(NunoA): Check if the default folder works in Windows
    # Query Firehose and get URLs for archives
    res <- queryFirehoseData(...)
    stop_for_status(res)
    url <- parseUrlsFromFirehoseResponse(res)
    
    # Ignore specific archives
    exclude <- paste(exclude, collapse = "|")
    url <- url[!grepl(exclude, url)]
    
    # Check which folders have already been downloaded to the given directory
    noMD5 <- gsub(".md5", "", url)
    base <- file_path_sans_ext(basename(noMD5), compression = TRUE)
    archives <- file.path(folder, base)
    missing <- url[!file.exists(archives)]
    
    # Evenly divide the progress bar in one (download) + number of files to load
    md5 <- grepl(".md5", url)
    archives <- split(archives[!md5], names(url[!md5]))
    progress(divisions = (1 + length(archives)))
    
    # Download and prepare archives not present in the given directory
    if (length(missing) > 0) {
        downloaded <- downloadFiles(missing, folder, progress)
        prepareFirehoseArchives(downloaded, folder, progress)
    } else {
        progress("Archives already downloaded")
    }
    
    ## TODO(NunoA): Check if it's possible to show READR progress in a Shiny app
    # Load the files using readr (faster and can show progress)
    loaded <- lapply(archives, loadFirehoseFolders, exclude, progress)
    return(loaded)
}
