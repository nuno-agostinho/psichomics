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
    return(warn_for_status(heartbeat))
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
#' queryFirehoseData(cohort = "ACC", data_type = "mRNASeq")
#' 
#' # Querying for data from a specific date
#' dates <- getFirehoseDates()
#' queryFirehoseData(date = dates[2], cohort = "ACC")
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
    response <- GET("http://firebrowse.org",
                    path = "api/v1/Archives/StandardData", query = query)
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
    response <- fromJSON(content(response, "text"))
    return(response)
}

#' Query the Firehose API for the datestamps of the data available and parse the
#' response
#'
#' @return Character with datestamps of the data available
#' @export
#'
#' @examples
#' getFirehoseDates()
getFirehoseDates <- function() {
    dates <- parseFirehoseMetadata("Dates")
    return(dates$Dates)
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
#' @param urls Chareacter: download links
#' @param folder Character: directory to store the downloaded archives
#' @param ... Extra parameters passed to the function download.file
#'
#' @return Invisible TRUE if every file was successfully downloaded
#' @export
#'
#' @examples
#' urls <- paste0("https://unsplash.it/400/300/?image=", 570:572)
#' downloadFiles(urls, "~/Pictures")
#' 
#' # Download without printing to console
#' downloadFiles(urls, "~/Pictures", quiet = TRUE)
downloadFiles <- function(urls, folder, ...) {
    destination <- file.path(folder, basename(urls))
    for (i in seq_along(urls))
        download.file(urls[i], destination[i], ...)
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
    md5sums <- tools::md5sum(filesToCheck)
    md5table <- read.table(md5file, stringsAsFactors = FALSE)[[1]]
    return(md5sums %in% md5table)
}

#' Downloads and processes data from the Firehose API and loads it into R
#' 
#' @param folder Character: directory to store the downloaded archives (by
#' default, it saves in the user's "Downloads" folder)
#' @param ... Extra parameters to be passed to \code{\link{queryFirehoseData}}
#' 
#' @export
#' @examples 
#' loadFirehoseData
loadFirehoseData <- function(folder = "~/Downloads", ...) {
    ## TODO(NunoA): Check if the default folder works in Windows
    # Query Firehose with the user query
    res <- queryFirehoseData(...)
    stop_for_status(res)
    
    # Parse the query
    res <- fromJSON(content(res, "text"))[[1]]
    id <- paste(res$cohort, res$data_type, res$protocol)
    urls <- res$urls
    names(urls) <- id
    
    # Download the appropriate archives into a specific location
    archives <- downloadFiles(unlist(urls), folder)
    
    # Check integrety of the downloaded archives
    filesOfInterest <- archives[tools::file_ext(archives) != "md5"]
    print(filesOfInterest)
    return(filesOfInterest)
    filesAreValid <- all(
        vapply(filesOfInterest,
               function(i) checkIntegrity(i, paste0(i, ".md5")),
               logical(1)))
    
    if (filesAreValid){
        print(filesOfInterest)
        # Extract the contents of the archives
        # Remove the original archives
        # Load the files into R using readr (faster and shows progress)
    } else stop("Error: files are not valid according to the MD5 hashes.")
}