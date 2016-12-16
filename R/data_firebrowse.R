#' Echo progress to console using \code{cat}
#' 
#' @param ... Strings to print to console
#' @param console Boolean: print to console? TRUE by default
#' @return NULL (this function is used to modify the Shiny session's state)
echoProgress <- function(..., console=TRUE) {
    if (console) cat(paste(...), fill=TRUE)
}

#' Returns the date format used by the Firehose API
#'
#' @return Named list with Firehose API's date formats
#'
#' @examples
#' format <- psichomics:::getFirehoseDateFormat()
#' 
#' # date format to use in a query to Firehose API
#' format$query
#' 
#' # date format to parse a date in a response from Firehose API
#' format$response
getFirehoseDateFormat <- function() {
    query <- "%Y_%m_%d"
    response <- "%d %b %Y" 
    return(list(query=query, response=response))
}

#' Parse the date from a response
#' @param string Character: dates
#' @return Parsed date
parseDateResponse <- function(string) {
    format <- getFirehoseDateFormat()$response
    date <- strsplit(string, " ")
    date <- lapply(date, function(i) paste(i[2:4], collapse=" "))
    date <- as.Date(unlist(date), format=format)
    return(date)
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
#' @importFrom httr GET warn_for_status http_error
#' @importFrom methods is
#'
#' @examples
#' isFirehoseUp()
isFirehoseUp <- function() {
    link <- paste0("http://firebrowse.org/api/v1/Metadata/HeartBeat")
    heartbeat <- tryCatch(GET(link, query=list(format="json")), error=return)
    if (is(heartbeat, "error")) {
        return(FALSE)
    } else if (http_error(heartbeat)) {
        warn_for_status(heartbeat, "reach Firehose API")
        return(FALSE)
    } else {
        return(TRUE)
    }
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
#'
#' @importFrom httr GET
#'
#' @examples
#' cohort <- psichomics:::getFirehoseCohorts()[1]
#' psichomics:::queryFirehoseData(cohort = cohort, data_type = "mRNASeq")
#' 
#' # Querying for data from a specific date
#' dates <- psichomics:::getFirehoseDates()
#' dates <- format(dates, psichomics:::getFirehoseDateFormat()$query)
#' 
#' psichomics:::queryFirehoseData(date = dates[2], cohort = cohort)
queryFirehoseData <- function(format = "json", date = NULL, cohort = NULL, 
                              data_type = NULL, tool = NULL, platform = NULL,
                              center = NULL, level = NULL, protocol = NULL,
                              page = NULL, page_size = NULL, sort_by = NULL) {
    # Only allow these response formats
    format <- match.arg(format, c("json", "csv", "tsv"))
    
    # Format date
    if (!is.null(date)) date <- gsub("-", "_", date, fixed=TRUE)
    
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
#'
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#'
#' @examples
#' psichomics:::parseFirehoseMetadata("Dates")
#' psichomics:::parseFirehoseMetadata("Centers")
#' psichomics:::parseFirehoseMetadata("HeartBeat")
#' 
#' # Get the abbreviation and description of all cohorts available
#' psichomics:::parseFirehoseMetadata("Cohorts")
#' # Get the abbreviation and description of the selected cohorts
#' psichomics:::parseFirehoseMetadata("Cohorts", cohort = c("ACC", "BRCA"))
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
#' if (isFirehoseUp()) getFirehoseDates()
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
#' if (isFirehoseUp()) getFirehoseCohorts()
getFirehoseCohorts <- function(cohort = NULL) {
    response <- parseFirehoseMetadata("Cohorts", cohort=cohort)
    cohorts <- response$Cohorts[[2]]
    names(cohorts) <- response$Cohorts[[1]]
    return(cohorts)
}

#' Download files to a given directory
#'
#' @param url Character: download links
#' @param folder Character: directory to store the downloaded archives
#' @param ... Extra parameters passed to the download function
#' @param download Function to use to download files
#' @param progress Function to show the progress (default is to print progress
#' to console)
#' 
#' @importFrom utils download.file
#' 
#' @return Invisible TRUE if every file was successfully downloaded
#'
#' @examples
#' \dontrun{
#' url <- paste0("https://unsplash.it/400/300/?image=", 570:572)
#' downloadFiles(url, "~/Pictures")
#' 
#' # Download without printing to console
#' downloadFiles(url, "~/Pictures", quiet = TRUE)
#' }
downloadFiles <- function(url, folder, progress = echoProgress,
                          download = download.file, ...) {
    destination <- file.path(folder, basename(url))
    for (i in seq_along(url)) {
        progress("Downloading file", detail = basename(url[i]), i, length(url))
        download(url[i], destination[i], ...)
    }
    cat("Downloading completed", fill=TRUE)
    return(destination)
}

#' Compute the 32-byte MD5 hashes of one or more files and check with given md5
#' file
#'
#' @param filesToCheck Character: files to calculate and match MD5 hashes
#' @param md5file Character: file containing correct MD5 hashes
#'
#' @importFrom digest digest
#' @importFrom utils read.table
#'
#' @return Logical vector showing TRUE for files with matching md5sums and FALSE
#' for files with non-matching md5sums
checkIntegrity <- function(filesToCheck, md5file) {
    if (is.na(md5file)) return(FALSE)
    md5sums <- digest(file = filesToCheck)
    md5table <- read.table(md5file, stringsAsFactors = FALSE)[[1]]
    return(md5sums %in% md5table)
}

#' Prepares Firehose archives in a given directory
#'
#' Checks Firehose archives' integrity using the MD5 files, extracts the content
#' of the archives, moves the content to newly-created folders and removes the 
#' original downloaded archives.
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
#'
#' @examples
#' file <- paste0(
#'     "~/Downloads",
#'     "ACC/20151101/gdac.broadinstitute.org_ACC.",
#'     "Merge_Clinical.Level_1.2015110100.0.0.tar.gz")
#' md5 <- paste0(file, ".md5")
#' \dontrun{
#' prepareFirehoseArchives(archive = file, md5 = paste0(file, ".md5"))
#' }
prepareFirehoseArchives <- function(archive, md5, folder, outdir) {
    # Check integrety of the downloaded archives with the MD5 files
    validFiles <- simplify2array(Map(checkIntegrity, archive, md5))
    
    ## TODO(NunoA): Should we try to download the invalid archives again?
    ## What if they're constantly invalid? Only try n times before giving up?
    if (!all(validFiles)) {
        warning("The MD5 hashes failed when checking the following files:\n",
                paste(archive[!validFiles], collapse = "\n\t"))
    }
    
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

#' Retrieve URLs from a response to a Firehose data query
#'
#' @param res Response from httr::GET to a Firehose data query
#'
#' @return Named character with URLs
#' 
#' @importFrom jsonlite fromJSON
#' @importFrom httr content
#'
#' @examples
#' res <- psichomics:::queryFirehoseData(cohort = "ACC")
#' url <- psichomics:::parseUrlsFromFirehoseResponse(res)
parseUrlsFromFirehoseResponse <- function(res) {
    # Parse the query response
    parsed <- content(res, "text", encoding = "UTF8")
    parsed <- fromJSON(parsed)[[1]]
    parsed$date <- parseDateResponse(parsed$date)
    
    # Get cohort names
    cohort <- getFirehoseCohorts()
    cohort <- cohort[parsed$cohort]
    
    ## TODO(NunoA): maybe this could be simplified?
    # Split URLs from response by cohort and datestamp
    url <- split(parsed$url, paste(cohort, format(parsed$date, "%Y-%m-%d")))
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
#' @include formats.R
#'
#' @param folder Character: folder(s) in which to look for Firehose files
#' @param exclude Character: files to exclude from the loading
#' @param progress Function to show the progress (default is to print progress
#' to console)
#' 
#' @return List with loaded data.frames
loadFirehoseFolders <- function(folder, exclude="", progress=echoProgress) {
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
        progress("Processing file", detail = basename(files[each]), each, 
                 length(files))
        loaded[[each]] <- parseValidFile(files[each], formats)
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    return(loaded)
}

#' Downloads and processes data from the Firehose API and loads it into R
#' 
#' @param folder Character: directory to store the downloaded archives (by
#' default, it saves in the user's "Downloads" folder)
#' @param data Character: data to load
#' @param exclude Character: files and folders to exclude from downloading and
#' from loading into R (by default, it excludes ".aux.", ".mage-tab." and
#' "MANIFEST.TXT" files)
#' @param ... Extra parameters to be passed to \code{\link{queryFirehoseData}}
#' @param progress Function to show the progress (default is to print progress
#' to console)
#' @param download Boolean: download missing files through the function
#' \code{download.file} (TRUE by default)
#' 
#' @include formats.R
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom httr stop_for_status
#' @importFrom utils download.file
#' 
#' @return URL of missing files ("missing" class) if files need to be downloaded
#' and if the argument \code{download} is \code{FALSE}; else, a list with loaded
#' data
#' @export
#' 
#' @examples 
#' \dontrun{
#' loadFirehoseData(cohort = "ACC", data_type = "Clinical")
#' }
loadFirehoseData <- function(folder=NULL, 
                             data=NULL, 
                             exclude=c(".aux.", ".mage-tab.", "MANIFEST.txt"),
                             ..., progress = echoProgress, download=TRUE) {
    args <- list(...)
    
    datasets <- unlist(getFirehoseDataTypes())
    # Data types to load
    args$data_type <- c(data[!data %in% datasets], "mRNASeq")
    # Datasets to ignore
    exclude <- c(exclude, datasets[!datasets %in% data])
    
    # Query Firehose and get URLs for archives
    res <- do.call(queryFirehoseData, args)
    stop_for_status(res)
    url <- parseUrlsFromFirehoseResponse(res)
    
    # Don't download specific items
    exclude <- paste(escape(exclude), collapse = "|")
    url <- url[!grepl(exclude, url)]
    
    # Get the file names without extensions
    md5  <- file_ext(url) == "md5"
    base <- basename(url)
    base[!md5] <- file_path_sans_ext(base[!md5], compression = TRUE)
    
    # Check which files are missing from the given directory
    if (is.null(folder)) folder <- getDownloadsFolder()
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
            progress(divisions = 1)
            cat("Triggered the download of files", fill=TRUE)
            
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
    
    # Check if there are folders to unarchive
    archives <- unlist(lapply(possibleExtensions, function (i)
        i[i %in% basename(downloadedFiles[!downloadedMD5])]))
    tar <- grepl(".tar", archives, fixed = TRUE)
    
    if (length(archives[tar]) > 0) {
        # Extract the content, check intergrity, move archives to newly-created
        # folders and remove original archives
        progress("Extracting archives...", divisions = 1 + length(folders))
        prepareFirehoseArchives(fullPath(archives[tar]), 
                                fullPath(base[md5][tar]), folder, folders)
        progress("Archives prepared")
    } else {
        # Set the progress bar to the number of folders to load
        progress(divisions = length(folders))   
    }
    
    # Get the full path of the files
    downloadedFiles <- list.files(folder, recursive=TRUE, full.names=TRUE, 
                                  include.dirs=TRUE)
    folders <- fullPath(base[!md5])
    folders <- split(folders, categories)
    
    # Load the files (but discard empty files with no rows)
    loaded <- lapply(folders, loadFirehoseFolders, exclude, progress)
    loaded <- lapply(loaded, function(i) Filter(nrow, i))
    loaded <- processDatasetNames(loaded)
    return(loaded)
}

#' Creates a UI set with options to add data from TCGA/Firehose
#' @param ns Namespace function
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny tagList uiOutput selectizeInput actionButton textAreaInput
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function(ns) {
    cohorts <- getFirehoseCohorts()
    acronyms <- names(cohorts)
    names(acronyms) <- sprintf("%s (%s)", cohorts, names(cohorts))
    
    dates <- as.character(getFirehoseDates())
    names(dates) <- dates
    names(dates)[1] <- paste(names(dates)[1], "(most recent)")
    
    dataTypes <- getFirehoseDataTypes()
    dataTypes[[1]] <- dataTypes[[1]][1]
    
    tagList(
        uiOutput(ns("firebrowseDataModal")),
        uiOutput(ns("pathSuggestions")),
        uiOutput(ns("iframeDownload")),
        selectizeInput(ns("firehoseCohort"), "Tumour type", acronyms,
                       multiple = TRUE,
                       options = list(placeholder = "Select cohort(s)")),
        selectizeInput(ns("firehoseDate"), "Date", dates, multiple = TRUE,
                       selected = NULL, options = list(
                           placeholder = "Select sample date")),
        selectizeInput(ns("firehoseData"), "Data type", multiple = TRUE,
                       c("Clinical data"="Clinical", dataTypes), 
                       options = list(placeholder = "Select data types")),
        textAreaInput(ns("dataFolder"), "Folder to store the data",
                      value = getDownloadsFolder(),
                      placeholder = "Insert data folder"),
        bsTooltip(ns("dataFolder"), placement = "right",
                  options = list(container = "body"),
                  "Data not available in this folder will be downloaded."),
        processButton(ns("getFirehoseData"), "Load data"))
}

#' User interface of the TCGA/Firebrowse loader
#' 
#' @param id Character: identifier
#' @param panel Function to enclose interface
#' 
#' @importFrom shiny NS
#' 
#' @return HTML of the interface
firebrowseUI <- function(id, panel) {
    ns <- NS(id)
    
    panel(style="info",
          title=list(icon("plus-circle"), "Load TCGA/Firehose data"),
          value="Load TCGA/Firehose data", uiOutput(ns("checkFirebrowse")))
}

#' Return an user interface depending on the status of the Firebrowse API
#' 
#' If the API is working, it'll be loaded. Else, a message will appear warning
#' the user that the API is down and that will let check again if the API is 
#' back online.
#' 
#' @param ns Namespace function
#' 
#' @importFrom shiny br icon tagList actionButton
#' 
#' @return HTML elements
checkFirebrowse <- function(ns) {
    startProgress("Checking Firebrowse API", 1)
    if (isFirehoseUp()) {
        updateProgress("Loading interface")
        ui <- addTCGAdata(ns)
    } else {
        ui <- tagList(icon("exclamation-circle"),
                      "Firebrowse seems to be offline at the moment.", br(), 
                      br(), actionButton(ns("refreshFirebrowse"),
                                         icon=icon("refresh"),
                                         "Check Firebrowse again", 
                                         class="btn-primary"))
    }
    closeProgress("Firebrowse interface loaded")
    return(ui)
}

#' Set data from Firehose
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param replace Boolean: replace loaded data? TRUE by default
#' 
#' @importFrom shinyjs disable enable
#' @importFrom shiny div fluidRow column icon
#' @importFrom shinyBS bsTooltip
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
setFirehoseData <- function(input, output, session, replace=TRUE) {
    ns <- session$ns
    time <- startProcess("getFirehoseData")
    
    # Load data from Firehose
    data <- loadFirehoseData(folder = input$dataFolder,
                             cohort = input$firehoseCohort,
                             date = gsub("-", "_", input$firehoseDate),
                             data = input$firehoseData,
                             progress = updateProgress,
                             download = FALSE)
    
    if (any(class(data) == "missing")) {
        updateProgress(divisions = 1)
        setURLtoDownload(data)
        
        infoModal(
            session, "Download requested data",
            "The requested data will be downloaded. When the downloads",
            "finish, click the button", tags$b("Load data"), 
            "again to process and load the downloaded data.", br(), br(), 
            tags$div(
                class="alert", class="alert-warning", role="alert",
                fluidRow(style="display: flex; align-items: center;", column(
                    10, "Confirm that files will be downloaded to the folder",
                    tags$b(input$dataFolder)),
                    column(2, tags$i(class="fa fa-question-circle", 
                                     id=ns("helpDownloadFolder"))))),
            bsTooltip(ns("helpDownloadFolder"), placement="right",
                      paste("This program checks for files in the given",
                            "folder. You can either:<br/>\u2022 Change the",
                            "folder where the downloaded items are located",
                            "<br/>\u2022 Move the downloaded items to the",
                            "given folder"),
                      options=list(container="body")),
            modalId="firebrowseDataModal",
            footer=actionButton(ns("acceptDownload"), "Download data",
                                class="btn-primary", "data-dismiss"="modal"))
    } else if (!is.null(data)) {
        if(replace)
            setData(data)
        else
            setData(c(getData(), data))
    }
    endProcess("getFirehoseData", time)
}

firebrowseServer <- function(input, output, session, active) {
    ns <- session$ns
    
    # If Firebrowse is unaccessible, allow user to try again
    output$checkFirebrowse <- renderUI(isolate(checkFirebrowse(ns)))
    observeEvent(input$refreshFirebrowse,
                 output$checkFirebrowse <- renderUI(checkFirebrowse(ns)))
    
    # # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))
    
    # Update available clinical data attributes to use in a formula
    output$pathSuggestions <- renderUI({
        checkInside <- function(path, showFiles=FALSE) {
            if (substr(path, nchar(path), nchar(path)) == "/")
                content <- list.files(path, full.names = TRUE)
            else
                content <- list.files(dirname(path), full.names = TRUE)
            
            # Show only directories if showFiles is FALSE
            if (!showFiles) content <- content[dir.exists(content)]
            return(basename(content))
        }
        
        textSuggestions(ns("dataFolder"), checkInside(input$dataFolder),
                        char=.Platform$file.sep)
    })
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input$getFirehoseData, {
        if (length(isolate(input$firehoseCohort)) == 0) {
            errorModal(session, "No tumour type",
                       "Please, input a tumour type.", 
                       modalId="firebrowseDataModal")
        } else if (length(isolate(input$firehoseDate)) == 0) {
            errorModal(session, "No date",
                       "Please, input date of samples of interest.",
                       modalId="firebrowseDataModal")
        } else if (length(isolate(input$firehoseData)) == 0) {
            errorModal(session, "No data types",
                       "Please, input data types of interest.",
                       modalId="firebrowseDataModal")
        } else if (!is.null(getData())) {
            loadedDataModal(session,
                            "firebrowseDataModal",
                            "firebrowseReplace",
                            "firebrowseAppend")
        } else {
            setFirehoseData(input, output, session)
        }
    })
    
    # Load data when the user presses to replace data
    observeEvent(input$firebrowseReplace,
                 setFirehoseData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$firebrowseAppend,
                 setFirehoseData(input, output, session, replace=FALSE))
    
    # Download data through the browser
    observeEvent(input$acceptDownload, {
        url <- getURLtoDownload()
        if (!is.null(url)) {
            cat("Triggered file downloads", fill=TRUE)
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