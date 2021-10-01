#' Get GTEx data information
#'
#' @family functions associated with GTEx data retrieval
#' @return GTEx data information
#' @export
#'
#' @examples
#' getGtexDataTypes()
getGtexDataTypes <- function() {
    c("Sample attributes"="sampleInfo",
      "Subject phenotypes"="subjectInfo",
      "Gene expression"="geneExpr",
      "Junction quantification"="junctionQuant")
}

#' @rdname getGtexDataTypes
#' @export
#'
#' @examples
#' getGtexReleases()
getGtexReleases <- function() {
    release <- c(8, 7, 6, 4)
    n <- max(release) + 1
    getReleaseIfNotNull <- function(n) if (!is.null(getGtexDataURL(n))) n
    release <- c(release, getReleaseIfNotNull(n), getReleaseIfNotNull(n + 1))
    return(release)
}

#' @rdname appUI
#'
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom shiny helpText
#' @importFrom shinyjs hidden
gtexDataUI <- function(id, panel) {
    ns <- NS(id)
    panel(style="info", title=list(icon("plus-circle"), "GTEx data loading"),
          value="Automatically load GTEx data",
          uiOutput(ns("modal")),
          helpText("GTEx data are downloaded from the",
                   a(href="http://www.gtexportal.org", target="_blank",
                     "GTEx Data Portal"), "website."),
          selectizeInput(ns("release"), "Version release", width="100%",
                         getGtexReleases(), selected=8),
          selectizeInput(ns("dataTypes"), "Data type", multiple=TRUE,
                         width="100%", getGtexDataTypes(),
                         selected=getGtexDataTypes(), options=list(
                             placeholder="Select data types",
                             plugins=list("remove_button"))),
          browseDownloadFolderInput(ns("folder")),
          bsCollapse(
              id=ns("filterCollapse"),
              bsCollapsePanel(
                  title=tagList(icon("lungs"), "Tissues to load",
                                contextUI(ns("filterText"))),
                  value="Load by tissue",
                  div(id=ns("loadingAvailableTissues"), class="progress",
                      div(class="progress-bar progress-bar-striped active",
                          role="progressbar", style="width:100%",
                          "Loading tissues from sample attributes...")),
                  hidden(
                      selectizeInput(ns("tissues"), label=NULL, width="100%",
                                     choices=c("Select one or more tissues"=""),
                                     multiple=TRUE, options=list(
                                         plugins=list("remove_button")))))),
          processButton(ns("load"), "Load data"))
}

#' Get GTEx tissues from given GTEx sample attributes
#'
#' @inheritParams loadGtexData
#'
#' @family functions associated with GTEx data retrieval
#' @return Character: available tissues
#' @export
#'
#' @examples
#' \dontrun{
#' getGtexTissues()
#' }
getGtexTissues <- function(folder=getDownloadsFolder(),
                           release=getGtexReleases()[[1]]) {
    metadata <- loadGtexData(folder, "sampleInfo", release=release,
                             progress=FALSE)[[1]][[1]]
    freq     <- table(metadata[["Tissue Type (area of retrieval)"]])
    tissues  <- names(freq)
    names(tissues) <- sprintf("%s (%s samples)", names(freq), as.vector(freq))
    if (!is.na(match("", tissues))) tissues <- tissues[-match("", tissues)]
    return(tissues)
}

#' Load GTEx file
#'
#' @param path Character: path to file
#' @param pattern Character: pattern of the format type to load file
#' @param samples Character: samples to filter datasets
#'
#' @return Loaded file as a data frame
#' @keywords internal
loadGtexFile <- function(path, pattern, samples=NULL) {
    if (!is.null(path)) {
        if (!is.character(path) && !is.null(path$datapath))
            path <- path$datapath
        else
            path <- path
    } else {
        return(NULL)
    }

    # Retrieve correct format to load GTEx file
    formats <- loadFileFormats()
    filterFormats <- function(i, pattern) {
        if (!is.null(i$filename)) {
            grepl(pattern, i$filename, fixed=TRUE) &&
                grepl("GTEx", i$filename, fixed=TRUE)
        } else {
            return(FALSE)
        }
    }
    format <- formats[sapply(formats, filterFormats, pattern)]

    select <- NULL
    if (pattern %in% c("junction", "gene") && !is.null(samples)) {
        # Parse all columns instead of specific ones
        customFormat <- format
        for (i in seq(customFormat)) customFormat[[i]]$ignoreCols <- NULL

        # Load GTEx data exclusively for matching samples
        allSamples <- colnames(loadFile(path, customFormat, nrows=0))
        select <- c(1, # Retrieve junction identifier
                    which(allSamples %in% samples))
    }
    parsed <- loadFile(path, format, select=select)

    if (!is.null(samples)) {
        if (pattern == "Sample") {
            # Retrieve samples based on tissues
            parsed <- parsed[samples, ]
        } else if (pattern == "Subject") {
            # Retrieve subjects for which samples are available
            subjects <- getSubjectFromSample(samples, rownames(parsed))
            subjects <- subjects[!is.na(subjects)]
            subjects <- sort(unique(subjects))
            parsed <- parsed[subjects, ]
        }
    }
    return(parsed)
}

downloadGtexFiles <- function(link, folder) {
    type            <- names(link)
    filepath        <- file.path(folder, basename(link))
    names(filepath) <- type
    toDownload <- !file.exists(gsub("\\.gz$", "", filepath))
    if (sum(toDownload) > 0) {
        updateProgress("Downoading data...", divisions=sum(toDownload))
        for (i in which(toDownload)) {
            updateProgress("Downloading file", detail=type[i])
            download.file(link[i], filepath[i])
        }
    }
    toDecompress <- grepl("\\.gz$", filepath) & file.exists(filepath)
    if (sum(toDecompress) > 0) {
        updateProgress("Extracting files...", divisions=sum(toDecompress))
        for (each in which(toDecompress)) {
            updateProgress("Extracting file", detail=type[each])
            gunzip(filepath[each])
        }
    }
    return(filepath)
}

#' Get links to download GTEx data
#'
#' @param release Numeric: GTEx data release
#' @param domain Character: GTEx data storage domain
#' @param offline Boolean: simulate offline behaviour
#'
#' @importFrom XML xmlParse xmlToDataFrame
#' @importFrom httr GET http_error
#'
#' @return Character with URLs to download GTEx data
#' @keywords internal
getGtexDataURL <- function(release, domain="https://storage.googleapis.com",
                           offline=FALSE) {
    path <- paste0("gtex_analysis_v", release)
    resp <- try(GET(domain, path=path, timeout(3)))
    date <- NULL
    if (!is(resp, "try-error") && !http_error(resp) && !offline) {
        doc <- xmlParse(resp)
        df  <- xmlToDataFrame(doc, nodes=xmlRoot(doc)[-c(seq(4))],
                              stringsAsFactors=FALSE)
        files  <- c("annotations/.*Annotations_SampleAttributesDS\\.txt",
                    "annotations/.*Annotations_SubjectPhenotypes.*DS\\.txt",
                    "rna_seq_data/.*_gene_reads\\.gct\\.gz",
                    "rna_seq_data/.*junction.*\\.gz")
        index <- sapply(files, grep, df[[1]])
        if (length(Filter(length, index)) == 0) return(NULL)
        res   <- df[index, "Key"]
        date  <- max(as.Date(df[index, "LastModified"]))
    } else if (release == 8) {
        res <- c(
            "annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
            "annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
            paste0("rna_seq_data/GTEx_Analysis_2017-06-05_v8_",
                   c("RNASeQCv1.1.9_gene_reads.gct.gz",
                     "STARv2.5.3a_junctions.gct.gz")))
    } else if (release == 7) {
        res <- c("annotations/GTEx_v7_Annotations_SampleAttributesDS.txt",
                 "annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
                 paste0("rna_seq_data/GTEx_Analysis_2016-01-15_v7_",
                        c("RNASeQCv1.1.8_gene_reads.gct.gz",
                          "STARv2.4.2a_junctions.gct.gz")))
    } else if (release == 6) {
        res <- c("annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
                 "annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",
                 paste0("rna_seq_data/GTEx_Analysis_v6_RNA-seq_",
                        c("RNA-SeQCv1.1.8_gene_reads.gct.gz",
                          "Flux1.6_junction_reads.txt.gz")))
    } else if (release == 4) {
        res <- c(
            "annotations/GTEx_Data_V4_Annotations_SampleAttributesDS.txt",
            "annotations/GTEx_Data_V4_Annotations_SubjectPhenotypes_DS.txt",
            paste0("rna_seq_data/GTEx_Analysis_V4_RNA-seq_",
                   c("RNA-SeQCv1.1.8_gene_reads.gct.gz",
                     "Flux1.6_junction_reads.txt.gz")))
    } else {
        res <- NULL
    }
    if (!is.null(res)) {
        res               <- file.path(domain, path, res)
        names(res)        <- getGtexDataTypes()
        attr(res, "date") <- date

    }
    return(res)
}

#' Download and load GTEx data
#'
#' @param data Character: data types to load (see \code{getGtexDataTypes})
#' @param folder Character: folder containing data
#' @param tissue Character: tissues to load (if \code{NULL}, load all); tissue
#' selection may speed up data loading
#' @param release Numeric: GTEx data release to load
#' @param progress Boolean: display progress?
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom R.utils gunzip
#'
#' @family functions associated with GTEx data retrieval
#' @family functions to load data
#' @return List with loaded data
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and load all available GTEx data
#' data <- loadGtexData()
#'
#' # Download and load only junction quantification and sample info from GTEx
#' getGtexDataTypes()
#' data <- loadGtexData(data=c("sampleInfo", "junctionQuant"))
#'
#' # Download and load only data for specific tissues
#' getGtexTissues()
#' data <- loadGtexData(tissue=c("Stomach", "Small Intestine"))
#'
#' # Download and load data from a specific GTEx data release
#' data <- loadGtexData(tissue=c("Stomach", "Small Intestine"), release=7)
#' }
loadGtexData <- function(folder=getDownloadsFolder(), data=getGtexDataTypes(),
                         tissue=NULL, release=getGtexReleases()[[1]],
                         progress=TRUE) {
    stopifnot("Argument 'data' cannot be NULL"   = !is.null(data),
              "Argument 'folder' cannot be NULL" = !is.null(folder))

    link <- getGtexDataURL(release)[data]
    if (is.null(link)) stop("No data available for GTEx V", release)

    folder <- file.path(folder, paste0("GTEx_V", release))
    if (!dir.exists(folder)) dir.create(folder)
    filepath <- downloadGtexFiles(link, folder)

    if (progress) updateProgress("Loading files...", divisions=length(data))
    loadThisGtexFile <- function(path, pattern, release, samples=NULL) {
        name <- ifelse(is.character(path), basename(path), path$name)
        if (progress) updateProgress("Processing file", detail=name)
        loaded <- loadGtexFile(path, pattern, samples)
        attr(loaded, "description") <- gsub(
            "GTEx", paste0("GTEx V", release), attr(loaded, "description"),
            fixed=TRUE)
        return(loaded)
    }
    filepath <- gsub("\\.gz$", "", filepath)
    sampleMetadata <- filepath["sampleInfo"]
    clinical       <- filepath["subjectInfo"]
    geneExpr       <- filepath["geneExpr"]
    junctionQuant  <- filepath["junctionQuant"]

    loaded <- list()
    samples <- NULL
    if (!is.na(sampleMetadata)) {
        sampleAttrs <- loadThisGtexFile(sampleMetadata, "Sample", release)
        if (!is.null(tissue)) {
            # Filter which samples match desired tissues
            tissue         <- tolower(unique(tissue))
            tissueCol      <- "Tissue Type (area of retrieval)"
            sampleTissues  <- tolower(sampleAttrs[[tissueCol]])
            matchedTissues <- sampleTissues %in% tissue
            samples        <- rownames(sampleAttrs)[matchedTissues]

            if (length(samples) == 0) {
                noMatch <- paste(tissue[-length(tissue)], collapse=", ")
                if (noMatch == "") {
                    noMatch <- tissue
                } else {
                    noMatch <- paste(noMatch, tissue[length(tissue)],
                                     sep=" or ")
                }
                stop("No tissues match ", noMatch, ". Run the function ",
                     "`getGtexTissues` to check available tissues.")
            }
            sampleAttrs <- sampleAttrs[samples, ]
        }
        loaded[[1]] <- sampleAttrs
    }

    if (!is.na(clinical)) {
        loaded[[2]] <- loadThisGtexFile(clinical, "Subject", release, samples)
    }
    if (!is.na(junctionQuant)) {
        loaded[[3]] <- loadThisGtexFile(junctionQuant, "junction", release,
                                        samples)
    }
    if (!is.na(geneExpr)) {
        loaded[[4]] <- loadThisGtexFile(geneExpr, "gene", release, samples)
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    attr(loaded, "source") <- paste0("GTEx v", release)

    gtex <- setNames(list(loaded), paste0("GTEx_V", release))
    gtex <- processDatasetNames(gtex)
    if (progress) closeProgress()
    return(gtex)
}

#' Shiny wrapper to load GTEx data
#'
#' @param session Shiny session
#' @inheritParams appServer
#' @param replace Boolean: replace loaded data?
#'
#' @inherit psichomics return
#' @keywords internal
loadGtexDataShiny <- function(session, input, replace=TRUE) {
    dataTypes <- input$dataTypes
    folder    <- input$folder
    tissue    <- input$tissues
    release   <- input$release

    time <- startProcess("load")
    data <- loadGtexData(folder, dataTypes, tissue, release)

    if (!is.null(data)) {
        if (!replace) data <- c(getData(), data)
        data <- processDatasetNames(data)
        setData(data)
    }
    endProcess("load", time)
}

#' @rdname appServer
#' @importFrom shinyjs show hide
gtexDataServer <- function(input, output, session) {
    prepareFileBrowser(session, input, "folder", directory=TRUE)

    observeEvent(input$load, {
        dataTypes <- input$dataTypes
        folder    <- input$folder
        release   <- input$release

        if (is.null(dataTypes)) {
            errorModal(session, "No data types selected",
                       "Please, select one or more data types.",
                       caller="Load GTEx data")
        } else if (is.null(folder) || folder == "") {
            errorModal(session, "No folder selected",
                       "Please, select a folder.",
                       caller="Load GTEx data")
        } else if (is.null(release) || release == "") {
            errorModal(session, "No GTEx release version selected",
                       "Please, select a GTEx release.",
                       caller="Load GTEx data")
        } else if (!is.null(getData())) {
            loadedDataModal(session, "modal", "replace", "append")
        } else {
            loadGtexDataShiny(session, input)
        }
    })

    # Select available tissues from GTEx
    showAvailableTissues <- reactive({
        folder       <- input$folder
        release      <- input$release
        progressBar  <- "loadingAvailableTissues"
        tissueSelect <- "tissues"

        fadeIn  <- function(id, ...) show(id, anim=TRUE, ...)
        fadeOut <- function(id, ...) hide(id, anim=TRUE, ...)

        fadeOut(progressBar)
        fadeOut(tissueSelect)
        tissues <- tryCatch(getGtexTissues(folder, release),
                            error=return, warning=return)
        suppressWarnings(closeProgress())
        fadeIn(tissueSelect, animType="fade")
        tissues <- c(tissues, "Select one or more tissues"="")
        return(tissues)
    })

    observe({
        if (!identical(input$filterCollapse, "Load by tissue")) return(NULL)
        tissues <- showAvailableTissues()
        updateSelectizeInput(session, "tissues", choices=tissues)
    })

    # Update number of tissues selected in context
    output$filterText <- renderText({
        tissues <- input$tissues
        if (is.null(tissues) || length(tissues) == 0) {
            text <- "Load all tissues"
        } else {
            len  <- length(tissues)
            text <- sprintf("Load %s tissue%s", len, ifelse(len == 1, "", "s"))
        }
        return(text)
    })

    # Replace or append data to existing data
    observeEvent(input$replace, loadGtexDataShiny(session, input, replace=TRUE))
    observeEvent(input$append, loadGtexDataShiny(session, input, replace=FALSE))
}

attr(gtexDataUI, "loader") <- "data"
attr(gtexDataServer, "loader") <- "data"
