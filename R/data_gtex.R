#' Get GTEx data types
#' 
#' @family functions associated with GTEx data retrieval
#' @return GTEx data types
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
          selectizeInput(ns("dataTypes"), "Data type", multiple=TRUE,
                         width="100%", getGtexDataTypes(), 
                         selected=names(getGtexDataTypes()), options=list(
                             placeholder="Select data types",
                             plugins=list("remove_button"))),
          fileBrowserInput(
              ns("folder"), "Folder where data is stored",
              value=getDownloadsFolder(),
              placeholder="No folder selected",
              info=TRUE, infoFUN=bsTooltip,
              infoTitle=paste("Data will be downloaded if not available in",
                              "this folder.")),
          bsCollapse(
              id=ns("filterCollapse"),
              bsCollapsePanel(
                  title=tagList(icon("filter"), "Filter tissues to load"), 
                  value="Load by tissue",
                  div(id=ns("loadingAvailableTissues"), class="progress",
                      div(class="progress-bar progress-bar-striped active",
                          role="progressbar", style="width:100%",
                          "Loading tissues from sample attributes...")),
                  hidden(
                      selectizeInput(ns("tissues"), label=NULL, width="100%",
                                     choices=c("Select one or more tissues"=""),
                                     multiple=TRUE)))),
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
getGtexTissues <- function(folder=getDownloadsFolder()) {
    sampleFile <- "GTEx_v7_Annotations_SampleAttributesDS.txt"
    filepath <- file.path(folder, sampleFile)
    names(filepath) <- sampleFile
    downloadGtexFiles(filepath, "Samples")
    
    filepath <- gsub("\\.gz$", "", filepath)
    sampleMetadata <- grep("Sample", filepath, value=TRUE)
    metadata <- loadGtexFile(sampleMetadata, "Sample")
    freq     <- table(metadata[["Tissue Type (area of retrieval)"]])
    
    tissues        <- names(freq)
    names(tissues) <- sprintf("%s (%s samples)", names(freq), as.vector(freq))
    if (!is.na(match("", tissues)))
        return(tissues[-match("", tissues)])
    else
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
        allSamples <- colnames(parseValidFile(path, customFormat, nrows=0))
        select <- c(1, # Retrieve junction identifier 
                    which(allSamples %in% samples))
    }
    parsed <- parseValidFile(path, format, select=select)
    
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

downloadGtexFiles <- function(filepath, dataTypes) {
    
    link <- paste0("https://storage.googleapis.com/gtex_analysis_v7/", 
                   rep(c("annotations/", "rna_seq_data/"), each=2),
                   names(filepath))
    toDownload <- !file.exists(gsub("\\.gz$", "", filepath))
    if (sum(toDownload) > 0) {
        updateProgress("Downoading data...", divisions=sum(toDownload))
        for (i in which(toDownload)) {
            updateProgress("Downloading file", detail=dataTypes[i])
            download.file(link[i], filepath[i])
        }
    }
    
    toDecompress <- grepl("\\.gz$", filepath) & file.exists(filepath)
    if (sum(toDecompress) > 0) {
        updateProgress("Extracting files...", divisions=sum(toDecompress))
        for (each in which(toDecompress)) {
            updateProgress("Extracting file", detail=dataTypes[each])
            gunzip(filepath[each])
        }
    }
}

#' Load GTEx data
#'
#' @param data Character: data types to load (see \code{getGtexDataTypes})
#' @param folder Character: folder containing data
#' @param tissue Character: tissues to load (if \code{NULL}, load all); tissue
#' selection may speed up data loading
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom R.utils gunzip
#' 
#' @family functions associated with GTEx data retrieval
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
#' }
loadGtexData <- function(folder=getDownloadsFolder(), data=getGtexDataTypes(), 
                         tissue=NULL) {
    if (is.null(data))   stop("Argument 'data' cannot be NULL.")
    if (is.null(folder)) stop("Argument 'folder' cannot be NULL.")
    
    files <- c(
        "GTEx_v7_Annotations_SampleAttributesDS.txt",
        "GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
        "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz",
        "GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz")
    names(files) <- getGtexDataTypes()
    files <- files[data]
    filepath <- file.path(folder, files)
    names(filepath) <- files
    
    downloadGtexFiles(filepath, data)
    
    updateProgress("Loading files...", divisions=length(data))
    
    loadThisGtexFile <- function(path, pattern, samples=NULL) {
        name <- ifelse(is.character(path), basename(path), path$name)
        updateProgress("Processing file", detail=name)
        loaded <- loadGtexFile(path, pattern, samples)
        return(loaded)
    }
    
    filepath <- gsub("\\.gz$", "", filepath)
    sampleMetadata <- grep("Sample",     filepath, value=TRUE)
    clinical       <- grep("Subject",    filepath, value=TRUE)
    junctionQuant  <- grep("junctions",  filepath, value=TRUE)
    geneExpr       <- grep("gene_reads", filepath, value=TRUE)
    
    loaded <- list()
    samples <- NULL
    if (length(sampleMetadata) > 0) {
        sampleAttrs <- loadThisGtexFile(sampleMetadata, "Sample")
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
    
    if (length(clinical) > 0)
        loaded[[2]] <- loadThisGtexFile(clinical, "Subject", samples)
    
    if (length(junctionQuant) > 0)
        loaded[[3]] <- loadThisGtexFile(junctionQuant, "junction", samples)
    
    if (length(geneExpr) > 0)
        loaded[[4]] <- loadThisGtexFile(geneExpr, "gene", samples)
    
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    
    gtex <- setNames(list(loaded), "GTEx")
    gtex <- processDatasetNames(gtex)
    closeProgress()
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
    dataTypes  <- input$dataTypes
    folder     <- input$folder
    tissue     <- input$tissues
    
    time <- startProcess("load")
    data <- loadGtexData(dataTypes, folder, tissue)
    
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
    observeEvent(input$load, {
        if (!is.null(getData()))
            loadedDataModal(session, "modal", "replace", "append")
        else
            loadGtexDataShiny(session, input)
    })
    
    # Select available tissues from GTEx
    showAvailableTissues <- reactive({
        folder         <- input$folder
        progressBar    <- "loadingAvailableTissues"
        tissueSelect   <- "tissues"
        
        fadeIn  <- function(id, ...) show(id, anim=TRUE, ...)
        fadeOut <- function(id, ...) hide(id, anim=TRUE, ...)
        
        fadeOut(progressBar)
        fadeOut(tissueSelect)
        tissues <- tryCatch(
          getGtexTissues(folder), error=return, warning=return)
        fadeIn(tissueSelect, animType="fade")
        tissues <- c(tissues, "Select available tissues"="")
        
        updateSelectizeInput(session, "tissues", choices=tissues)
    })
    
    observe({
        if (!identical(input$filterCollapse, "Load by tissue")) return(NULL)
        showAvailableTissues()
    })
    
    # Replace or append data to existing data
    observeEvent(input$replace, loadGtexDataShiny(session, input, replace=TRUE))
    observeEvent(input$append, loadGtexDataShiny(session, input, replace=FALSE))
}

attr(gtexDataUI, "loader") <- "data"
attr(gtexDataServer, "loader") <- "data"