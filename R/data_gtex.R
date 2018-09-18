#' @rdname appUI
#' 
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom shiny helpText
#' @importFrom shinyjs hidden
gtexDataUI <- function(id, panel) {
    ns <- NS(id)
    
    panel(style="info", title=list(icon("plus-circle"), "Load GTEx files"),
          value="Load GTEx files",
          uiOutput(ns("modal")),
          helpText("Please download files from the",
                   a(href="http://www.gtexportal.org", target="_blank",
                     "GTEx Data Portal"), "and load them here."),
          fileBrowserInput(ns("sampleInfo"), "Sample attributes (TXT file)",
                           placeholder="No file selected"),
          fileBrowserInput(ns("subjectInfo"), "Subject phenotypes (TXT file)",
                           placeholder="No file selected"),
          fileBrowserInput(ns("junctionQuant"), "Junction read counts",
                           placeholder="No file selected"),
          fileBrowserInput(ns("geneExpr"), "Gene expression", 
                           placeholder="No file selected"),
          bsCollapse(
              id=ns("filterCollapse"),
              bsCollapsePanel(
                  title=tagList(icon("filter"), "Filter tissue(s) to load"), 
                  value="Load by tissue",
                  div(id=ns("loadingAvailableTissues"), class="progress",
                      div(class="progress-bar progress-bar-striped active",
                          role="progressbar", style="width:100%",
                          "Loading tissues from sample attributes")),
                  hidden(errorDialog(
                      paste(
                          "Please select a file containing GTEx sample",
                          "attributes. This file is required to obtain",
                          "available tissues in GTEx."), 
                      id=ns("missingData"), style="margin: 10px;")),
                  hidden(warningDialog(
                      paste(
                          "An issue occurred while reading the file",
                          "containing sample attributes. Please confirm",
                          "if the input file is the correct one."), 
                      id=ns("fileReadWarning"), style="margin: 10px;")),
                  hidden(
                      selectizeInput(ns("tissues"), label=NULL, width="100%",
                                     choices=c("Select one or more tissues"=""),
                                     multiple=TRUE)))),
          processButton(ns("load"), "Load data"))
}

#' Get GTEx tissues from given GTEx sample attributes
#'
#' @param sampleMetadata Character: path to sample attributes
#'
#' @return Character: available tissues
#' @export
getGtexTissues <- function(sampleMetadata) {
    tissueCol <- "Tissue Type (area of retrieval)"
    metadata  <- loadGtexFile(sampleMetadata, "Sample")
    freq      <- table(metadata[[tissueCol]])
    
    tissues        <- names(freq)
    names(tissues) <- sprintf("%s (%s samples)", names(freq), as.vector(freq))
    if (!is.na(match("", tissues)))
        return(tissues[-match("", tissues)])
    else
        return(tissues)
}

#' @rdname getGtexTissues
getGTExTissues <- getGtexTissues

#' Load GTEx file
#'
#' @param path Character: path to file
#' @param pattern Character: pattern of the format type to load file
#' @param samples Character: samples to filter datasets
#'
#' @return Loaded file as a data frame
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
        } else if (pattern=="Subject") {
            # Retrieve patients for which samples are available
            patients <- getPatientFromSample(samples, rownames(parsed))
            patients <- patients[!is.na(patients)]
            patients <- sort(unique(patients))
            parsed <- parsed[patients, ]
        }
    }
    return(parsed)
}

#' @rdname loadGtexFile
loadGTExFile <- loadGtexFile

#' Load GTEx data
#'
#' @param clinical Character: path to subject information (the TXT file)
#' @param sampleMetadata Character: path to sample metadata (the TXT file)
#' @param geneExpr Character: path to gene read counts, RPKMs or TPMs
#' @param junctionQuant Character: path to junction quantification
#' @param tissue Character: tissue(s) of interest when loading data (all tissues
#' are loaded by default); if only some tissue(s) are of interest, this may
#' speed up loading the data
#'
#' @importFrom tools file_path_sans_ext
#'
#' @return List with loaded data
#' @export
loadGtexData <- function(clinical=NULL, sampleMetadata=NULL, junctionQuant=NULL,
                         geneExpr=NULL, tissue=NULL) {
    if (is.null(clinical) && is.null(sampleMetadata) && is.null(junctionQuant) 
        && is.null(geneExpr))
        stop("No input data was given.")
    
    loaded <- list()
    validFiles <- !vapply(
        list(clinical, sampleMetadata, junctionQuant, geneExpr),
        is.null, logical(1))
    updateProgress("Loading files...", divisions=sum(validFiles))
    
    loadThisGtexFile <- function(path, pattern, samples=NULL) {
        name <- ifelse(is.character(path), basename(path), path$name)
        updateProgress("Processing file", detail=name)
        loaded <- loadGtexFile(path, pattern, samples)
        return(loaded)
    }
    
    if (!is.null(tissue) && is.null(sampleMetadata))
        stop("Filtering by tissue requires sample metadata as input.")
    
    samples <- NULL
    if (!is.null(sampleMetadata)) {
        sampleAttrs <- loadThisGtexFile(sampleMetadata, "Sample")
        if (!is.null(tissue)) {
            if (is.null(sampleAttrs)) {
                stop("No GTEx samples were found in the given file. Confirm ",
                     "if this is the correct file.")
            }
            
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
    
    if (!is.null(clinical))
        loaded[[2]] <- loadThisGtexFile(clinical, "Subject", samples)
    
    if (!is.null(junctionQuant))
        loaded[[3]] <- loadThisGtexFile(junctionQuant, "junction", samples)
    
    if (!is.null(geneExpr))
        loaded[[4]] <- loadThisGtexFile(geneExpr, "gene", samples)
    
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    
    data <- setNames(list(loaded), "GTEx")
    data <- processDatasetNames(data)
    return(data)
}

#' Shiny wrapper to load GTEx data
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param replace Boolean: replace loaded data? TRUE by default
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
loadGtexDataShiny <- function(session, input, replace=TRUE) {
    tissue <- input$tissues
    
    subjectInfo <- input$subjectInfo
    if (is.null(subjectInfo) || identical(subjectInfo, "")) 
        subjectInfo <- NULL
    
    sampleInfo <- input$sampleInfo
    if (is.null(sampleInfo) || identical(sampleInfo, "")) 
        sampleInfo <- NULL
    
    junctionQuant <- input$junctionQuant
    if (is.null(junctionQuant) || identical(junctionQuant, ""))
        junctionQuant <- NULL
    
    geneExpr <- input$geneExpr
    if (is.null(geneExpr) || identical(geneExpr, ""))
        geneExpr <- NULL
    
    files <- c("Sample attributes"=sampleInfo, "Subject phenotypes"=subjectInfo,
               "Junction read counts"=junctionQuant, "Gene expression"=geneExpr)
    
    if (all(is.null(files))) {
        errorModal(session, "No file provided",
                   "Please input at least one GTEx file.", modalId="modal",
                   caller="Load GTEx data")
    } else if (any(!isFile(files))) {
        formatFileInfo <- function(item, files) {
            filepath <- prepareWordBreak(files[[item]])
            tagList(tags$b(names(files[item])), tags$br(),
                    tags$kbd(filepath), tags$br(), tags$br())
        }
        
        nonExisting   <- files[!isFile(files)]
        filesNotFound <- do.call(
            tagList, lapply(names(nonExisting), formatFileInfo, nonExisting))
        
        errorModal(session, "Files not found",
                   "The following files were not found:", tags$br(), tags$br(),
                   filesNotFound, modalId="modal", caller="Load GTEx data")
    } else {
        time <- startProcess("load")
        data <- loadGtexData(clinical=subjectInfo, sampleMetadata=sampleInfo,
                             junctionQuant=junctionQuant, geneExpr=geneExpr,
                             tissue=tissue)
        
        if (!is.null(data)) {
            if (!replace) data <- c(getData(), data)
            data <- processDatasetNames(data)
            setData(data)
        }
        endProcess("load", time)
    }
}

#' @rdname appServer
#' @importFrom shinyjs show hide
gtexDataServer <- function(input, output, session) {
    prepareFileBrowser(session, input, "sampleInfo")
    prepareFileBrowser(session, input, "subjectInfo")
    prepareFileBrowser(session, input, "junctionQuant")
    prepareFileBrowser(session, input, "geneExpr")
    
    observeEvent(input$load, {
        if (!is.null(getData()))
            loadedDataModal(session, "modal", "replace", "append")
        else
            loadGtexDataShiny(session, input)
    })
    
    # Select available tissues from GTEx
    showAvailableTissues <- reactive({
        sampleMetadata <- input$sampleInfo
        progressBar    <- "loadingAvailableTissues"
        tissueSelect   <- "tissues"
        alert          <- "missingData"
        warning        <- "fileReadWarning"
        
        fadeIn  <- function(id, ...) show(id, anim=TRUE, ...)
        fadeOut <- function(id, ...) hide(id, anim=TRUE, ...)
        
        fadeOut(progressBar)
        if (!is.null(sampleMetadata) && !identical(sampleMetadata, "")) {
            fadeOut(tissueSelect)
            
            tissues <- tryCatch(getGtexTissues(sampleMetadata), error=return,
                                warning=return)
            
            if (is.null(tissues) || is(tissues, "error") ||
                is(tissues, "warning")) {
                # Warn the user when having an issue reading the file
                fadeOut(tissueSelect, animType="fade")
                fadeIn(warning)
                tissues <- character(0)
            } else {
                fadeIn(tissueSelect, animType="fade")
                fadeOut(warning)
                tissues <- c(tissues, "Select available tissues"="")
            }
            fadeOut(alert)
        } else {
            fadeOut(tissueSelect, animType="fade")
            fadeIn(alert)
            fadeOut(warning)
            tissues <- character(0)
        }
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