#' Interface to load GTEx data
#' 
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom shiny helpText fileInput
#' 
#' @param id Character: namespace identifier
#' @param panel Function to deal with the interface
#' 
#' @return NULL (this function is used to modify the Shiny session's state) 
gtexDataUI <- function(id, panel) {
    ns <- NS(id)
    
    panel(style="info", title=list(icon("plus-circle"), "Load GTEx files"),
          value="Load GTEx files",
          uiOutput(ns("modal")),
          helpText("Please, download files from the",
                   a(href="http://www.gtexportal.org", target="_blank",
                     "GTEx Data Portal.")),
          fileInput(ns("inputSampleInfo"),
                    "Choose file with GTEx sample attributes"),
          fileInput(ns("inputSubjectInfo"),
                    "Choose file with GTEx subject phenotypes"),
          fileInput(ns("inputJunctionQuant"), 
                    "Choose file with GTEx junction read counts"),
          processButton(ns("load"), "Load data"))
}

#' Load GTEx data given input
#' 
#' @param input Shiny input
#' @param replace Boolean: replace loaded data? TRUE by default
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
loadGtexData <- function(input, replace=TRUE) {
    time <- startProcess("load")
    # Retrieve formats to parse files
    formats <- loadFileFormats()
    
    loaded <- list()
    index <- 1
    
    startProgress("Loading files...", 3)
    
    loadGtexFile <- function(element, pattern) {
        format <- formats[sapply(formats, function(i)
            grepl(pattern, i$filename, fixed=TRUE) &&
                grepl("GTEx", i$filename, fixed=TRUE))]
        path <- input[[element]]
        if ( !is.null(path) ) {
            updateProgress("Processing file", detail=path$name)
            parsed <- parseValidFile(path$datapath, format)
            return(parsed)
        } else {
            return(NULL)
        }
    }
    
    loaded[[1]] <- loadGtexFile("inputSubjectInfo", "Subject")
    loaded[[2]] <- loadGtexFile("inputSampleInfo", "Sample")
    loaded[[3]] <- loadGtexFile("inputJunctionQuant", "junction")
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    
    data <- setNames(list(loaded), "GTEx")
    data <- processDatasetNames(data)
    
    if (!is.null(data)) {
        if(replace) {
            setData(data)
        } else {
            data <- processDatasetNames(c(getData(), data))
            setData(data)
        }
    }
    endProcess("load", time)
}

#' Server logic to load GTEx data
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
gtexDataServer <- function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$load, {
        if (!is.null(getData()))
            loadedDataModal(session, "modal", "replace", "append")
        else
            loadGtexData(input)
    })
    
    # Replace or append data to existing data
    observeEvent(input$replace, loadGtexData(input, replace=TRUE))
    observeEvent(input$append, loadGtexData(input, replace=FALSE))
}

attr(gtexDataUI, "loader") <- "data"
attr(gtexDataServer, "loader") <- "data"