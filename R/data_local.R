#' Get Firebrowse choices of data types
#' @importFrom R.utils capitalize
#' 
#' @return Named character vector
getFirebrowseDataChoices <- function() {
    choices <- c(paste0(c("junction", "exon"),
                        "_quantification"), "Preprocess",
                 paste0("RSEM_", c("isoforms", "genes")),
                 paste0(c("junction", "gene", "exon"),
                        "_expression"), "genes_normalized")
    names(choices) <- capitalize(gsub("_", " ", choices))
    return(choices)
}

#' Creates a UI set with options to add a file from the local storage
#' 
#' @param ns Namespace function
#' 
#' @importFrom shiny tagList uiOutput textInput selectizeInput actionButton
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function(ns) {
    tagList(
        uiOutput(ns("localDataModal")),
        uiOutput(ns("pathSuggestions")),
        textAreaInput(ns("localFolder"), "Folder where data is stored",
                      value="~/Downloads/", placeholder="Insert local folder"),
        textInput(ns("localCategory"), label="Data category name"),
        selectizeInput(ns("localIgnore"), "Files/directories to ignore",
                       choices=getFirebrowseDataChoices(),
                       selected=c("RSEM_isoforms", "exon_quantification"),
                       multiple=TRUE, options=list(
                           # Allow to add new items
                           create=TRUE, createOnBlur=TRUE,
                           placeholder="Input files to exclude")),
        actionButton(ns("acceptFile"), class="btn-primary", "Load files")
    ) # end of list
}

#' @importFrom shinyBS bsCollapse bsCollapsePanel
localDataUI <- function(id, panel) {
    ns <- NS(id)
    
    panel(
        style="info", title=list(icon("plus-circle"), "Add local files"),
        value="Add local files", addLocalFile(ns))
}

#' Load local files
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param replace Boolean: replace loaded data? TRUE by default
#' 
#' @importFrom shinyjs disable enable
setLocalData <- function(input, output, session, replace=TRUE) {
    disable("acceptFile")
    
    folder <- input$localFolder
    category <- input$localCategory
    ignore <- c(".aux.", ".mage-tab.", input$localIgnore)
    
    # Get all files in the specified directory and subdirectories
    files <- list.files(folder, recursive = TRUE, full.names = TRUE)
    
    # Exclude undesired subdirectories or files
    files <- files[!dir.exists(files)]
    ignore <- paste(ignore, collapse = "|")
    if (ignore != "") files <- files[!grepl(ignore, files)]
    
    startProgress("Searching inside the folder...", divisions=length(files))
    
    # Try to load files and remove those with 0 rows
    loaded <- list()
    formats <- loadFileFormats()
    for (each in seq_along(files)) {
        updateProgress("Processing file", detail = basename(files[each]))
        loaded[[each]] <- parseValidFile(files[each], formats)
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    
    data <- setNames(list(loaded), category)
    if (!is.null(data)) {
        data <- processDatasetNames(data)
        if(replace)
            setData(data)
        else
            setData(c(getData(), data))
    }
    
    closeProgress()
    enable("acceptFile")
}

#' Server logic to load local data
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny updateTextInput
localDataServer <- function(input, output, session) {
    ns <- session$ns
    
    # # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))
    
    # Update available clinical data attributes to use in a formula
    output$pathSuggestions <- renderUI({
        checkInside <- function(path, showFiles=FALSE) {
            if (substr(path, nchar(path), nchar(path)) == "/") {
                content <- list.files(path, full.names=TRUE)
            } else {
                content <- list.files(dirname(path), full.names=TRUE)
            }
            
            # Show only directories if showFiles is FALSE
            if (!showFiles) content <- content[dir.exists(content)]
            return(basename(content))
        }
        
        textSuggestions(ns("localFolder"), checkInside(input$localFolder),
                        char=.Platform$file.sep)
    })
    
    # If data is loaded, let user replace or append to loaded data
    observeEvent(input$acceptFile, {
        folder <- input$localFolder
        if (!dir.exists(folder)) {
            if (file.exists(folder)){
                # Asking for a folder, not a file
                errorModal(session, "Folder not found",
                           "The path is directing to a file, but only folders are",
                           "accepted. Please, insert the path to a folder.",
                           modalId="localDataModal")
            } else {
                # Folder not found
                errorModal(session, "Folder not found",
                           "Check if path is correct.",
                           modalId="localDataModal")
            }
            enable("acceptFile")
            return(NULL)
        }
        
        if (!is.null(getData()))
            loadedDataModal(session,
                            "localDataModal",
                            "localReplace",
                            "localAppend")
        else
            setLocalData(input, output, session)
    })
    
    # Load data when the user presses to replace data
    observeEvent(input$localReplace,
                 setLocalData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$localAppend,
                 setLocalData(input, output, session, replace=FALSE))
    
    # Update category name input based on given folder
    observe({
        folder <- input$localFolder
        updateTextInput(session, "localCategory", value=basename(folder))
    })
}

attr(localDataUI, "loader") <- "data"
attr(localDataServer, "loader") <- "data"