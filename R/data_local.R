#' @importFrom R.utils capitalize
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
#' @return A UI set that can be added to a UI definition
addLocalFile <- function(ns) {
    tagList(
        uiOutput(ns("localDataModal")),
        uiOutput(ns("pathSuggestions")),
        textAreaInput(ns("localFolder"), "Folder where data is stored",
                      value="~/Downloads/", placeholder="Insert local folder"),
        textInput(ns("localCategory"), label="Category name",
                  value="Adenoid cystic carcinoma (ACC) 2016"),
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
    
    sub <- dir(folder, full.names=TRUE)[dir.exists(
        dir(folder, full.names=TRUE))]
    
    startProgress("Searching inside the folder...",
                  divisions=1 + length(sub))
    loaded <- loadFirehoseFolders(sub, ignore, updateProgress)
    data <- setNames(list(loaded), category)
    
    if (!is.null(data)) {
        if(replace)
            setData(data)
        else
            setData(c(getData(), data))
    }
    
    closeProgress()
    enable("acceptFile")
}

localDataServer <- function(input, output, session, active) {
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
}

attr(localDataUI, "loader") <- "data"
attr(localDataServer, "loader") <- "data"