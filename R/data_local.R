#' Interface to load local data
#'
#' @importFrom shiny textInput
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' 
#' @param id Character: namespace identifier
#' @param panel Function to deal with the interface
#' 
#' @return NULL (this function is used to modify the Shiny session's state) 
localDataUI <- function(id, panel) {
    ns <- NS(id)
    
    addLocalFile <- tagList(
        uiOutput(ns("localDataModal")),
        uiOutput(ns("pathSuggestions")),
        textAreaInput(ns("localFolder"), "Folder where data is stored",
                      value=getDownloadsFolder(), 
                      placeholder="Insert local folder"),
        textInput(ns("localCategory"), label="Data category name"),
        selectizeInput(ns("localIgnore"), "Files/directories to ignore",
                       choices=getFirehoseDataTypes(),
                       multiple=TRUE, options=list(
                           # Allow to add new items
                           create=TRUE, createOnBlur=TRUE,
                           placeholder="Input files to exclude")),
        processButton(ns("acceptFile"), "Load files"))
    
    panel(
        style="info", title=list(icon("plus-circle"), "Load local files"),
        value="Load local files", addLocalFile)
}

#' Load local files
#' 
#' @param folder Character: path to folder containing files of interest
#' @param name Character: name of the category containing all loaded datasets
#' @param ignore Character: skip folders and filenames that match the expression
#' @param progress Function to keep track of the progress
#' 
#' @importFrom stats setNames
#' 
#' @return List of data frames from valid files
#' @export
#' 
#' @examples
#' \dontrun{
#' folder <- "~/Downloads/ACC 2016"
#' data <- loadLocalFiles(folder)
#' 
#' ignore <- c(".aux.", ".mage-tab.", "junction quantification")
#' loadLocalFiles(folder, ignore)
#' }
loadLocalFiles <- function(folder, ignore=c(".aux.", ".mage-tab."), name="Data",
                           progress=echoProgress) {
    # Get all files in the specified directory and subdirectories
    files <- list.files(folder, recursive=TRUE, full.names=TRUE)
    
    # Exclude undesired subdirectories or files
    files <- files[!dir.exists(files)]
    ignore <- paste(ignore, collapse = "|")
    if (ignore != "") files <- files[!grepl(ignore, files)]
    
    progress("Searching inside the folder...", divisions=length(files))
    
    loaded <- list()
    formats <- loadFileFormats()
    for (each in seq_along(files)) {
        progress("Processing file", detail = basename(files[each]))
        loaded[[each]] <- parseValidFile(files[each], formats)
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    
    data <- setNames(list(loaded), name)
    data <- processDatasetNames(data)
    return(data)
}

#' Load local files
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param replace Boolean: replace loaded data? TRUE by default
#' 
#' @importFrom shinyjs disable enable
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
setLocalData <- function(input, output, session, replace=TRUE) {
    time <- startProcess("acceptFile")
    
    folder <- input$localFolder
    category <- input$localCategory
    ignore <- c(".aux.", ".mage-tab.", input$localIgnore)
    
    # Load valid local files
    progress <- updateProgress
    data <- loadLocalFiles(folder, name=category, ignore, progress)
    
    if (!is.null(data)) {
        if(replace) {
            setData(data)
        } else {
            data <- processDatasetNames(c(getData(), data))
            setData(data)
        }
    }
    endProcess("acceptFile", time)
}

#' Server logic to load local data
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny updateTextInput
#' @return NULL (this function is used to modify the Shiny session's state)
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
                           "The path is directing to a file, but only folders",
                           "are accepted. Please, insert the path to a folder.",
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