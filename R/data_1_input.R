#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function(ns) {
    tagList(
        textAreaInput(ns("localFolder"), "Folder where data is stored",
                      value = "~/Downloads/", placeholder = "Insert local folder"),
        textInput(ns("localCategory"), label = "Category name",
                  value = "Adenoid cystic carcinoma (ACC) 2016"),
        selectizeInput(ns("localIgnore"), "Files/directories to ignore",
                       choices = c(".aux.", ".mage-tab.",
                                   paste0(c("junction", "exon"),
                                          "_quantification"), "Preprocess",
                                   paste0("RSEM_", c("isoforms", "genes")),
                                   paste0(c("junction", "gene", "exon"),
                                          "_expression"), "genes_normalized"),
                       selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                    "exon_quantification"),
                       multiple = TRUE, options = list(
                           # Allow to add new items
                           create = TRUE, createOnBlur=TRUE,
                           placeholder = "Input files to exclude")),
        actionButton(ns("acceptFile"), class = "btn-primary", "Load files")
    ) # end of list
}

#' Creates a UI set with options to add data from TCGA/Firehose
#' @importFrom shinyBS bsTooltip
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function(ns) {
    if (isFirehoseUp()) {
        cohorts <- getFirehoseCohorts()
        acronyms <- names(cohorts)
        names(acronyms) <- sprintf("%s (%s)", cohorts, names(cohorts))
        
        dates <- as.character(getFirehoseDates())
        
        tagList(
            selectizeInput(ns("firehoseCohort"), "Cohort", acronyms,
                           multiple = TRUE, selected = c("ACC", "BLCA"),
                           options = list(placeholder = "Select cohort(s)")),
            selectizeInput(ns("firehoseDate"), "Date", dates, multiple = TRUE,
                           selected = dates[1], options = list(
                               placeholder = "Select sample date")),
            selectizeInput(ns("dataType"), "Data type",
                           c("Clinical", "mRNASeq"), 
                           multiple = TRUE, selected = "Clinical",
                           options = list(
                               placeholder = "Select data types")),
            selectizeInput(ns("firehoseIgnore"), "Files/archives to ignore",
                           choices = c(".aux.", ".mage-tab.",
                                       paste0(c("junction", "exon"),
                                              "_quantification"), "Preprocess",
                                       paste0("RSEM_", c("isoforms", "genes")),
                                       paste0(c("junction", "gene", "exon"),
                                              "_expression"), "genes_normalized"),
                           selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                        "exon_quantification"),
                           multiple = TRUE, options = list(
                               # Allow to add new items
                               create = TRUE, createOnBlur=TRUE,
                               placeholder = "Input files to exclude")),
            bsTooltip(ns("firehoseIgnore"), placement = "right",
                      options = list(container = "body"),
                      paste("Files which contain these terms won\\'t be",
                            "either downloaded or loaded.")),
            textAreaInput(ns("dataFolder"), "Folder to store the data",
                          value = "~/Downloads/",
                          placeholder = "Insert data folder"),
            bsTooltip(ns("dataFolder"), placement = "right",
                      options = list(container = "body"),
                      "Data not available in this folder will be downloaded."),
            actionButton(class = "btn-primary", type = "button",
                         ns("getFirehoseData"), "Get data"))
    } else {
        list(icon("exclamation-circle"),
             "Firehose seems to be offline at the moment.")
    }
}

#' Create a modal warning the user of already loaded data
#' @param modalId Character: identifier of the modal
#' @param replaceButtonId Character: identifier of the button to replace data
#' @param keepButtonId Character: identifier of the button to append data
loadedDataModal <- function(session, modalId, replaceButtonId, keepButtonId) {
    ns <- session$ns
    warningModal(session, "Data already loaded",
                 "Would you like to", tags$b("replace"), "the loaded data or",
                 tags$b("keep"), "both the previous and new data?",
                 footer = tagList(
                     actionButton(ns(keepButtonId), "data-dismiss"="modal",
                                  label="Keep both"),
                     actionButton(ns(replaceButtonId), class = "btn-warning",
                                  "data-dismiss"="modal", label="Replace")),
                 modalId=modalId)
}

#' @importFrom shinyBS bsCollapse bsCollapsePanel
inputUI <- function(id, tab) {
    ns <- NS(id)
    tab("Input", br(),
        uiOutput(ns("localDataModal")),
        uiOutput(ns("firebrowseDataModal")),
        uiOutput(ns("pathAutocomplete")),
        uiOutput(ns("iframeDownload")),
        bsCollapse(
            id = ns("addData"),
            open = "Add TCGA/Firehose data",
            bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"), "Add local files"),
                value = "Add local files",
                addLocalFile(ns)),
            bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"),
                             "Add TCGA/Firehose data"),
                value = "Add TCGA/Firehose data",
                addTCGAdata(ns)))
    )
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
    ignore <- input$localIgnore
    
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

#' Set data from Firehose
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param replace Boolean: replace loaded data? TRUE by default
#' @importFrom shinyjs disable enable 
setFirehoseData <- function(input, output, session, replace=TRUE) {
    disable("getFirehoseData")
    
    # Load data from Firehose
    data <- loadFirehoseData(
        folder = input$dataFolder,
        cohort = input$firehoseCohort,
        date = gsub("-", "_", input$firehoseDate),
        data_type = input$dataType,
        exclude = input$firehoseIgnore,
        progress = updateProgress,
        output = output)
    
    if (!is.null(data)) {
        if(replace)
            setData(data)
        else
            setData(c(getData(), data))
    }
    
    closeProgress()
    enable("getFirehoseData")
}

inputServer <- function(input, output, session, active) {
    ns <- session$ns
    
    # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))
    
    # Update available clinical data attributes to use in a formula
    output$pathAutocomplete <- renderUI({
        checkInside <- function(path, showFiles=FALSE) {
            if (substr(path, nchar(path), nchar(path)) == "/") {
                content <- list.files(path, full.names = TRUE)
            } else {
                content <- list.files(dirname(path), full.names = TRUE)
            }
            
            # Show only directories if showFiles is FALSE
            if (!showFiles) content <- content[dir.exists(content)]
            return(basename(content))
        }

        tagList(
            textComplete(ns("localFolder"),
                         checkInside(input$localFolder),
                         char=.Platform$file.sep),
            textComplete(ns("dataFolder"),
                         checkInside(input$dataFolder),
                         char=.Platform$file.sep)
        )
    })
    
    # Check if data is already loaded and ask the user if it should be replaced
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
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input$getFirehoseData, {
        if (!is.null(getData()))
            loadedDataModal(session,
                            "firebrowseDataModal",
                            "firebrowseReplace",
                            "firebrowseAppend")
        else
            setFirehoseData(input, output, session)
    })
    
    # Load data when the user presses to replace data
    observeEvent(input$firebrowseReplace,
                 setFirehoseData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$firebrowseAppend,
                 setFirehoseData(input, output, session, replace=FALSE))
}

attr(inputUI, "loader") <- "data"
attr(inputServer, "loader") <- "data"