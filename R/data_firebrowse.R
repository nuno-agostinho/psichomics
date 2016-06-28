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
            uiOutput(ns("firebrowseDataModal")),
            uiOutput(ns("pathAutocomplete")),
            uiOutput(ns("iframeDownload")),
            selectizeInput(ns("firehoseCohort"), "Cohort", acronyms,
                           multiple = TRUE, selected = c("ACC"),
                           options = list(placeholder = "Select cohort(s)")),
            selectizeInput(ns("firehoseDate"), "Date", dates, multiple = TRUE,
                           selected = dates[1], options = list(
                               placeholder = "Select sample date")),
            selectizeInput(ns("firehoseData"), "Data type",
                           c("Clinical", getFirebrowseDataChoices()), 
                           multiple = TRUE, selected = "Clinical",
                           options = list(
                               placeholder = "Select data types")),
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

#' @importFrom shinyBS bsCollapse bsCollapsePanel
firebrowseUI <- function(id, panel) {
    ns <- NS(id)
    
    panel(
        style = "info",
        title = list(icon("plus-circle"), "Add TCGA/Firehose data"),
        value = "Add TCGA/Firehose data", addTCGAdata(ns))
}

#' Set data from Firehose
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param replace Boolean: replace loaded data? TRUE by default
#' @importFrom shinyjs disable enable 
setFirehoseData <- function(input, output, session, replace=TRUE) {
    disable("getFirehoseData")
    
    data <- input$firehoseData
    datasets <- getFirebrowseDataChoices()
    # Data types to load
    data_type <- c(data[!data %in% datasets], "mRNASeq")
    # Datasets to ignore
    ignore <- datasets[!datasets %in% data]
    
    # Load data from Firehose
    data <- loadFirehoseData(
        folder = input$dataFolder,
        cohort = input$firehoseCohort,
        date = gsub("-", "_", input$firehoseDate),
        data_type = data_type,
        exclude = c(".aux.", ".mage-tab.", ignore),
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

firebrowseServer <- function(input, output, session, active) {
    ns <- session$ns
    
    # # The button is only enabled if it meets the conditions that follow
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
        
        textComplete(ns("dataFolder"), checkInside(input$dataFolder),
                     char=.Platform$file.sep)
    })
    
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

attr(firebrowseUI, "loader") <- "data"
attr(firebrowseServer, "loader") <- "data"