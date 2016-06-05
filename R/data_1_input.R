name <- "Input"

#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function() {
    list(
        textAreaInput(id("localFolder"), "Folder where data is stored",
                      value = "~/Downloads/", placeholder = "Insert local folder"),
        textInput(id("localCategory"), label = "Category name", 
                  value = "Adenoid cystic carcinoma (ACC) 2016"),
        selectizeInput(id("localIgnore"), "Files/directories to ignore",
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
        actionButton(id("acceptFile"), class = "btn-primary", "Load files")
    ) # end of list
}

#' Creates a UI set with options to add data from TCGA/Firehose
#' @importFrom shinyBS bsTooltip
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function() {
    if (isFirehoseUp()) {
        cohorts <- getFirehoseCohorts()
        acronyms <- names(cohorts)
        names(acronyms) <- sprintf("%s (%s)", cohorts, names(cohorts))
        
        dates <- as.character(getFirehoseDates())
        
        list(
            selectizeInput(id("firehoseCohort"), "Cohort", acronyms,
                           multiple = TRUE, selected = c("ACC", "BLCA"),
                           options = list(placeholder = "Select cohort(s)")),
            selectizeInput(id("firehoseDate"), "Date", dates, multiple = TRUE,
                           selected = dates[1], options = list(
                               placeholder = "Select sample date")),
            selectizeInput(id("dataType"), "Data type",
                           c("Clinical", "mRNASeq"), 
                           multiple = TRUE, selected = "Clinical",
                           options = list(
                               placeholder = "Select data types")),
            selectizeInput(id("firehoseIgnore"), "Files/archives to ignore",
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
            bsTooltip(id("firehoseIgnore"), placement = "right",
                      options = list(container = "body"),
                      paste("Files which contain these terms won\\'t be",
                            "either downloaded or loaded.")),
            textAreaInput(id("dataFolder"), "Folder to store the data",
                          value = "~/Downloads/",
                          placeholder = "Insert data folder"),
            bsTooltip(id("dataFolder"), placement = "right",
                      options = list(container = "body"),
                      "Data not available in this folder will be downloaded."),
            actionButton(class = "btn-primary", type = "button",
                         id("getFirehoseData"), "Get data"))
    } else {
        list(icon("exclamation-circle"),
             "Firehose seems to be offline at the moment.")
    }
}

#' Create a modal warning the user of loaded data
loadedDataModal <- function(modalId, replaceButtonId, keepButtonId) {
    bsModal2(modalId,
             div(icon("exclamation-triangle"), "Data already loaded"),
             NULL, size = "small", style = "warning",
             "Would you like to", tags$b("replace"), "the loaded data or",
             tags$b("keep"), "both the previous and new data?",
             footer = tagList(
                 actionButton(keepButtonId, "data-dismiss"="modal",
                              label="Keep both"),
                 actionButton(replaceButtonId, class = "btn-warning",
                              "data-dismiss"="modal", label="Replace")))
}

#' @importFrom shinyBS bsCollapse bsCollapsePanel bsAlert
ui <- function() {
    list(
        # TODO(NunoA): Show alerts from renderUI
        bsAlert(anchorId = id("alert2")),
        loadedDataModal(id("localDataModal"), id("localReplace"),
                        id("localAppend")),
        loadedDataModal(id("firebrowseDataModal"), id("firebrowseReplace"),
                        id("firebrowseAppend")),
        uiOutput(id("pathAutocomplete")),
        uiOutput("iframeDownload"),
        bsCollapse(
            id = id("addData"),
            open = "Add TCGA/Firehose data",
            bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"), "Add local files"),
                value = "Add local files",
                addLocalFile()),
            bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"),
                             "Add TCGA/Firehose data"),
                value = "Add TCGA/Firehose data",
                addTCGAdata()))
    )
}

#' Load local files
#' @param replace Boolean: replace loaded data? TRUE by default
#' @importFrom shinyjs disable enable
setLocalData <- function(input, output, session, replace=TRUE) {
    disable(id("acceptFile"))
    
    folder <- input[[id("localFolder")]]
    category <- input[[id("localCategory")]]
    ignore <- input[[id("localIgnore")]]
    
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
    enable(id("acceptFile"))
}

#' Set data from Firehose
#' @param replace Boolean: replace loaded data? TRUE by default
#' @importFrom shinyjs disable enable 
setFirehoseData <- function(input, output, session, replace=TRUE) {
    disable(id("getFirehoseData"))
    
    # Load data from Firehose
    data <- loadFirehoseData(
        folder = input[[id("dataFolder")]],
        cohort = input[[id("firehoseCohort")]],
        date = gsub("-", "_", input[[id("firehoseDate")]]),
        data_type = input[[id("dataType")]],
        exclude = input[[id("firehoseIgnore")]],
        progress = updateProgress,
        output = output)
    
    if (!is.null(data)) {
        if(replace)
            setData(data)
        else
            setData(c(getData(), data))
    }
    
    closeProgress()
    enable(id("getFirehoseData"))
}

#' @importFrom shinyBS toggleModal
server <- function(input, output, session) {
    # The button is only enabled if it meets the conditions that follow
    # observe(toggleState(id("acceptFile"), input[[id("species")]] != ""))
    
    # Update available clinical data attributes to use in a formula
    output[[id("pathAutocomplete")]] <- renderUI({
        checkInside <- function(path) {
            if (substr(path, nchar(path), nchar(path)) == "/") {
                list.files(path)
            } else {
                list.files(dirname(path))
            }
        }
        
        tagList(
            textComplete(id("localFolder"), 
                         checkInside(input[[id("localFolder")]]),
                         char=.Platform$file.sep),
            textComplete(id("dataFolder"), 
                         checkInside(input[[id("dataFolder")]]),
                         char=.Platform$file.sep)
        )
    })
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input[[id("acceptFile")]], {
        if (!is.null(getData()))
            toggleModal(session, id("localDataModal"), "open")
        else
            setLocalData(input, output, session)
    })
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("localReplace")]],
                 setLocalData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input[[id("localAppend")]],
                 setLocalData(input, output, session, replace=FALSE))
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input[[id("getFirehoseData")]], {
        if (!is.null(getData()))
            toggleModal(session, id("firebrowseDataModal"), "open")
        else
            setFirehoseData(input, output, session)
    })
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("firebrowseReplace")]],
                 setFirehoseData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input[[id("firebrowseAppend")]],
                 setFirehoseData(input, output, session, replace=FALSE))
}