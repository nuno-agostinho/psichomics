## TODO(NunoA): correct progress bar when loading local files

name <- "Input"

#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function() {
    list(
        textAreaInput(id("localFolder"), "Folder where data is stored",
                      value = "~/Downloads"),
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
#' 
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function() {
    if (isFirehoseUp()) {
        cohorts <- getFirehoseCohorts()
        names(cohorts) <- sprintf("%s (%s)", names(cohorts), cohorts)
        dates <- as.character(getFirehoseDates())
        
        list(
            selectizeInput(id("firehoseCohort"), "Cohort", cohorts,
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
                          value = "~/Downloads",
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
loadedDataModal <- function(modalId, replaceButtonId) {
    bsModal2(modalId,
             div(icon("exclamation-triangle"), "Data already loaded"),
             NULL, size = "small", style = "warning",
             "Would you like to replace the loaded data?",
             footer = list(actionButton(replaceButtonId,
                                        class = "btn-warning",
                                        "data-dismiss"="modal", 
                                        label = "Replace")))
}
    
ui <- function() {
    list(
        # TODO(NunoA): Show alerts from renderUI
        bsAlert(anchorId = id("alert2")),
        loadedDataModal(id("localDataModal"), id("localReplace")),
        loadedDataModal(id("firebrowseDataModal"), id("firebrowseReplace")),
        uiOutput(id("pathAutocomplete")),
        uiOutput("iframeDownload"),
        shinyBS::bsCollapse(
            id = id("addData"),
            open = "Add TCGA/Firehose data",
            shinyBS::bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"), "Add local files"),
                value = "Add local files",
                addLocalFile()),
            shinyBS::bsCollapsePanel(
                style = "info",
                title = list(icon("plus-circle"),
                             "Add TCGA/Firehose data"),
                value = "Add TCGA/Firehose data",
                addTCGAdata()))
    )
}

server <- function(input, output, session) {
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
            loadLocalData()
    })
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("localReplace")]], loadLocalData())
    
    # Load local files
    loadLocalData <- function() {
        shinyjs::disable(id("acceptFile"))
        
        folder <- input[[id("localFolder")]]
        category <- input[[id("localCategory")]]
        ignore <- input[[id("localIgnore")]]
        
        sub <- dir(folder, full.names = T)[dir.exists(dir(folder, full.names = T))]
        
        startProgress("Searching inside the folder...", divisions = 1 + length(sub))
        loaded <- loadFirehoseFolders(sub, ignore, updateProgress)
        data <- setNames(list(loaded), category)
        setData(data)
        
        closeProgress()
        shinyjs::enable(id("acceptFile"))
    }
    
    # The button is only enabled if it meets the conditions that follow
    # observe(toggleState(id("acceptFile"), input[[id("species")]] != ""))
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input[[id("getFirehoseData")]], {
        if (!is.null(getData()))
            toggleModal(session, id("firebrowseDataModal"), "open")
        else
            loadAllData()
    })
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("firebrowseReplace")]], loadAllData())
    
    # Load Firehose data
    loadAllData <- function() {
        shinyjs::disable(id("getFirehoseData"))
        
        # Load data from Firehose
        data <- loadFirehoseData(
            folder = input[[id("dataFolder")]],
            cohort = input[[id("firehoseCohort")]],
            date = gsub("-", "_", input[[id("firehoseDate")]]),
            data_type = input[[id("dataType")]],
            exclude = input[[id("firehoseIgnore")]],
            progress = updateProgress,
            output = output)
        
        if (!is.null(data))
            setData(data)
        
        closeProgress()
        shinyjs::enable(id("getFirehoseData"))
    }
}