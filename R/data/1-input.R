name <- "Input"

#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function() {
    list(
        textInput(id("localFolder"), "Folder where data is stored",
                  value = "~/Downloads"),
        textInput(id("localCategory"), label = "Category name", 
                  value = "Adenoid cystic carcinoma (ACC) 2016"),
        selectizeInput(id("localIgnore"), "Files/directories to ignore",
                       choices = c(".aux.", ".mage-tab.", "MANIFEST.txt",
                                   "exon_quantification", "Preprocess",
                                   paste0("RSEM_", c("isoforms", "genes")),
                                   paste0(c("junction", "gene", "exon"),
                                          "_expression"), "genes_normalized"),
                       selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                    "MANIFEST.txt", "exon_quantification"),
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
                           choices = c(".aux.", ".mage-tab.", "MANIFEST.txt", 
                                       "exon_quantification", "Preprocess",
                                       paste0("RSEM_", c("isoforms", "genes")),
                                       paste0(c("junction", "gene", "exon"),
                                              "_expression"), "genes_normalized"),
                           selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                        "MANIFEST.txt", "exon_quantification"),
                           multiple = TRUE, options = list(
                               # Allow to add new items
                               create = TRUE, createOnBlur=TRUE,
                               placeholder = "Input files to exclude")),
            bsTooltip(id("firehoseIgnore"), placement = "right",
                      options = list(container = "body"),
                      paste("Files which contain these terms won\\'t be",
                            "either downloaded or loaded.")),
            textInput(id("dataFolder"), "Folder to store the data",
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

ui <- function() {
    list(
        # TODO(NunoA): Show alerts from renderUI
        bsAlert(anchorId = id("alert2")),
        bsModal2(id("dataReplace"),
                 div(icon("exclamation-triangle"), "Data already loaded"),
                 NULL, size = "small", style = "warning",
                 "Would you like to replace the loaded data?",
                 footer = list(
                     actionButton(id("replace"),
                                  class = "btn-warning",
                                  "data-dismiss"="modal", 
                                  label = "Replace"))),
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
    # Load user files
    observeEvent(input[[id("acceptFile")]], {
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
        
        # error <- function(msg) { print(msg); return(NULL) }
        # if(is.null(input[[id("dataFile")]]))
            # error("No data input selected")
        # if(input[[id("species")]] == "")
            # error("Species field can't be empty")
        
        # inFile <- input$dataFile
        # info <- read.table(inFile$datapath, sep = input$sep,
        #                    header = input$header)
    })
    
    # The button is only enabled if it meets the conditions that follow
    observe(toggleState(id("acceptFile"), input[[id("species")]] != ""))
    
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
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("replace")]], loadAllData())
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input[[id("getFirehoseData")]], {
        if (!is.null(getData()))
            toggleModal(session, id("dataReplace"), "open")
        else
            loadAllData()
    })
}