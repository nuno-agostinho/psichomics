#' @rdname appUI
#' 
#' @importFrom shiny textInput
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsPopover
localDataUI <- function(id, panel) {
    ns <- NS(id)

    addMultipleFiles <- tagList(
        helpText("All fields below are optional."),
        fileBrowserInput(
            ns("sampleInfo"),
            "File with sample information",
            placeholder="No file selected",
            info=TRUE, infoFUN=bsPopover, infoTitle=paste(
                "File containing sample identifiers as rows and their",
                "attributes as columns."),
            infoContent=paste(
                tags$ul(
                    class="popover-list",
                    tags$li("The first column must contain sample identifiers",
                            "and be named", tags$kbd("Sample ID")),
                    tags$li("Optionally, indicate the subject associated to",
                            "each sample in a column named",
                            tags$kbd("Subject ID"))),
                tags$hr(), helpText("Example:"), tags$table(
                    class="table table-condensed",
                    tags$thead(
                        tableRow("Sample ID", "Type", "Tissue", "Subject",
                                 th=TRUE)),
                    tags$tbody(
                        tableRow("SMP-01", "Tumour", "Lung", "SUBJ-03"),
                        tableRow("SMP-02", "Normal", "Blood", "SUBJ-12"),
                        tableRow("SMP-03", "Normal", "Blood", "SUBJ-25"))))),
        fileBrowserInput(
            ns("subjectInfo"),
            "File with subject information",
            placeholder="No file selected",
            info=TRUE, infoFUN=bsPopover, infoTitle=paste(
                "File containing subject identifiers as rows and their",
                "attributes as columns."),
            infoContent=paste(
                "The first column must contain subject identifiers and be",
                "named", tags$kbd("Subject ID"), tags$hr(),
                helpText("Example:"), tags$table(
                    class="table table-condensed",
                    tags$thead(
                        tableRow("Subject ID", "Age", "Gender", "Race", 
                                 th=TRUE)),
                    tags$tbody(
                        tableRow("SUBJ-01", "4", "Female", "Black"),
                        tableRow("SUBJ-02", "12", "Female", "Black"),
                        tableRow("SUBJ-03", "8", "Female", "Asian"))))),
        geneExprFileInput(ns("geneExpr")),
        fileBrowserInput(
            ns("junctionQuant"),
            "File with exon-exon junction read counts",
            placeholder="No file selected",
            info=TRUE, infoFUN=bsPopover, infoTitle=paste(
                "File containing the read counts of each exon-exon junction",
                "(rows) per sample (columns)."),
            infoContent=paste(
                tags$ul(
                    class="popover-list",
                    tags$li(
                        "The first column must contain junction identifiers",
                        "and be named", tags$kbd("Junction ID")),
                    tags$li(
                        "Only numbers (or X and Y) are extracted from the",
                        "junction identifier. Acceptable junction identifiers",
                        "include: ", tags$kbd("10_18748_21822"), ", ",
                        tags$kbd("chromosome 10 (18748 to 21822)"), " and ", 
                        tags$kbd("chr10:18748-21822")),
                    tags$li(
                        "The strand is optional. If desired, it needs to come",
                        "in the end of the string. For instance,", 
                        tags$kbd("10:3213:9402:+"), "and",
                        tags$kbd("chr10:3213-9402 -"))),
                tags$hr(), helpText("Example:"), tags$table(
                    class="table table-condensed",
                    tags$thead(
                        tableRow("Junction ID", "SMP-18", "SMP-03", th=TRUE)),
                    tags$tbody(
                        tableRow("10:6752-7393", "4", "0"),
                        tableRow("10:18748-21822", "8", "46"),
                        tableRow("10:24257-25325", "83", "65"))))),
        bsCollapse(
            id=ns("ASquantLoadCollapse"),
            bsCollapsePanel(
                title=tagList(icon("plus-circle"),
                              "Load alternative splicing quantification"),
                value="Load alternative splicing quantification",
                ASquantFileInput(ns("ASquant"), ns("customSpecies"),
                                 ns("customAssembly")))),
        textInput(ns("userFilesCategory"), label="Dataset name", width = "100%",
                  value="User dataset", placeholder="Name to identify dataset"),
        processButton(ns("loadMultipleFiles"), "Load files"))
    
    addFolder <- tagList(
        helpText("For your convenience, move all files to a single folder and",
                 "load them by locating that folder in the field below."),
        fileBrowserInput(ns("localFolder"), "Folder where data is stored",
                         placeholder="No folder selected",
                         value=getDownloadsFolder()),
        textInput(ns("localCategory"), label="Dataset name",
                  placeholder="Name to identify dataset", width = "100%"),
        selectizeInput(ns("localIgnore"), "Files/directories to ignore",
                       choices=getFirebrowseDataTypes(), width = "100%",
                       multiple=TRUE, options=list(
                           # Allow to add new items
                           create=TRUE, createOnBlur=TRUE,
                           placeholder="These files will not be loaded")),
        processButton(ns("acceptFile"), "Load files"))
    
    panel(style="info", title=list(icon("plus-circle"), "Load user files"),
          value="Load local files",
          uiOutput(ns("localDataModal")),
          tabsetPanel(
              tabPanel("Stepwise file input", addMultipleFiles),
              tabPanel("Folder input", addFolder)))
}

#' Load local files
#' 
#' @param folder Character: path to folder containing files of interest
#' @param name Character: name of the category containing all loaded datasets
#' @param ignore Character: skip folders and filenames that match the expression
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
loadLocalFiles <- function(folder, ignore=c(".aux.", ".mage-tab."), 
                           name="Data") {
    # Get all files in the specified directory and subdirectories
    files <- list.files(folder, recursive=TRUE, full.names=TRUE)
    
    # Exclude undesired subdirectories or files
    files <- files[!dir.exists(files)]
    ignore <- paste(ignore, collapse = "|")
    if (ignore != "") files <- files[!grepl(ignore, files)]
    
    updateProgress("Searching inside the folder...", divisions=length(files))
    
    loaded <- list()
    formats <- loadFileFormats()
    for (each in seq_along(files)) {
        updateProgress("Processing file", detail = basename(files[each]))
        loadedFile <- tryCatch(parseValidFile(files[each], formats),
                               warning=return, error=return)
        if (!is(loadedFile, "warning") && !is(loadedFile, "error"))
            loaded[[each]] <- loadedFile
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    loaded <- loadTCGAsampleMetadata(loaded)
    
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
    if (identical(category, "")) category <- "User dataset"
    ignore <- c(".aux.", ".mage-tab.", input$localIgnore)
    
    # Load valid local files
    data <- loadLocalFiles(folder, name=category, ignore)
    
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

#' @rdname setLocalData
setMultipleFilesData <- function(input, output, session, replace=TRUE) {
    time <- startProcess("loadMultipleFiles")
    category <- input$userFilesCategory
    if (identical(category, "")) category <- "User dataset"
    
    # Load files
    files <- c("Sample metadata"        =input$sampleInfo,
               "Clinical data"          =input$subjectInfo, 
               "Gene expression"        =input$geneExpr,
               "Junction quantification"=input$junctionQuant,
               "Inclusion levels"       =input$ASquant)
    files <- files[files != ""]
    ASquant <- input$ASquant
    
    # Check if at least one input file was given
    if (length(files) == 0) {
        errorModal(session, "No files selected",
                   "Please provide at least one file.",
                   modalId="localDataModal")
        endProcess("loadMultipleFiles", time, closeProgressBar=FALSE)
        return(NULL)
    }
    
    loaded <- list()
    allFormats <- loadFileFormats()
    updateProgress("Loading files...", divisions=length(files))
    for (each in seq(files)) {
        file <- files[[each]]
        if (file == "") next
        
        # Get appropriate format for each type
        dataType <- names(files[each])
        formats <- allFormats[sapply(allFormats, "[[", "dataType") == dataType]
        
        updateProgress("Processing file", detail=basename(file))
        loadedFile <- tryCatch(parseValidFile(file, formats),
                               warning=return, error=return)
        if (!is(loadedFile, "warning") && !is(loadedFile, "error"))
            loaded[[each]] <- loadedFile
    }
    
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    loaded <- loadTCGAsampleMetadata(loaded)
    
    data <- setNames(list(loaded), category)
    data <- processDatasetNames(data)
    
    if (!is.null(data)) {
        if(replace) {
            setData(data)
        } else {
            data <- processDatasetNames(c(getData(), data))
            setData(data)
        }
        
        if (all(ASquant != "")) {
            # Set species and genome assembly for loaded dataset
            setCategory(category)
            setSpecies(input$customSpecies)
            setAssemblyVersion(input$customAssembly)
        }
    }
    endProcess("loadMultipleFiles", time)
}

#' @rdname appServer
#' @importFrom shiny updateTextInput
localDataServer <- function(input, output, session) {
    ns <- session$ns
    
    ## Load data based on individual files
    prepareFileBrowser(session, input, "sampleInfo")
    prepareFileBrowser(session, input, "subjectInfo")
    prepareFileBrowser(session, input, "geneExpr")
    prepareFileBrowser(session, input, "junctionQuant")
    prepareFileBrowser(session, input, "ASquant")
    
    ## Load data based on folder
    prepareFileBrowser(session, input, "localFolder", directory=TRUE)
    
    # # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))
    
    # Update category name input based on given folder
    observe({
        folder <- input$localFolder
        if (!is.null(folder))
            updateTextInput(session, "localCategory", value=basename(folder))
    })
    
    # If data is loaded, let user replace or append to loaded data
    observeEvent(input$loadMultipleFiles, {
        if (!is.null(getData())) {
            loadedDataModal(session,
                            "localDataModal",
                            "multipleFilesReplace",
                            "multipleFilesLocalAppend")
        } else {
            setMultipleFilesData(input, output, session)
        }
    })
    
    # Load data when the user presses to replace data
    observeEvent(input$multipleFilesReplace,
                 setMultipleFilesData(input, output, session, replace=TRUE))
    
    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$multipleFilesLocalAppend,
                 setMultipleFilesData(input, output, session, replace=FALSE))
    
    # If data is loaded, let user replace or append to loaded data
    observeEvent(input$acceptFile, {
        folder <- input$localFolder
        if (!dir.exists(folder)) {
            # Folder not found
            errorModal(session, "Folder not found", "Check if path is correct.",
                       modalId="localDataModal")
            enable("acceptFile")
            return(NULL)
        } else if (!is.null(getData())) {
            loadedDataModal(session,
                            "localDataModal",
                            "localReplace",
                            "localAppend")
        } else {
            setLocalData(input, output, session)
        }
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