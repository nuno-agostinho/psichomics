#' @import shiny shinyBS shinyjs
name <- "Data"

#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function() {
    list(fileInput("dataFile", "Choose folder", multiple = T),
         textInput("species", label = "Species", placeholder = "Required"),
         textInput("common.name", label = "Common name"),
         uiOutput("testing"),
         actionButton("acceptFile", "Send file")
    ) # end of list
}

#' Creates a UI set with options to add data from TCGA/Firehose
#' 
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function() {
    cohorts <- getFirehoseCohorts()
    names(cohorts) <- sprintf("%s (%s)", names(cohorts), cohorts)
    
    if (isFirehoseUp()) {
        list(selectizeInput("firehoseCohort", "Cohort", cohorts,
                            multiple = TRUE, selected = c("ACC", "BLCA"),
                            options = list(placeholder = "Select cohort(s)")),
             selectizeInput("firehoseDate", "Date",
                            as.character(getFirehoseDates()), multiple = TRUE,
                            selected = "2015-11-01", options = list(
                                placeholder = "Select sample date")),
             selectizeInput("dataType", "Data type",
                            c("Clinical", "mRNASeq"), multiple = TRUE,
                            selected = "Clinical", options = list(
                                placeholder = "Select data types")),
             selectizeInput("firehoseExclude",
                            "Files/archives to exclude", multiple = TRUE,
                            choices = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                        "MANIFEST.txt", "exon_quantification"),
                            selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                         "MANIFEST.txt", "exon_quantification"),
                            # Allow to add new items
                            options = list(
                                create = TRUE, createOnBlur=TRUE,
                                placeholder = "Input files to exclude")),
             textInput("dataFolder", "Folder to store the data",
                       value = "~/Downloads",
                       placeholder = "Insert data folder"),
             bsTooltip("dataFolder", placement = "right", 
                       options = list(container = "body"),
                       "Data not available in this folder will be downloaded."),
             actionButton("getFirehoseData", "Get data"))
    } else {
        list(p("Not able to reach Firehose."))
    }
}

ui <- function(tab) {
    tab(name,
        sidebarLayout(
            sidebarPanel(
                h3("Data input"),
                shinyBS::bsCollapse(
                    id = "addData",
                    open = "Add data",
                    shinyBS::bsCollapsePanel(
                        style = "info",
                        title = "Add local files",
                        #list(img(src = "add-button-16px.png"))
                        addLocalFile())
                ),
                shinyBS::bsCollapse(
                    id = "addTCGA",
                    open = "Add TCGA/Firehose data",
                    shinyBS::bsCollapsePanel(
                        style = "info",
                        title = "Add TCGA/Firehose data",
                        #list(img(src = "add-button-16px.png"))
                        addTCGAdata())
                )
            ),
            mainPanel(
                # TODO(NunoA): Show alerts from renderUI
                bsAlert(anchorId = "alert2"),
                bsModal2("dataReplace", "Data already loaded", NULL,
                         size = "small",
                         "Would you like to replace the loaded data?",
                         footer = list(
                             actionButton("replace",
                                          "data-dismiss"="modal", 
                                          label = "Replace"))),
                uiOutput("tablesOrAbout")
            )
        )
    )
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param title Character: tab title
#' @param tableId Character: id of the datatable
#' @param description Character: description of the table (optional)
#' @param ... Extra arguments to pass to the function dataTableOutput
#'
#' @return The HTML code for a tabPanel template
#' @export
tabTable <- function(title, id, description = NULL, ...) {
    if(!is.null(description))
        d <- p(tags$strong("Table description:"), description, hr())
    else
        d <- NULL
    tablename <- paste("table", id, sep = ".")
    tabPanel(title, br(), d, dataTableOutput(tablename, ...))
}

#' Server logic
#' 
#' @return Part of the server logic related to this tab
server <- function(input, output, session){
    observeEvent(input$category, setCategory(input$category))
    
    observe({
        # The button is only enabled if it meets the conditions that follow
        toggleState("acceptFile", input$species != "")
    })
    
    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(getData())) {
            includeMarkdown("about.md")
        } else {
            list(selectInput("category", "Select category:",
                             choices = names(getData())),
                 uiOutput("datatabs"))
        }
    }) # end of renderUI
    
    observeEvent(input$acceptFile, {
        error <- function(msg) { print(msg); return(NULL) }
        if(is.null(input$dataFile)) error("No data input selected")
        if(input$species == "") error("Species field can't be empty")
        
        # inFile <- input$dataFile
        # info <- read.table(inFile$datapath, sep = input$sep,
        #                    header = input$header)
    }) # end of observeEvent
    
    loadAllData <- reactive({
        shinyjs::disable("getFirehoseData")
        
        # Load data from Firehose
        setData(
            loadFirehoseData(
                folder = input$dataFolder,
                cohort = input$firehoseCohort,
                date = gsub("-", "_", input$firehoseDate),
                data_type = input$dataType,
                exclude = input$firehoseExclude,
                progress = updateProgress))
        
        closeProgress()
        shinyjs::enable("getFirehoseData")
    })
    
    # Load data when the user presses to replace data
    observeEvent(input$replace, loadAllData())
    
    # Load Firehose data
    observeEvent(input$getFirehoseData, {
        if (!is.null(getData())) {
            toggleModal(session, "dataReplace", "open")
        } else {
            loadAllData()
        }
    }) # end of observeEvent
    
    # Render tabs with data tables
    output$datatabs <- renderUI({
        data <- getCategoryData()
        do.call(
            tabsetPanel,
            lapply(seq_along(names(data)),
                   function(i) {
                       tabTable(names(data)[i],
                                id = paste(input$category, i, sep = "."),
                                description = attr(data[[i]], "description"))
                   })
        )
    }) # end of renderUI
    
    # Render a specific data table from sharedData
    renderData <- function(index, data, group) {
        tablename <- paste("table", names(getData())[group],
                           index, sep = ".")
        
        if (isTRUE(attr(data[[index]], "rowNames")))
            table <- cbind(names = rownames(data[[index]]), data[[index]])
        else
            table <- data[[index]]
        
        # Subset to show default columns if any
        if (!is.null(attr(table, "show")))
            table <- subset(table, select = attr(table, "show"))
        
        output[[tablename]] <- renderDataTable(
            table, options = list(pageLength = 10, scrollX=TRUE))
    }
    
    # Render data tables every time the data changes
    observe({
        for (group in seq_along(getData())) {
            data <- getData()[[group]]
            lapply(seq_along(data), renderData, data, group)
        }
    })
}