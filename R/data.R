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
        list(selectInput("firehoseCohort", "Cohort", cohorts,
                         multiple = TRUE, selected = c("ACC", "BLCA")),
             selectInput("firehoseDate", "Date",
                         as.character(getFirehoseDates()), multiple = TRUE,
                         selected = "2015-11-01"),
             selectInput("dataType", "Data type",
                         c("Clinical", "mRNASeq"), multiple = TRUE,
                         selected = "Clinical"),
             textInput("firehoseExclude",
                       "Files/archives to exclude (separated by comma)",
                       value = "RSEM_isoforms, .aux., .mage-tab., MANIFEST.txt, exon_quantification"),
             textInput("dataFolder", "Folder where the data is located",
                       value = "~/Downloads"),
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
                bsModal2("modalExample", "Data already loaded",
                         "tabBut", size = "small",
                         "Would you like to replace the loaded data?",
                         footer = list(
                             actionButton("replace", "data-dismiss"="modal", label = "Replace"))),
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
    observe({
        # The button is only enabled if it meets the conditions that follow
        toggleState("acceptFile",
                    input$species != "")
    })
    
    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(shared.data$data)) {
            includeMarkdown("about.md")
        } else {
            list(selectInput("category", "Select category:",
                             choices = names(shared.data$data)),
                 uiOutput("datatabs"))
        }
    }) # end of renderUI
    
    observeEvent(input$acceptFile, {
        #         output$testing <- renderUI({
        #             list(
        #                 badge(inputId="badge1", Sys.time()),
        #                 buttonGroups(actionButton("test111", "Left"),
        #                              actionButton("test222", "Middle"),
        #                              actionButton("test333", "Right")),
        #                 progressbar(sample(1:100, 1)),
        #                 dropdown(inputId="dropdownMenu1"),
        #                 alertNew(progressbar(sample(1:100, 1)))
        #             )
        #         })
        
        error <- function(msg) { print(msg); return(NULL) }
        if(is.null(input$dataFile)) error("No data input selected")
        if(input$species == "") error("Species field can't be empty")
        
        # inFile <- input$dataFile
        # info <- read.table(inFile$datapath, sep = input$sep,
        #                    header = input$header)
        # createAlert(session, anchorId = "alert2", title = "Yay!",
        #             content = list(progressbar(sample(1:100, 1))),
        #             style = "success", append = FALSE)
    }) # end of observeEvent
    
    loadAllData <- reactive({
        shinyjs::disable("getFirehoseData")
        
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Hang in there...", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        
        updateProgress <- function(message, value = NULL, max = NULL,
                                   divisions = NULL, detail = NULL) {
            if (!is.null(divisions)) {
                shared.data$progress.divisions <- divisions
                return(NULL)
            }
            divisions <- shared.data$progress.divisions
            if (is.null(value)) {
                value <- progress$getValue()
                value <- value + (progress$getMax() - value)
            }
            if (is.null(max)) {
                # print(paste(message, value))
                progress$inc(amount = value/divisions, message = message,
                             detail = detail)
            } else {
                # print(paste(message, value, max))
                progress$inc(amount = 1/max/divisions, message = message,
                             detail = detail)
            }
            # print(progress$getValue())
        }
        
        # Parse exclude
        exclude <- strsplit(input$firehoseExclude, ",")[[1]]
        exclude <- trimWhitespace(exclude)
        
        # Load data from Firehose
        shared.data$data <- loadFirehoseData(
            folder = input$dataFolder,
            cohort = input$firehoseCohort,
            date = gsub("-", "_", input$firehoseDate),
            data_type = input$dataType,
            exclude = exclude,
            progress = updateProgress)
        shared.data$progress.divisions <- NULL
        shinyjs::enable("getFirehoseData")
    })
    
    observeEvent(input$replace, { loadAllData() })
    
    # Load Firehose data
    observeEvent(input$getFirehoseData, {
        if (!is.null(shared.data$data)) {
            toggleModal(session, "modalExample", "open")
        } else {
            loadAllData()
        }
    }) # end of observeEvent
    
    # Render tabs with data tables
    output$datatabs <- renderUI({
        data <- shared.data$data[[input$category]]
        do.call(tabsetPanel,
                lapply(seq_along(names(data)), function(i) {
                    tabTable(names(data)[i], id = paste(input$category, i, sep = "."),
                             description = attr(data[[i]], "description"))
                })
        )
    }) # end of renderUI
    
    # Render data tables every time the data changes
    observe({
        # For better performance, try to use conditional panels
        if (!is.null(shared.data$data)) {
            for (k in seq_along(shared.data$data)) {
                data <- shared.data$data[[k]]
                lapply(seq_along(data), function(i) {
                    tablename <- paste("table", names(shared.data$data)[k],
                                       i, sep = ".")
                    print(attributes(data[[i]]))
                    if (isTRUE(attr(data[[i]], "rowNames")))
                        d <- cbind(names = rownames(data[[i]]), data[[i]])
                    else
                        d <- data[[i]]
                    output[[tablename]] <- renderDataTable(d,
                        options = list(pageLength = 10, scrollX=TRUE))
                })
            }
        }
    }) # end of observe
}