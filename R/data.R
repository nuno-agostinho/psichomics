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
                bsModal("modalExample", "There's already data, dumb!",
                        "tabBut", size = "large", p("Yeah! Dumb!")),
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
tabTable <- function(title, tableId, description = NULL, ...) {
    if(!is.null(description))
        d <- p(tags$strong("Table description:"), description)
    else
        d <- NULL
    tabPanel(title, br(), d, dataTableOutput(tableId, ...))
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
    
    output$tablesOrAbout <- renderUI({
        # If no data file is loaded, show welcome screen
        if(!is.null(shared.data$a)) {
            includeMarkdown("about.md")
        } else {
            tabsetPanel(
                tabTable("Test", "dataTable2",
                         "this is a pretty table with interesting info."),
                tabTable("Another", "dataTable"),
                tabTable("More", "dataTable3")
            )
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
        
        inFile <- input$dataFile
        info <- read.table(inFile$datapath, sep = input$sep,
                           header = input$header)
        shared.data$a <<- new("Classification",
                              species = input$species,
                              common.name = input$common.name,
                              inclusion.levels = info)
        createAlert(session, anchorId = "alert2", title = "Yay!",
                    content = list(progressbar(sample(1:100, 1))),
                    style = "success", append = FALSE)
    }) # end of observeEvent
    
    output$dataTable <- renderDataTable(mtcars)
}