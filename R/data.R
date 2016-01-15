name <- "Data"

#' Creates a collapsible UI set with options to add a file from the local
#' storage or from TCGA.
#' 
#' @param number Integer: A number to identify the set
#' 
#' @return A collapse panel that can be added to a UI definition.
addFileCollapse <- function(number) {
    shinyBS::bsCollapsePanel(
        style = "info",
        title = "Add file", #list(img(src = "add-button-16px.png"),
        #     paste("Add file", number)),
        h4("Add local files"),
        fileInput(paste0("dataFile", number), "Choose folder", multiple = T),
        selectInput("sep", "Choose separator", 
                    choices = list("Tab" = "\t", "Comma"=",", "Space"=" "),
                    selected = "\t"),
        uiOutput("testing"),
        actionButton(paste0("acceptFile", number), "Send file")
    ) # end of bsCollapsePanel
}

ui <- function() {
    tabPanel(
        name,
        sidebarLayout(
            sidebarPanel(
                h3("File input"),
                shinyBS::bsCollapse(
                    id = "addFiles",
                    addFileCollapse(1)
                    # addFileCollapse(2)
                )
            ),
            mainPanel(
                # TODO(NunoA): Show alerts from renderUI
                bsAlert(anchorId = "alert2"),
                uiOutput("tableOrAbout")
            )
        )
    )
}

#' Server logic
#' 
#' @return Part of the server logic related to this tab
server <- function(input, output, session){
    thisData <- reactiveValues(data=NULL)
    
    observeEvent(input$acceptFile1, {
        output$testing <- renderUI({
            list(
                badge(inputId="badge1", Sys.time()),
                buttonGroups(actionButton("test111", "Left"),
                             actionButton("test222", "Middle"),
                             actionButton("test333", "Right")),
                progressbar(sample(1:100, 1)),
                dropdown(inputId="dropdownMenu1")
            )
        })
        validate(need(input$dataFile1, label="No file selected"))
        inFile <- input$dataFile1
        thisData$data <- read.table(inFile$datapath, sep = input$sep)
        createAlert(session, anchorId = "alert2", title = "Yay!",
                    content = "File successfully imported",
                    style = "success", append = FALSE)
    }) # end of observeEvent
    
    output$tableOrAbout <- renderUI({
        # If no data file is loaded, show welcome screen
        if(is.null(thisData$data))
            includeMarkdown("about.md")
        else
            dataTableOutput("dataTable")
    })
    
    output$dataTable <- renderDataTable({
        validate(if(is.null(thisData$data)) return(NULL))
        thisData$data
    })
}
