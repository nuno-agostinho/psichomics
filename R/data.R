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
        checkboxInput("header", "Has header", value = FALSE),
        textInput("species", label = "Species", placeholder = "Required"),
        textInput("common.name", label = "Common name"),
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
                    open = "Add file",
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
    observeEvent(input$acceptFile1, {
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
        if(is.null(input$dataFile1)) return(error("No data input selected"))
        if(input$species == "") return(error("Species field can't be empty"))
        
        inFile <- input$dataFile1
        info <- read.table(inFile$datapath, sep = input$sep,
                           header = input$header)
        shared.data$a <<- new("Organism",
                              species = input$species,
                              common.name = input$common.name,
                              inclusion.levels = info)
        createAlert(session, anchorId = "alert2", title = "Yay!",
                    content = list(progressbar(sample(1:100, 1))),
                    style = "success", append = FALSE)
    }) # end of observeEvent
    
    output$tableOrAbout <- renderUI({
        # If no data file is loaded, show welcome screen
        if(is.null(shared.data$a))
            includeMarkdown("about.md")
        else
            dataTableOutput("dataTable")
    })
    
    output$dataTable <- renderDataTable(inclusion.levels(data$a))
}