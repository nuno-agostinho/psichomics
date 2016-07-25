## TODO(NunoA): should default columns be a perfect match or just a partial
## match? If only a partial match... that would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

#' Create a modal warning the user of already loaded data
#' @param modalId Character: identifier of the modal
#' @param replaceButtonId Character: identifier of the button to replace data
#' @param keepButtonId Character: identifier of the button to append data
#' @param session Shiny session
loadedDataModal <- function(session, modalId, replaceButtonId, keepButtonId) {
    ns <- session$ns
    warningModal(session, "Data already loaded",
                 "Would you like to", tags$b("replace"), "the loaded data or",
                 tags$b("keep"), "both the previous and new data?",
                 footer = tagList(
                     actionButton(ns(keepButtonId), "data-dismiss"="modal",
                                  label="Keep both"),
                     actionButton(ns(replaceButtonId), class = "btn-warning",
                                  "data-dismiss"="modal", label="Replace")),
                 modalId=modalId)
}

#' User interface of the data module
#' @param id Character: identifier
#' @param tab Function to create tab
#' @return HTML elements
dataUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "data", bsCollapsePanel,
                             priority="localDataUI")
    
    tab(title=div(icon("table"), "Data"),
        sidebarLayout(
            sidebarPanel(
                do.call(bsCollapse, uiList)
            ),
            mainPanel( uiOutput(ns("tablesOrAbout")) ))
    )
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param ns Namespace function
#' @param title Character: tab title
#' @param tableId Character: id of the datatable
#' @param description Character: description of the table (optional)
#' @param columns Character: column names of the datatable
#' @param colsToShow Boolean: columns to show
#' @param data Data frame: dataset of interest
#'
#' @importFrom shinyBS bsTooltip
#' @importFrom DT dataTableOutput
#' @importFrom shiny hr br tabPanel selectizeInput column fluidRow p mainPanel
#' downloadButton
#'
#' @return The HTML code for a tabPanel template
tabDataset <- function(ns, title, tableId, columns, colsToShow, data,
                       description=NULL) {
    tablename <- ns(paste("table", tableId, sep="-"))
    
    downloadId <- paste(tablename, "download", sep="-")
    download <- downloadButton(downloadId, "Download dataset", "pull-right")
    
    if(!is.null(description)) {
        description <- p(tags$strong("Table description:"), description)
        download <- fluidRow(
            column(10, description),
            column(2, download)
        )
    }
    
    # Get class of each column
    colType <- sapply(1:ncol(data), function(i) class(data[[i]]))
    colType[colType == "character"] <- "string"
    
    # Show class of each column
    choices <- columns
    names(choices) <- sprintf("%s (%s class)", columns, colType)
    
    visibleColumns <- selectizeInput(
        paste(tablename, "columns", sep="-"), label="Columns to show", 
        choices=choices, selected=colsToShow, multiple=TRUE, width="auto", 
        options=list(plugins=list('remove_button', 'drag_drop'),
                     render=I("{ item: function(item, escape) {
                                 return '<div>' + escape(item.value) + '</div>';
                               } }")))
    
    tooltip <- bsTooltip(paste(tablename, "columns", sep="-"),
                         paste("Dataset columns to show; empty this input to",
                               "show all columns (it can be slow for large",
                               "datasets)"),
                         placement = "top", options = list(container = "body"))
    
    tabPanel(title, br(), download, visibleColumns, tooltip, hr(),
             dataTableOutput(tablename))
}

#' Render a specific data tab (including data table and related interface)
#' @param index Integer: index of the data to load
#' @param data Data frame: data with everything to load
#' @param name Character: name of the dataset
#' @param input Shiny session input
#' @param output Shiny session output
#' 
#' @importFrom DT renderDataTable
#' @importFrom shiny downloadHandler br
createDataTab <- function(index, data, name, input, output) {
    tablename <- paste("table", name, index, sep="-")
    
    table <- data[[index]]
    # Only show default columns if they are defined (don't cause problems)
    subsetToShow <- table
    
    colsToShow <- input[[paste(tablename, "columns", sep="-")]]
    if (!is.null(colsToShow)) {
        match <- colsToShow %in% colnames(table)
        subsetToShow <- subset(table, select=colsToShow[match])
    }
    
    output[[tablename]] <- renderDataTable(
        subsetToShow, style="bootstrap", selection='none', filter="top",
        server=TRUE, options=list(pageLength=10))
    
    output[[paste(tablename, "download", sep="-")]] <- downloadHandler(
        filename = paste(name, attr(table, "tablename")),
        content = function(file) write.table(table, file, quote = FALSE,
                                             row.names = TRUE, sep = "\t"))
}

#' Server logic of the data module
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#'
#' @importFrom shiny selectInput tabsetPanel tags h1 h2
#'
#' @return Part of the server logic related to this tab
dataServer <- function(input, output, session) {
    ns <- session$ns
    
    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(getData())) {
            tagList(
                h1("Welcome"),
                "This app provides an effortless analysis of tumour cells with",
                "a focus on alternative splicing.",
                h2("Instructions"),
                tags$ol(id="list",
                    tags$li("Load tumour transcriptomic data from TCGA",
                            "(clinical data and junction quantification)"),
                    tags$li("Quantify the alternative splicing events"),
                    tags$li("Explore statistically significant events or",
                            "individual events of interest")
                ), br(), br(),
                p(style="text-align:right",
                    tags$a(href="http://imm.medicina.ulisboa.pt/group/compbio/",
                           target="_blank", "Nuno Morais Lab, iMM"), 
                    "(Nuno Agostinho, 2015-2016)", br(),
                    "Special thanks to my lab colleagues for their work-related",
                    br(), "support and supporting chatter."
                )
            )
        } else
            list(selectInput(ns("category"), "Select category:", width="auto",
                             choices=names(getData())),
                 uiOutput(ns("datatabs")))
    })
    
    # Set the category of the data when possible
    observeEvent(input$category, setCategory(input$category))
    
    # Render tables when data changes
    observe({
        data <- getData()
        categoryData <- getCategoryData()
        for (category in seq_along(data)) {
            name <- getCategories()[category]
            # Create data tab for each dataset in a data category
            lapply(seq_along(categoryData), createDataTab,
                   data=data[[category]], name, input, output)
        }
    })
    
    # Render tabs with data tables
    output$datatabs <- renderUI({
        categoryData <- getCategoryData()
        category <- getCategory()
        
        dataTablesUI <- lapply(
            seq_along(categoryData),
            function(i) {
                data <- categoryData[[i]]
                tabDataset(ns, names(categoryData)[i], 
                           paste(category, i, sep="-"), names(data),
                           attr(data, "show"), data,
                           description=attr(data, "description"))
            })
        do.call(tabsetPanel, c(id=ns("datasetTab"), dataTablesUI))
    })
    
    # Change the active dataset
    observe( setActiveDataset(input$datasetTab) )
    
    # Run server logic from the scripts
    getServerFunctions("data", priority="localDataServer")
}

attr(dataUI, "loader") <- "app"
attr(dataServer, "loader") <- "app"