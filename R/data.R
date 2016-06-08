## TODO(NunoA): should default columns be a perfect match or just a partial
## match? If only a partial match... that would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

#' User interface
#' @param tab Function to create tab
dataUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "data")
    
    tab(title=div(icon("table"), "Data"),
        sidebarLayout(sidebarPanel(
            do.call(tabsetPanel,
                    c(type="pill", uiList)
            )
        ),
        mainPanel( uiOutput(ns("tablesOrAbout")) ))
    )
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param title Character: tab title
#' @param tableId Character: id of the datatable
#' @param description Character: description of the table (optional)
#' @param columns Character: column names of the datatable
#' @param colsToShow Boolean: columns to show
#'
#' @importFrom shinyBS bsTooltip
#'
#' @return The HTML code for a tabPanel template
tabDataset <- function(ns, title, tableId, columns, colsToShow,
                       description=NULL) {
    tablename <- ns(paste("table", tableId, sep="-"))
    
    if(!is.null(description))
        description <- p(tags$strong("Table description:"), description)
    
    downloadId <- paste(tablename, "download", sep="-")
    download <- downloadButton(downloadId, "Download this dataset")
    
    visibleColumns <- selectizeInput(paste(tablename, "columns", sep="-"),
                                     label="Columns to show", choices=columns,
                                     selected=colsToShow, multiple=TRUE,
                                     width="100%")
    
    tooltip <- bsTooltip(paste(tablename, "columns", sep="-"),
                         paste("Dataset columns to show; empty this input to",
                               "show all columns (it can be slow for large",
                               "datasets)"),
                         placement = "top", options = list(container = "body"))
    
    tabPanel(title, br(), description, download, br(), visibleColumns, tooltip,
             hr(), dataTableOutput(tablename))
}

#' Render a specific data tab (including data table and related interface)
#' @param index Integer: index of the data to load
#' @param data Data frame: data with everything to load
#' @param name Character: name of the dataset
#' @param input Shiny session input
#' @param output Shiny session output
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
    
    # Show row names if there are any
    if (isTRUE(attr(table, "rowNames")))
        subsetToShow <- cbind(Row=rownames(subsetToShow), subsetToShow)
    
    output[[tablename]] <- renderDataTable(
        subsetToShow, options=list(pageLength=10, scrollX=TRUE))
    
    output[[paste(tablename, "download", sep="-")]] <- downloadHandler(
        filename = paste(name, attr(table, "tablename")),
        content = function(file) write.table(table, file, quote = FALSE,
                                             row.names = TRUE, sep = "\t"))
}

#' Server logic
#'
#' @return Part of the server logic related to this tab
dataServer <- function(input, output, session) {
    ns <- session$ns
    
    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(getData())) {
            includeMarkdown(insideFile("shiny", "about.md"))
        } else
            list(selectInput(ns("category"), "Select category:",
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
                tabDataset(ns,
                           names(categoryData)[i],
                           paste(category, i, sep="-"),
                           names(categoryData[[i]]),
                           attr(categoryData[[i]], "show"),
                           description=attr(categoryData[[i]], "description"))
            })
        do.call(tabsetPanel, c(id=ns("datasetTab"), dataTablesUI))
    })
    
    # Run server logic from the scripts
    getServerFunctions("data", reactive(input$datasetTab))
}

attr(dataUI, "loader") <- "app"
attr(dataServer, "loader") <- "app"