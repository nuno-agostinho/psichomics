## TODO(NunoA): should default columns be a perfect match or just a partial
## match? If only a partial match... that would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

name <- "Data"
id <- function(value) objectId(name, value)

## Loads valid scripts from the indicated folder
dataName <- "data"
dataEnvs <- sourceScripts(folder=paste0(tabsFolder, "data/"),
                          check=c(dataName, "ui"),
                          parentEnv=environment())
dataEnvs.server <- lapply(dataEnvs, "[[", "server")

# Get name of the loaded scripts
names <- sapply(dataEnvs, "[[", dataName)

ui <- function(tab) {
    tab(title=div(icon("table"), name), sidebarLayout(sidebarPanel(
        do.call(
            tabsetPanel,
            c(type="pill",
              lapply(dataEnvs, function(pill)
                  tabPanel(pill$name, br(), pill$ui()))
            ))),
        mainPanel( uiOutput(id("tablesOrAbout")) )))
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param title Character: tab title
#' @param id Character: id of the datatable
#' @param description Character: description of the table (optional)
#'
#' @return The HTML code for a tabPanel template
tabTable <- function(title, id, description=NULL) {
    tablename <- id(paste("table", id, sep="-"))
    if(!is.null(description))
        descr <- p(tags$strong("Table description:"), description)
    else
        descr <- NULL
    
    download <- downloadButton(paste(tablename, "download", sep="-"),
                               "Download")
    tabPanel(title, br(), descr, download, hr(), dataTableOutput(tablename))
}

#' Server logic
#'
#' @return Part of the server logic related to this tab
server <- function(input, output, session) {
    # Runs server logic from the scripts
    lapply(dataEnvs.server, do.call, list(input, output, session))
    
    # Show welcome screen when there's no data loaded
    output[[id("tablesOrAbout")]] <- renderUI({
        if(is.null(getData()))
            includeMarkdown("about.md")
        else
            list(selectInput(id("category"), "Select category:",
                             choices=names(getData())),
                 uiOutput(id("datatabs")))
    })
    
    # Set the category of the data when possible
    observeEvent(input[[id("category")]], setCategory(input[[id("category")]]))
    
    # Render tables when data changes
    observe({
        data <- getData()
        categoryData <- getCategoryData()
        for (group in seq_along(data)) {
            lapply(seq_along(categoryData), renderDataTab, data=data[[group]],
                   group)
        }
    })
    
    # Render tabs with data tables
    output[[id("datatabs")]] <- renderUI({
        categoryData <- getCategoryData()
        category <- getCategory()
        
        dataTablesUI <- lapply(
            seq_along(names(categoryData)),
            function(i)
                tabTable(names(categoryData)[i],
                         # columns=names(categoryData[[i]]),
                         id=paste(category, i, sep="-"),
                         description=attr(categoryData[[i]], "description")))
        do.call(tabsetPanel, c(id=id("dataTypeTab"), dataTablesUI))
    }) # end of renderUI
    
    # Render a specific data tab (including data table and related interface)
    renderDataTab <- function(index, data, group) {
        tablename <- id(paste("table", getCategories()[group], index, sep="-"))
        
        table <- data[[index]]
        # Only show default columns if they are defined (don't cause problems)
        subsetToShow <- table
        colsToShow <- attr(table, "show")
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
            filename = paste(getCategories()[group], attr(table, "tablename")),
            content = function(file) write.table(table, file, quote = FALSE,
                                                 row.names = TRUE, sep = "\t"))
    } # end of renderDataTab
}