name <- "Data"
id <- function(value) objectId(name, value)

# # Loads valid scripts from the indicated folder
dataName <- "data"
dataEnvs <- sourceScripts(folder = paste0(tabsFolder, "data/"),
                          check = c(dataName, "ui"),
                          parentEnv = environment())
dataEnvs.server <- lapply(dataEnvs, "[[", "server")

# Get name of the loaded scripts
names <- sapply(dataEnvs, "[[", dataName)

ui <- function(tab) {
    tab(name,
        sidebarLayout(
            sidebarPanel(
                tabsetPanel(
                    type = "pill",
                    tabPanel(dataEnvs[[1]]$name, br(), dataEnvs[[1]]$ui()),
                    tabPanel(dataEnvs[[2]]$name, br(), dataEnvs[[2]]$ui()),
                    tabPanel("Exon/intron inclusion", br(), p("Test"))
                )),
            mainPanel( uiOutput(id("tablesOrAbout")) )))
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param title Character: tab title
#' @param tableId Character: id of the datatable
#' @param description Character: description of the table (optional)
#' @param ... Extra arguments to pass to the function dataTableOutput
#'
#' @return The HTML code for a tabPanel template
tabTable <- function(title, id, columns, description = NULL) {
    tablename <- id(paste("table", id, sep = "-"))
    if(!is.null(description))
        d <- p(tags$strong("Table description:"), description, hr())
    else
        d <- NULL
    tabPanel(title, br(), d, dataTableOutput(tablename))
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
                             choices = names(getData())),
                 uiOutput(id("datatabs")))
    })

    # Set the category of the data when possible
    observeEvent(input[[id("category")]], setCategory(input[[id("category")]]))

    # Render tables when data changes
    observe({
        data <- getData()
        categoryData <- getCategoryData()
        for (group in seq_along(data)) {
            lapply(seq_along(categoryData), renderDataTab,
                   data = data[[group]], group)
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
                         columns = names(categoryData[[i]]),
                         id = paste(category, i, sep = "-"),
                         description = attr(categoryData[[i]], "description")))
        do.call(tabsetPanel, c(id = id("dataTypeTab"), dataTablesUI))
    }) # end of renderUI

    # Render a specific data tab (including data table and related interface)
    renderDataTab <- function(index, data, group) {
        tablename <- id(paste("table", getCategories()[group],
                              index, sep = "-"))

        table <- data[[index]]
        if (isTRUE(attr(table, "rowNames")))
            table <- cbind(Row = rownames(table), table)

        # Only show default columns if they are defined
        subsetToShow <- table
        colsToShow <- attr(table, "show")
        if (!is.null(colsToShow)) {
            match <- colsToShow %in% colnames(table)
            subsetToShow <- subset(table, select = colsToShow[match])
        }

        output[[tablename]] <- renderDataTable(
            subsetToShow, options = list(pageLength = 10, scrollX=TRUE))
    } # end of renderDataTab
}