#' @import shiny shinyBS shinyjs
name <- "Data"
id <- function(id) paste(name, id, sep = "_")

#' Creates a UI set with options to add a file from the local storage
#' 
#' @return A UI set that can be added to a UI definition
addLocalFile <- function() {
    list(fileInput(id("dataFile"), "Choose folder", multiple = T),
         textInput(id("species"), label = "Species", placeholder = "Required"),
         textInput(id("commonName"), label = "Common name"),
         uiOutput(id("testing")),
         actionButton(id("acceptFile"), "Send file")
    ) # end of list
}

groupsUI <- function(tablename, columns) {
    tabId <- function(value) id(paste(tablename, value, sep = "_"))
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", tabId("subsetBy"), sign, what)
    
    list(fluidRow(
        column(2,
               selectizeInput(tabId("subsetBy"), "Subset by",
                              c("Column", "Rows", "Expression", "Grep"))),
        conditionalPanel(
            checkId("==", "Column"),
            column(4, selectizeInput(tabId("groupColumn"), "Select column",
                                     choices = columns)),
            column(1, icon2("info-circle", id = tabId("info-circle")),
                   bsTooltip(tabId("info-circle"),
                             "Groups will be created automatically depending on the given column."))),
        conditionalPanel(
            checkId("==", "Rows"),
            column(3, textInput(tabId("groupRows"), "Select rows"),
                   bsTooltip(tabId("groupRows"),
                             "Select rows like in R. To create a group with rows 1 to 6, 8, 10 to 19, but not 17, insert 1:6, 8, 10:19, -17"))),
        conditionalPanel(
            checkId("==", "Expression"),
            column(3, textInput(tabId("groupExpression"), "Subset expression"),
                   bsTooltip(tabId("groupExpression"),
                             'To select rows where column X4 is higher than 8 and "alive" in X7, type X4 > 8 & X7 == "alive"'))),
        conditionalPanel(
            checkId("==", "Grep"),
            column(3, textInput(tabId("grepExpression"), "GREP expression")),
            column(3, selectizeInput(tabId("grepColumn"),
                                     "Select column to GREP",
                                     choices = columns))),
        conditionalPanel(checkId("!=", "Column"),
                         column(2, textInput(tabId("groupName"), "Group name"))),
        column(2, actionButton(tabId("createGroup"), "Create group")),
        # Align the "create group" button with other inputs
        tags$style(type='text/css', paste0(
            "#", tabId("createGroup"), " { width:100\\%; margin-top: 25px;}")),
        tags$style(type='text/css', paste0(
            "#", tabId("info-circle"),
            " { width:100\\%; margin-top: 35px;}"))),
        dataTableOutput(tabId("groupsTable"))
    )
}

#' Creates a UI set with options to add data from TCGA/Firehose
#' 
#' @return A UI set that can be added to a UI definition
addTCGAdata <- function() {
    cohorts <- getFirehoseCohorts()
    names(cohorts) <- sprintf("%s (%s)", names(cohorts), cohorts)
    
    if (isFirehoseUp()) {
        list(selectizeInput(id("firehoseCohort"), "Cohort", cohorts,
                            multiple = TRUE, selected = c("ACC", "BLCA"),
                            options = list(placeholder = "Select cohort(s)")),
             selectizeInput(id("firehoseDate"), "Date",
                            as.character(getFirehoseDates()), multiple = TRUE,
                            selected = "2015-11-01", options = list(
                                placeholder = "Select sample date")),
             selectizeInput(id("dataType"), "Data type",
                            c("Clinical", "mRNASeq"), 
                            multiple = TRUE, selected = "Clinical",
                            options = list(
                                placeholder = "Select data types")),
             selectizeInput(id("firehoseExclude"),
                            "Files/archives to exclude", multiple = TRUE,
                            choices = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                        "MANIFEST.txt", "exon_quantification"),
                            selected = c("RSEM_isoforms", ".aux.", ".mage-tab.",
                                         "MANIFEST.txt", "exon_quantification"),
                            # Allow to add new items
                            options = list(
                                create = TRUE, createOnBlur=TRUE,
                                placeholder = "Input files to exclude")),
             textInput(id("dataFolder"), "Folder to store the data",
                       value = "~/Downloads",
                       placeholder = "Insert data folder"),
             bsTooltip(id("dataFolder"), placement = "right", 
                       options = list(container = "body"),
                       "Data not available in this folder will be downloaded."),
             actionButton(id("getFirehoseData"), "Get data"))
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
                    id = id("addData"),
                    open = "Add TCGA/Firehose data",
                    shinyBS::bsCollapsePanel(
                        style = "info",
                        title = list(icon("plus-circle"), "Add local files"),
                        value = "Add local files",
                        addLocalFile()),
                    shinyBS::bsCollapsePanel(
                        style = "info",
                        title = list(icon("plus-circle"),
                                     "Add TCGA/Firehose data"),
                        value = "Add TCGA/Firehose data",
                        addTCGAdata())),
                h3("Data grouping")#,
                # shinyBS::bsCollapse(id = "testttt",
                                    # )
            ),
            mainPanel(
                # TODO(NunoA): Show alerts from renderUI
                bsAlert(anchorId = id("alert2")),
                bsModal2(id("dataReplace"), "Data already loaded", NULL,
                         size = "small",
                         "Would you like to replace the loaded data?",
                         footer = list(
                             actionButton(id("replace"),
                                          "data-dismiss"="modal", 
                                          label = "Replace"))),
                uiOutput(id("tablesOrAbout")),
                uiOutput(id("iframe")))))
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
tabTable <- function(title, id, columns, description = NULL) {
    tablename <- id(paste("table", id, sep = "-"))
    if(!is.null(description))
        d <- p(tags$strong("Table description:"), description, hr(),
               groupsUI(tablename, columns),
               dataTableOutput(paste0(tablename, "-groupsList")), hr())
    else
        d <- NULL
    tabPanel(title, br(), d,
             dataTableOutput(tablename))
}

#' Server logic
#' 
#' @return Part of the server logic related to this tab
server <- function(input, output, session){
    observeEvent(input[[id("category")]], setCategory(input[[id("category")]]))
    
    # Show welcome screen when there's no data loaded
    output[[id("tablesOrAbout")]] <- renderUI({
        if(is.null(getData()))
            includeMarkdown("about.md")
        else
            list(selectInput(id("category"), "Select category:",
                             choices = names(getData())),
                 uiOutput(id("datatabs")))
    })
    
    # User files
    observeEvent(input[[id("acceptFile")]], {
        error <- function(msg) { print(msg); return(NULL) }
        if(is.null(input[[id("dataFile")]])) error("No data input selected")
        if(input[[id("species")]] == "") error("Species field can't be empty")
        
        # inFile <- input$dataFile
        # info <- read.table(inFile$datapath, sep = input$sep,
        #                    header = input$header)
    })
    
    # The button is only enabled if it meets the conditions that follow
    observe(toggleState(id("acceptFile"), input[[id("species")]] != ""))
    
    # Load Firehose data
    loadAllData <- reactive({
        shinyjs::disable(id("getFirehoseData"))
        
        # # Direct download by the browser
        # source = paste0("https://support.apple.com/library/APPLE/APPLECARE_ALLGEOS/HT1425/",
        #                 c("sample_iTunes.mov.zip", "sample_iPod.m4v.zip",
        #                   "sample_mpeg4.mp4.zip"))
        # iframe <- function(s) tags$iframe(width=1, height=1, frameborder=0,
        #                                   src=s)
        # output$iframe <- renderUI(lapply(source, iframe))
        
        # Load data from Firehose
        setData(
            loadFirehoseData(
                folder = input[[id("dataFolder")]],
                cohort = input[[id("firehoseCohort")]],
                date = gsub("-", "_", input[[id("firehoseDate")]]),
                data_type = input[[id("dataType")]],
                exclude = input[[id("firehoseExclude")]],
                progress = updateProgress))
        
        closeProgress()
        shinyjs::enable(id("getFirehoseData"))
    }) # end of reactive
    
    # Load data when the user presses to replace data
    observeEvent(input[[id("replace")]], loadAllData())
    
    # Render tables when data changes
    observe({
        data <- getData()
        categoryData <- getCategoryData()
        for (group in seq_along(data))
            lapply(seq_along(categoryData), renderDataTab,
                   data = data[[group]], group)
    })
    
    # Check if data is already loaded and ask the user if it should be replaced
    observeEvent(input[[id("getFirehoseData")]], {
        if (!is.null(getData()))
            toggleModal(session, id("dataReplace"), "open")
        else
            loadAllData()
    })
    
    # Render tabs with data tables
    output[[id("datatabs")]] <- renderUI({
        categoryData <- getCategoryData()
        category <- getCategory()
        do.call(
            tabsetPanel,
            c(id = id("dataTypeTab"),
              lapply(seq_along(names(categoryData)),
                     function(i)
                         tabTable(names(categoryData)[i],
                                  columns = names(categoryData[[i]]),
                                  id = paste(category, i, sep = "-"),
                                  description = attr(categoryData[[i]],
                                                     "description"))
                     )
              )
        )
    }) # end of renderUI
    
    # Render a specific data tab (including data table and related interface)
    renderDataTab <- function(index, data, group) {
        tablename <- id(paste("table", getCategories()[group], 
                              index, sep = "-"))
        
        # group UI
        tabId <- function(value) id(paste(tablename, value, sep = "_"))
        observeEvent(input[[tabId("createGroup")]], { setGroups(input, tabId) })
        output[[tabId("groupsTable")]] <- renderDataTable(
            getGroupsFrom(input[[id("dataTypeTab")]]),
            options = list(pageLength = 10, scrollX = TRUE))
        
        table <- data[[index]]
        if (isTRUE(attr(table, "rowNames")))
            table <- cbind(Row = rownames(table), table)
        
        # Subset to show default columns if any
        if (!is.null(attr(table, "show")))
            showTable <- subset(table, select = attr(table, "show"))
        else
            showTable <- table
        
        output[[tablename]] <- renderDataTable(
            showTable, options = list(pageLength = 10, scrollX=TRUE))
    } # end of renderDataTab
}

setGroups <- function (input, tabId) {
    # Get groups for the visible and active data table
    active <- input[[id("dataTypeTab")]]
    groups <- getGroupsFrom(active)
    
    if (is.null(groups)) {
        # Define initial matrix
        groups <- matrix(ncol = 3)
        colnames(groups) <- c("Names", "Subset", "Input")
        groups <- groups[-1, ]
    }
    
    subsetInput <- input[[tabId("subsetBy")]]
    if (subsetInput == "Column") {
        groupColumnInput <- input[[tabId("groupColumn")]]
        categoryData <- getCategoryData()
        names <- unique(categoryData[[active]][[groupColumnInput]])
        each <- cbind(names, subsetInput, groupColumnInput)
        groups <- rbind(groups, each)
    } else {
        elem <- switch(subsetInput,
                       "Rows" = input[[tabId("groupRows")]],
                       "Expression" = input[[tabId("groupExpression")]])
        row <- c(input[[tabId("groupName")]], subsetInput, elem)
        groups <- rbind(groups, row)
        rownames(groups) <- NULL
    }
    setGroupsFrom(active, groups)
}