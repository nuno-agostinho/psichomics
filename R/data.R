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

groupsUI <- function(table, tablename) {
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", id("subsetBy"), sign, what)
    
    list(
        hr(),
        fluidRow(
            column(2, selectizeInput(id("subsetBy"), "Subset by",
                                     c("Column", "Rows", "Expression"))),
            conditionalPanel(
                checkId("==", "Column"),
                column(4, selectizeInput(id("groupColumn"), "Select column",
                                         choices = names(table))),
                column(1, icon2("info-circle", id = id("info-circle")),
                       bsTooltip(id("info-circle"),
                                 "Groups will be created automatically depending on the given column."))),
            conditionalPanel(
                checkId("==", "Rows"),
                column(3, textInput(id("groupRows"), "Select rows"),
                       bsTooltip(id("groupRows"),
                                 "Select rows like in R. To create a group with rows 1 to 6, 8, 10 to 19, but not 17, insert 1:6, 8, 10:19, -17"))),
            conditionalPanel(
                checkId("==", "Expression"),
                column(3, textInput(id("groupExpression"), "Subset expression"),
                       bsTooltip(id("groupExpression"),
                                 'To select rows where column X4 is higher than 8 and "alive" in X7, type X4 > 8 & X7 == "alive"'))),
            column(2, conditionalPanel(checkId("!=", "Column"),
                                       textInput(id("groupName"), "Group name"))),
            column(2, actionButton(id("createGroup"), "Create group"))),
        # Align the "create group" button with other inputs
        tags$style(type='text/css', paste0(
            "#", id("createGroup"), " { width:100\\%; margin-top: 25px;}")),
        tags$style(type='text/css',
                   "#", id("icon-circle"), " { width:100\\%; margin-top: 35px;}"),
        tableOutput(id("groupsTable"))
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
                    open = "Add data",
                    shinyBS::bsCollapsePanel(
                        style = "info",
                        title = "Add local files",
                        #list(img(src = "add-button-16px.png"))
                        addLocalFile())
                ),
                shinyBS::bsCollapse(
                    id = id("addTCGA"),
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
tabTable <- function(title, id, description = NULL, ...) {
    tablename <- id(paste("table", id, sep = "-"))
    if(!is.null(description))
        d <- p(tags$strong("Table description:"), description,
               uiOutput(paste0(tablename, "-groups")),
               dataTableOutput(paste0(tablename, "-groupsList")), hr())
    else
        d <- NULL
    tabPanel(title, br(), d,
             dataTableOutput(tablename, ...))
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
            lapply(seq_along(categoryData), renderData,
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
            lapply(seq_along(names(categoryData)),
                   function(i) {
                       tabTable(names(categoryData)[i],
                                id = paste(category, i, sep = "-"),
                                description = attr(categoryData[[i]],
                                                   "description"))
                   })
        )
    }) # end of renderUI
    
    # Render a specific data table
    renderData <- function(index, data, group) {
        tablename <- id(paste("table", getCategories()[group], 
                              index, sep = "-"))
        
        table <- data[[index]]
        if (isTRUE(attr(table, "rowNames")))
            table <- cbind(Row = rownames(table), table)
        
        # output[[paste0(tablename, "-groups")]] <- renderUI(
            # groupsUI(table, tablename))
        
        # Subset to show default columns if any
        if (!is.null(attr(table, "show")))
            showTable <- subset(table, select = attr(table, "show"))
        else
            showTable <- table
        
        output[[tablename]] <- renderDataTable(
            showTable, options = list(pageLength = 10, scrollX=TRUE))
    } # end of renderData
    
    ## TODO(NunoA): currently, groups are implemented for each data category,
    ## but they should be implemented for each data type in a data category
    # Create groups from a data table
    observeEvent(input$createGroup, {
        if (is.null(getGroups())) {
            # Define initial matrix
            groups <- matrix(ncol = 3)
            colnames(groups) <- c("Names", "Subset", "Input")
            setGroups(groups[-1, ])
        }
        
        subsetInput <- input[[id("subsetBy")]]
        groupColumnInput <- input[[id("groupColumn")]]
        if (subsetInput == "Column") {
            categoryData <- getCategoryData()
            names <- unique(categoryData[[groupColumnInput]])
            groups <- cbind(names, subsetInput, groupColumnInput)
        } else {
            elem <- switch(subsetInput,
                           "Rows" = input$groupRows,
                           "Expression" = input$groupExpression)
            row <- c(groupColumnInput, subsetInput, elem)
            groups <- rbind(getGroups(), row)
            rownames(groups) <- NULL
        }
        setGroups(groups)
    })
    
    output[[id("groupsTable")]] <- renderTable(getGroups())
}