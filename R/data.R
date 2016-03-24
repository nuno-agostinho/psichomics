#' @import shiny shinyBS shinyjs
name <- "Data"
id <- function(...) paste(name, list(...), sep = "_")

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

groupsUI <- function() {
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", id("subsetBy"), sign, what)
    
    list(
        selectizeInput(id("subsetBy"), "Subset by",
                       c("Column", "Rows", "Expression", "Grep")),
        conditionalPanel(
            checkId("==", "Column"),
            selectizeInput(id("groupColumn"), "Select column", choices = NULL),
            icon2("info-circle", id = id("info-circle")),
            bsTooltip(id("info-circle"),
                      paste("Groups will be created automatically depending on",
                            "the given column."))),
        conditionalPanel(
            checkId("==", "Rows"),
            selectizeInput(id("groupRows"), "Select rows", choices = NULL,
                           multiple = TRUE,
                           # Allow to add new items
                           options = list(
                               create = TRUE, createOnBlur=TRUE,
                               ## TODO(NunoA): only allow numbers (use selectize.js REGEX option)
                               # Hide discarded user-created items in the dropdown
                               persist = FALSE),
                           bsTooltip(id("groupRows"),
                                     paste("Select rows like in R. To create a group with rows",
                                           "1 to 6, 8 and 10 to 19, insert 1:6, 8, 10:19")))),
        conditionalPanel(
            checkId("==", "Expression"),
            textInput(id("groupExpression"), "Subset expression"),
            bsTooltip(id("groupExpression"),
                      paste('To select rows where column X4 is higher than 8',
                            'and "alive" in X7, type X4 > 8 & X7 == "alive"'))),
        conditionalPanel(
            checkId("==", "Grep"),
            textInput(id("grepExpression"), "GREP expression"),
            selectizeInput(id("grepColumn"),
                           "Select column to GREP",
                           choices = NULL)),
        conditionalPanel(checkId("!=", "Column"),
                         textInput(id("groupName"), "Group name")),
        actionButton(id("createGroup"), "Create group"),
        uiOutput(id("groupsList"))
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
        tags$script("
                    Shiny.addCustomMessageHandler('getCheckedBoxes',
                        function(variable) {
                            var selected = [];
                            $(\"input[name='checkGroups']:checked\").each(function() {
                                selected.push($(this).attr('number'));
                            });
                            Shiny.onInputChange(variable, selected);
                        });
                    Shiny.addCustomMessageHandler('setZero',
                        function(variable) {
                            Shiny.onInputChange(variable, 0);
                        });"),
        sidebarLayout(
            sidebarPanel(
                tabsetPanel(type = "pill",
                            tabPanel(
                                "Data input", br(),
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
                                        addTCGAdata()))
                            ),
                            tabPanel("Data grouping", br(),
                                     shinyBS::bsCollapse(
                                         id = id("groupingData"),
                                         open = "Grouping",
                                         shinyBS::bsCollapsePanel(title = "Grouping", style = "info",
                                                                  groupsUI()))
                            )
                )),
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
        d <- p(tags$strong("Table description:"), description, hr())
    else
        d <- NULL
    tabPanel(title, br(), d, dataTableOutput(tablename))
}

#' Set new groups according to the user input
#' 
#' The groups are inserted in a matrix
createGroupFromInput <- function (input) {
    active <- input[[id("dataTypeTab")]]
    type <- input[[id("subsetBy")]]
    
    columnInput <- input[[id("groupColumn")]]
    data <- getCategoryData()[[active]]
    
    if (type == "Column") {
        # Get all unique values for given column and its respective rows
        colData <- data[[columnInput]]
        
        # Replace NAs for "NA" so they aren't discarded when splitting the data
        colData[is.na(colData)] <- "NA"
        
        # Split data according to the chosen column
        set <- split(data, colData, drop = FALSE)
        groupNames <- names(set)
        whichRows <- lapply(set, rownames)
        group <- cbind(groupNames, type, columnInput, whichRows)
    } else if (type == "Rows") {
        # Convert the given string into a sequence of numbers
        rows <- input[[id("groupRows")]]
        strRows <- paste(rows, collapse = ", ")
        rows <- unlist(lapply(rows, function(row) eval(parse(text = row))))
        whichRows <- list(rownames(data[rows, ]))
        group <- cbind(input[[id("groupName")]], type, strRows, whichRows)
    } else if (type == "Expression") {
        # Subset data using the given expression
        expr <- input[[id("groupExpression")]]
        set <- subset(data, eval(parse(text = expr)))
        whichRows <- list(rownames(set))
        group <- c(input[[id("groupName")]], type, expr, whichRows)
    } else if (type == "Grep") {
        ## TODO(NunoA): Subset data with the GREP expression for the given column
        group <- rep(NA, 4)
    } 
    # Standarise rows
    colnames(group) <- c("Names", "Subset", "Input", "Rows")
    rownames(groups) <- NULL
    return(group)
}

#' Server logic
#' 
#' @return Part of the server logic related to this tab
server <- function(input, output, session) {
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
    
    # Load user files
    observeEvent(input[[id("acceptFile")]], {
        error <- function(msg) { print(msg); return(NULL) }
        if(is.null(input[[id("dataFile")]]))
            error("No data input selected")
        if(input[[id("species")]] == "")
            error("Species field can't be empty")
        
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
        
        table <- data[[index]]
        if (isTRUE(attr(table, "rowNames")))
            table <- cbind(Row = rownames(table), table)
        
        # Only show default columns if they are defined
        if (!is.null(attr(table, "show")))
            showTable <- subset(table, select = attr(table, "show"))
        else
            showTable <- table
        
        output[[tablename]] <- renderDataTable(
            showTable, options = list(pageLength = 10, scrollX=TRUE))
    } # end of renderDataTab
    
    # Update columns available for creating groups when there's loaded data
    observeEvent(input[[id("dataTypeTab")]], {
        active <- input[[id("dataTypeTab")]]
        for (i in id("groupColumn", "grepColumn")) {
            updateSelectizeInput(session, i,
                                 choices = names(getCategoryData()[[active]]))
        }
    })
    
    # Create a new group when clicking on the createGroup button
    observeEvent(input[[id("createGroup")]], {
        # Create new group(s) from user input
        new <- createGroupFromInput(input)
        
        # Get groups for the data table that is visible and active
        active <- input[[id("dataTypeTab")]]
        groups <- getGroupsFrom(active)
        
        # Append the new group(s) to the groups already created
        groups <- rbind(groups, new)
        setGroupsFrom(active, groups)
    })
    
    # Render groups list and add features to merge/remove groups
    output[[id("groupsTable")]] <- renderDataTable({
        ## TODO(NunoA): Allow to remove and merge selected rows from the groups
        ## This could be done using checkboxes; how to retrieve which checkboxes
        ## were checked? Possible with data table or javascript?
        active <- input[[id("dataTypeTab")]]
        groups <- getGroupsFrom(active)
        
        # Don't show anything if there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            # Show number of rows for each group
            rows <- lapply(groups[ , 4], length)
            groups[ , 4] <- unlist(rows)
            # Add checkboxes
            pick <- paste("<input number=", 1:nrow(groups),
                          " name='checkGroups' type='checkbox'></input>")
            return(cbind("Pick" = pick, groups))
        }
    }, options = list(pageLength = 10, lengthChange = FALSE, scrollX = TRUE, 
                      filter = FALSE, info = FALSE, paginationType = "simple"),
    escape = FALSE)
    
    # Remove selected groups when pressing the button
    observeEvent(input[[id("removeGroups")]], {
        session$sendCustomMessage(type = "getCheckedBoxes", "removeGroups")
        sharedData$removeTime <- TRUE
    })
    
    # Remove selected groups if there are groups to be removed
    observe({
        if (!is.null(input$removeGroups) && all(input$removeGroups > 0) &&
            isTRUE(sharedData$removeTime)) {
            # Set groups to remove to 0 and flag to FALSE
            session$sendCustomMessage(type = "setZero", "removeGroups")
            sharedData$removeTime <- FALSE
            
            # Get groups for the data table that is visible and active
            active <- input[[id("dataTypeTab")]]
            groups <- getGroupsFrom(active)
            
            # Remove selected groups
            selected <- as.numeric(input$removeGroups)
            setGroupsFrom(active, groups[-selected, , drop=FALSE])
        }
    })
    
    # Merge selected groups when pressing the button
    observeEvent(input[[id("mergeGroups")]], {
        session$sendCustomMessage(type = "getCheckedBoxes", "mergeGroups")
        sharedData$mergeTime <- TRUE
    })
    
    # Merge selected groups if there are groups to be merged
    observe({
        if (!is.null(input$mergeGroups) && all(input$mergeGroups > 0) &&
            isTRUE(sharedData$mergeTime)) {
            # Set groups to remove to 0 and flag to FALSE
            session$sendCustomMessage(type = "setZero", "mergeGroups")
            sharedData$mergeTime <- FALSE
            
            # Get groups for the data table that is visible and active
            active <- input[[id("dataTypeTab")]]
            groups <- getGroupsFrom(active)
            
            # Create merged group
            selected <- as.numeric(input$mergeGroups)
            mergedFields <- lapply(1:3, function(i)
                paste(groups[selected, i], collapse = " + "))
            rowNumbers <- sort(Reduce(union, groups[selected, 4]))
            new <- matrix(c(mergedFields, list(rowNumbers)), ncol = 4)
            
            # Remove selected groups and add new merged group
            groups <- groups[-selected, , drop=FALSE]
            groups <- rbind(groups, new)
            setGroupsFrom(active, groups)
        }
    })
    
    # Merge selected groups when pressing the button
    observeEvent(input[[id("intersectGroups")]], {
        session$sendCustomMessage(type = "getCheckedBoxes", "intersectGroups")
        sharedData$intersectTime <- TRUE
    })
    
    # Merge selected groups if there are groups to be merged
    observe({
        if (!is.null(input$intersectGroups) && all(input$intersectGroups > 0) &&
            isTRUE(sharedData$intersectTime)) {
            # Set groups to remove to 0 and flag to FALSE
            session$sendCustomMessage(type = "setZero", "intersectGroups")
            sharedData$intersectTime <- FALSE
            
            # Get groups for the data table that is visible and active
            active <- input[[id("dataTypeTab")]]
            groups <- getGroupsFrom(active)
            
            # Create merged group
            selected <- as.numeric(input$intersectGroups)
            mergedFields <- lapply(1:3, function(i)
                paste(groups[selected, i], collapse = " ∩ "))
            rowNumbers <- sort(Reduce(intersect, groups[selected, 4]))
            new <- matrix(c(mergedFields, list(rowNumbers)), ncol = 4)
            
            # Remove selected groups and add new merged group
            groups <- groups[-selected, , drop=FALSE]
            groups <- rbind(groups, new)
            setGroupsFrom(active, groups)
        }
    })
    
    # Render groups interface only if any group exists
    output[[id("groupsList")]] <- renderUI({
        active <- input[[id("dataTypeTab")]]
        groups <- getGroupsFrom(active)
        
        # Don't show anything if there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            list(
                hr(),
                dataTableOutput(id("groupsTable")),
                actionButton(id("mergeGroups"), "Merge"),
                actionButton(id("intersectGroups"), "Intersect"),
                # actionButton(id("complementGroups"), "Complement"),
                # actionButton(id("subtractGroups"), "Subtract"),
                actionButton(id("removeGroups"), "Remove", icon = icon("times"))
            )
        }
    })
}