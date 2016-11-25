#' Group selection interface
#' 
#' @param id Character: identifier of the group selection
#' @param label Character: selectize label
#' @param placeholder Character: selectize placeholder
#' 
#' @importFrom shiny fluidRow column uiOutput selectizeInput actionButton
#' 
#' @return Interface for group selection
selectGroupsUI <- function (
    id, label, placeholder="Click on 'Groups' to create or edit groups") {
    editId <- paste0(id, "Edit")
    modalId <- paste0(id, "Modal")
    groupSelect <- selectizeInput(id, label, choices=NULL, multiple=TRUE, 
                                  width="auto",
                                  options=list(placeholder=placeholder))
    
    if (!is.null(label)) {
        class <- "inline_selectize"
    } else {
        class <- NULL
    }
    
    fluidRow(uiOutput(modalId),
             column(10, groupSelect),
             column(2, actionButton(editId, "Groups", class=class,
                                    class="pull-right", class="btn-info",
                                    style="z-index: 1; position: relative;")))
}

#' Group selection logic
#' 
#' @param session Shiny session
#' @param id Character: identifier of the group selection
#' @param datasetName Character: name of the dataset of interest
#' 
#' @importFrom shinyjs enable disable onclick toggleClass
#' 
#' @return Server logic for group selection
selectGroupsServer <- function(session, id, datasetName) {
    ns <- session$ns
    input <- session$input
    output <- session$output
    
    editId <- paste0(id, "Edit")
    modalId <- paste0(id, "Modal")
    showId <- uId <- paste0(id, "Show")
    uId <- paste0(id, "Call")
    
    output[[modalId]] <- renderUI({
        bsModal2(ns(showId), style="info", trigger=NULL, size=NULL,
                 div(icon("object-group"), "Groups"), 
                 groupsUI(ns(uId), getCategoryData()[[datasetName]]))
    })
    
    # Toggle group selection interface when clicking the "Edit" button
    observeEvent(input[[editId]],
                 toggleModal(session, showId, toggle="open"))
    
    callModule(groupsServer, uId, datasetName)
    
    # Update groups shown in the interface
    observe({
        groups <- names(getGroupsFrom(datasetName))
        if (is.null(groups)) {
            # Disable selection and animate button when clicking disabled input
            groups <- list()
            disable(id)
            onclick(id, runjs(paste0("$('#", ns(editId), 
                                     "').animateCss('rubberBand');")))
        } else {
            enable(id)
            onclick(id, NULL)
        }
        updateSelectizeInput(session, id, choices=groups, selected=groups)
    })
}

#' User interface to group by column
#' 
#' @param ns Namespace function
#' @param dataset Data frame: dataset of interest
#' 
#' @return HTML elements
groupByColumn <- function(ns, dataset) {
    tagList(
        helpText(
            "Automatically create groups according to the unique values of the",
            "selected column. For instance, to create groups by tumour stage,",
            "type", tags$b("tumor_stage"), "and select the first",
            "suggestion that appears."),
        selectizeInput(ns("groupColumn"), "Select column", width="auto", 
                       choices=c("Start typing to search for columns"="", 
                                 names(dataset)))
    )}

#' User interface to group by row
#' 
#' @param ns Namespace function
#' 
#' @return HTML elements
groupByRow <- function(ns) {
    tagList(
        selectizeInput(
            ns("groupRows"), "Row indexes", choices=NULL, multiple=TRUE,
            # Allow to add new items
            width="auto", options=list(
                create=TRUE, createOnBlur=TRUE,
                # Hide discarded user-created items in the dropdown
                persist=FALSE)),
        helpText("Type ", tags$kbd("1:6, 8, 10:19"), "to create a group with",
                 "rows 1 to 6, 8 and 10 to 19.")
    )
}

#' User interface to group by subset expression
#' 
#' @param ns Namespace function
#' 
#' @return HTML elements
groupByExpression <- function(ns) {
    tagList (
        textInput(ns("groupExpression"), "Subset expression", width="auto"),
        helpText('Type ', tags$kbd('X > 8 & Y == "alive"'), ' to select rows',
                 'with values higher than 8 for column X and "alive" for',
                 'column Y.'),
        uiOutput(ns("groupExpressionSuggestions"))
    )
}

#' User interface to group by grep expression
#' 
#' @param ns Namespace function
#' @param dataset Data frame: dataset of interest
#' 
#' @return HTML elements
groupByGrep <- function(ns, dataset) {
    tagList (
        textInput(ns("grepExpression"), "Regular expression", width="auto"),
        selectizeInput(ns("grepColumn"), "Select column to GREP",
                       choices=c("Start typing to search for columns"="",
                                 names(dataset)), width="auto")
    )}

#' Creates UI elements for the grouping feature
#' @param id Character: identifier
#' @param dataset Data frame or matrix: dataset of interest
#' @return HTML elements
groupsUI <- function(id, dataset) {
    ns <- NS(id)
    
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", ns("subsetBy"), sign, what)
    tagList(
        uiOutput(ns("alert")),
        radioButtons(ns("subsetBy"), "Subset by", inline=TRUE,
                     c("Column", "Rows", "Subset expression",
                       "Regular expression")),
        conditionalPanel(checkId("==", "Column"), groupByColumn(ns, dataset)),
        conditionalPanel(checkId("==", "Rows"), groupByRow(ns)),
        conditionalPanel(checkId("==", "Subset expression"), 
                         groupByExpression(ns)),
        conditionalPanel(checkId("==", "Regular expression"), 
                         groupByGrep(ns, dataset)),
        conditionalPanel(checkId("!=", "Column"),
                         textInput(ns("groupName"), "Group name", width="auto",
                                   placeholder="Unnamed")),
        actionButton(ns("createGroup"), "Create group", class ="btn-primary"),
        uiOutput(ns("groupsList"))
    )
}

#' Create groups with the indexes from the unique values of a given column from
#' a dataset
#' 
#' @param col Character: column name
#' @param dataset Matrix or data frame: dataset
#' 
#' @return Named list with the indexes of each unique value from a given column
#' @export
#' 
#' @examples 
#' df <- data.frame(gender=c("male", "female"),
#'                  stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))
#' createGroupByColumn(col="stage", dataset=df)
createGroupByColumn <- function(col, dataset) {
    colData <- as.character(dataset[[col]])
    
    # Replace missing values for "NA" so they are found using the `which` function
    colData[is.na(colData)] <- "NA"
    
    # Create groups according to the chosen column
    groupNames <- sort(unique(colData))
    group <- lapply(lapply(groupNames, `==`, colData), which)
    names(group) <- groupNames
    return(group)
}

#' Create groups from a given string of rows
#' 
#' @param session Shiny session
#' @param rows Character: rows separated by a comma
#' @param dataset Matrix or data frame: dataset
#' 
#' @importFrom shiny tags
#' @return NULL (this function is used to modify the Shiny session's state)
createGroupByRows <- function(session, rows, dataset) {
    # Convert the given string into a sequence of numbers
    rows <- unlist(lapply(rows, function(row) eval(parse(text=row))))
    rows <- sort(unique(rows))
    
    # Remove and warn if selected rows are greater than the rows number
    gtRows <- rows > nrow(dataset)
    if (any(gtRows)) {
        removed <- paste(rows[gtRows], collapse=" ")
        warningAlert(
            session, sum(gtRows), " indexes were above the number of rows ",
            "of the dataset (which is ", nrow(dataset), ").", br(),
            "The following numbers were discarded:", tags$code(removed))
        rows <- rows[!gtRows]
    }
    return(rows)
}

#' Set new groups according to the user input
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' @param dataset Data frame or matrix: dataset of interest
#' @param datasetName Character: name of the dataset
#' 
#' @return Matrix with the group names and respective indexes
createGroupFromInput <- function (session, input, output, dataset,
                                  datasetName) {
    type <- input$subsetBy
    
    if (type == "Column") {
        col <- input$groupColumn
        if (col == "") return(NULL)
        group <- createGroupByColumn(col, dataset)
        group <- cbind(names(group), "Column", col, group)
    } else if (type == "Rows") {
        rows <- input$groupRows
        strRows <- paste(rows, collapse=", ")
        allRows <- createGroupByRows(session, rows, dataset)
        group <- cbind(input$groupName, type, strRows, list(allRows))
    } else if (type == "Subset expression") {
        # Subset dataset using the given expression
        expr <- input$groupExpression
        
        # Test expression before running
        set <- tryCatch(subset(dataset, eval(parse(text=expr))),
                        error=return)
        
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorAlert(session, "Error in the subset expression.",
                       "Check if column names are correct.", br(),
                       "The following error was raised:",
                       tags$code(set$message))
            return(NULL)
        }
        
        rows <- match(rownames(set), rownames(dataset))
        group <- cbind(input$groupName, type, expr, list(rows))
    } else if (type == "Regular expression") {
        # Subset dataset column using given regular expression
        col <- input$grepColumn
        colData <- as.character(dataset[[col]])
        expr <- input$grepExpression
        
        # Test expression before running
        set <- tryCatch(grep(expr, colData), error=return)
        
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorAlert(session, "GREP expression error",
                       "The following error was raised:", br(),
                       tags$code(set$message))
            return(NULL)
        }
        
        strRows <- sprintf('"%s" in %s', expr, col)
        group <- cbind(input$groupName, "GREP", strRows, list(set))
    } 
    # Name group if empty
    if (group[[1]] == "") group[[1]] <- "Unnamed"
    # Standarise rows
    ns <- c("Names", "Subset", "Input", "Rows")
    if (is.matrix(group))
        colnames(group) <- ns
    else
        names(group) <- ns
    rownames(group) <- NULL
    return(group)
}

#' Rename duplicated names from a new group
#' 
#' @note The names of pre-existing groups are not modified.
#' 
#' @param new Matrix: new groups
#' @param old Matrix: pre-existing groups
#' 
#' @return Character with no duplicated group names
renameGroups <- function(new, old) {
    groupNames <- 1
    
    # Rename duplicated names from a new group
    if (!is.null(old) && nrow(old) != 0) {
        newNames <- unlist(new[ , groupNames])
        oldNames <- unlist(old[ , groupNames])
        new[ , groupNames] <- renameDuplicated(newNames, oldNames)
    }
    
    rownames(new) <- new[ , groupNames]
    return(new)
}

#' Set operations on groups
#' 
#' This function can be used on groups to merge, intersect, subtract, etc.
#' 
#' @param input Shiny input
#' @param session Shiny session
#' @param FUN Function: operation to set
#' @param buttonId Character: ID of the button to trigger operation
#' @param symbol Character: operation symbol
#' @param datasetName Character: name of dataset
#' @param sharedData Shiny app's global variable
#' @return NULL (this function is used to modify the Shiny session's state)
operateOnGroups <- function(input, session, FUN, buttonId, symbol=" ",
                            datasetName, sharedData=sharedData) {
    ns <- session$ns
    # Operate on selected groups when pressing the corresponding button
    observeEvent(input[[paste(buttonId, "button", sep="-")]], {
        # Get groups from the dataset
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        
        # Create new set
        new <- NULL
        selected <- input$groupsTable_rows_selected
        if (!identical(FUN, "remove")) {
            mergedFields <- lapply(1:3, function(i) {
                names <- paste(groups[selected, i], collapse=symbol)
                # Add parenthesis around new expression
                names <- paste0("(", names, ")")
                return(names)
            })
            rowNumbers <- sort(as.numeric(Reduce(FUN, groups[selected, 4])))
            new <- matrix(c(mergedFields, list(rowNumbers)), ncol=4)
        }
        
        # Remove selected groups
        if (identical(FUN, "remove"))# || input$removeSetsUsed)
            groups <- groups[-selected, , drop=FALSE]
        
        # Add new groups to top (if there are any)
        if (!is.null(new)) {
            new <- renameGroups(new, groups)
            groups <- rbind(new, groups)
        }
        setGroupsFrom(datasetName, groups)
    })
}

#' Server function for data grouping
#'
#' @inheritParams operateOnGroups
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom shinyjs disabled enable disable
#' @return NULL (this function is used to modify the Shiny session's state)
groupsServer <- function(input, output, session, datasetName) {
    ns <- session$ns
    
    # Update available attributes to suggest in the subset expression
    output$groupExpressionSuggestions <- renderUI({
        if (!is.null(datasetName))
            textSuggestions(ns("groupExpression"),
                            names(getCategoryData()[[datasetName]]))
    })
    
    # Create a new group when clicking on the createGroup button
    observeEvent(input$createGroup, {
        removeAlert(output)
        
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        if (is.null(datasetName)) {
            errorAlert(session, "Data missing", "Load some data first.")
            return(NULL)
        }
        
        new <- createGroupFromInput(session, input, output, 
                                    getCategoryData()[[datasetName]],
                                    datasetName)
        if (!is.null(new)) {
            # Rename duplicated group names
            new <- renameGroups(new, groups)
            
            # Append the new group(s) to the groups already created
            groups <- rbind(new, groups)
            setGroupsFrom(datasetName, groups)
        }
        
        updateSelectizeInput(session, "groupColumn", selected=character())
    })
    
    # Render groups list and show interface to manage groups
    output$groupsTable <- renderDataTable({
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        
        # Show groups only if there is at least one group
        if (!is.null(groups) && nrow(groups) > 0) {
            # Show number of rows for each group
            rows <- lapply(groups[ , 4], length)
            groups[ , 4] <- unlist(rows)
            
            # Ordering the groups (plus safety net for cases with one row)
            ord <- c(1, 4, 2, 3)
            ordered <- groups[ , ord]
            if (!is.matrix(ordered)) {
                ordered <- matrix(ordered, ncol=4)
                colnames(ordered) <- colnames(groups)[ord]
            }
            return(ordered)
        }
    }, style="bootstrap", escape=FALSE, server=TRUE, rownames=FALSE,
    options=list(pageLength=10, lengthChange=FALSE, scrollX=TRUE,
                 ordering=FALSE,
                 # Stack DataTable elements so they fit in the container
                 dom=paste0(
                     '<"row view-filter"<"col-sm-12"<"pull-left"l>',
                     '<"pull-right"f><"clearfix">>>',
                     'rt<"row view-pager"<"col-sm-12"<"text-center"ip>>>')))
    
    # Disable buttons if there's no row selected
    observe({
        if (!is.null(input$groupsTable_rows_selected)) {
            enable("setOperations")
        } else {
            disable("setOperations")
        }
    })
    
    # Remove selected groups
    removeId <- "removeGroups"
    operateOnGroups(input, session, FUN="remove", buttonId=removeId,
                    datasetName=datasetName)
    
    # Merge selected groups
    mergeId <- "mergeGroups"
    operateOnGroups(input, session, FUN=union, buttonId=mergeId, 
                    symbol=" \u222A ", datasetName=datasetName)
    
    # Intersect selected groups
    intersectId <- "intersectGroups"
    operateOnGroups(input, session, FUN=intersect, datasetName=datasetName,
                    buttonId=intersectId, symbol=" \u2229 ")
    
    # Render groups interface only if at least one group exists
    output$groupsList <- renderUI({
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        
        operationButton <- function(operation, operationId, ...) {
            actionButton(paste(operationId, "button", sep="-"), operation, ...)
        }
        
        # Don't show anything when there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            operations <- div(id=ns("setOperations"), class="btn-group",
                              operationButton("Merge", ns(mergeId)),
                              operationButton("Intersect", ns(intersectId)),
                              # actionButton("complementGroups", ns("Complement")),
                              # actionButton("subtractGroups", ns("Subtract")),
                              operationButton("Remove", ns(removeId),
                                              icon=icon("times")))
            tagList(
                hr(),
                dataTableOutput(ns("groupsTable")),
                helpText("Select groups by clicking on them in order to merge,",
                         "intersect or remove selected groups."),
                disabled(operations),
                actionButton(ns("removeAll"), class="btn-danger", 
                             "Remove all groups")#,
                #checkboxInput(ns("removeSetsUsed"), "Remove original groups",
                #              value=TRUE)
            )
        }
    })
    
    observeEvent(input$removeAll, setGroupsFrom(datasetName, NULL))
}

# attr(groupsUI, "loader") <- "data"
# attr(groupsServer, "loader") <- "data"