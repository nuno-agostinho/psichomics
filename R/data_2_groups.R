## TODO(NunoA): Allow to select the checkboxes in groups from different "pages"
## in the table; maybe by memorising the ones selected in each page? Then always 
## show the ones selected as checked... try to do it by memorising when clicking 
## to go on next page or something and clean checkboxes when 
## merging/intersect/removing groups

#' User interface to group by column
groupByColumn <- function(ns) {
    list(
        selectizeInput(ns("groupColumn"), "Select column", choices = NULL,
                       options = list(
                           placeholder = "Start typing to search for columns")),
        helpText("Groups will be created automatically depending on the",
                 "given column.")
    )}

#' User interface to group by row
groupByRow <- function(ns) { list(
    selectizeInput(
        ns("groupRows"), "Select rows", choices = NULL,
        multiple = TRUE,
        # Allow to add new items
        options = list(
            create = TRUE, createOnBlur=TRUE,
            ## TODO(NunoA): only allow numbers (use selectize.js REGEX option)
            # Hide discarded user-created items in the dropdown
            persist = FALSE)),
    helpText("Select rows like in R. Type ", tags$kbd("1:6, 8, 10:19"),
             "to create a group with rows 1 to 6, 8 and 10 to 19.")
)}

#' User interface to group by subset expression
groupByExpression <- function(ns) { list (
    textInput(ns("groupExpression"), "Subset expression"),
    helpText('Type ', tags$kbd('X > 8 & Y == "alive"'), ' to select rows with',
             'values higher than 8 for column X and "alive" for column Y.'),
    uiOutput(ns("groupExpressionAutocomplete"))
)}

#' User interface to group by grep expression
groupByGrep <- function(ns) { list (
    textInput(ns("grepExpression"), "GREP expression"),
    selectizeInput(ns("grepColumn"),
                   "Select column to GREP",
                   choices = NULL,
                   options = list(
                       placeholder = "Start typing to search for columns"))
)}

#' Creates UI elements for the grouping feature
groupsUI <- function(id, tab) {
    ns <- NS(id)
    
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", ns("subsetBy"), sign, what)
    tab("Groups", br(),
        uiOutput(ns("modal")),
        radioButtons(ns("subsetBy"), "Subset by", inline = TRUE,
                     c("Column", "Rows", "Expression", "Grep")),
        conditionalPanel(checkId("==", "Column"), groupByColumn(ns)),
        conditionalPanel(checkId("==", "Rows"), groupByRow(ns)),
        conditionalPanel(checkId("==", "Expression"), groupByExpression(ns)),
        conditionalPanel(checkId("==", "Grep"), groupByGrep(ns)),
        conditionalPanel(checkId("!=", "Column"),
                         textInput(ns("groupName"), "Group name",
                                   placeholder = "Unnamed")),
        actionButton(ns("createGroup"), "Create group", class ="btn-primary"),
        uiOutput(ns("groupsList"))
    )
}

#' Set new groups according to the user input
#' 
#' The groups are inserted in a matrix
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
createGroupFromInput <- function (input, output, session) {
    if (is.null(getActiveDataset())) {
        errorModal(session, "Data missing", "Load some data first.")
        return(NULL)
    }
    
    type <- input$subsetBy
    
    columnInput <- input$groupColumn
    data <- getCategoryData()[[getActiveDataset()]]
    
    if (type == "Column") {
        colData <- data[[columnInput]]
        
        # Replace NAs for "NA" so they can be find using the `which` function
        colData[is.na(colData)] <- "NA"
        
        # Create groups according to the chosen column
        groupNames <- sort(unique(colData))
        rows <- lapply(lapply(groupNames, `==`, colData), which)
        group <- cbind(groupNames, type, columnInput, rows)
    } else if (type == "Rows") {
        # Convert the given string into a sequence of numbers
        rows <- input$groupRows
        strRows <- paste(rows, collapse = ", ")
        rows <- unlist(lapply(rows, function(row) eval(parse(text = row))))
        rows <- sort(unique(rows))
        
        # Remove and warn if selected rows are greater than the rows number
        gtRows <- rows > nrow(data)
        if (any(gtRows)) {
            removed <- paste(rows[gtRows], collapse = " ")
            warningModal(
                session, "Selected rows don't exist",
                paste0(sum(gtRows), " numbers were above the number of rows ",
                       "of the active dataset (which is", nrow(data), ")."), 
                br(), br(), "The following numbers were discarded:", 
                tags$code(removed))
            rows <- rows[!gtRows]
        }
        group <- cbind(input$groupName, type, strRows, list(rows))
    } else if (type == "Expression") {
        # Subset data using the given expression
        expr <- input$groupExpression
        
        # Test expression before running
        set <- tryCatch(subset(data, eval(parse(text = expr))),
                        error = return)
        
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorModal(session, "Subset expression error",
                       "Check if column names are correct.", br(), br(), 
                       "The following error was raised:", br(),
                       tags$code(set$message))
            return(NULL)
        }
        
        rows <- match(rownames(set), rownames(data))
        group <- cbind(input$groupName, type, expr, list(rows))
    } else if (type == "Grep") {
        ## TODO(NunoA): Subset data with the GREP expression for the given column
        group <- rep(NA, 4)
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
#' @param sharedData Shiny app's global variable
#' @param FUN Function: operation to set
#' @param buttonId Character: ID of the button to trigger operation
#' @param symbol Character: operation symbol
operateOnGroups <- function(input, session, sharedData, FUN, buttonId, 
                            symbol=" ") {
    ns <- session$ns
    # Operate on selected groups when pressing the corresponding button
    observeEvent(input[[paste(buttonId, "button", sep="-")]], {
        session$sendCustomMessage(type="getCheckedBoxes", "selectedGroups")
        sharedData$javascriptSent <- TRUE
        sharedData$groupsFUN <- FUN
        sharedData$groupSymbol <- symbol
        # appServer saves the result to the R variable sharedData$selectedGroups
    })
}

#' Server function for data grouping
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
groupsServer <- function(input, output, session) {
    ns <- session$ns
    
    # Update available attributes to suggest in the group expression
    output$groupExpressionAutocomplete <- renderUI({
        if (!is.null(getActiveDataset())) {
            attributes <- names(getCategoryData()[[getActiveDataset()]])
            textComplete(ns("groupExpression"), attributes)
        }
    })
    
    # Update columns available for creating groups when there's loaded data
    observeEvent(getActiveDataset(), {
        for (i in c("groupColumn", "grepColumn")) {
            updateSelectizeInput(session, i, selected = NULL,
                                 choices = names(getCategoryData()[[getActiveDataset()]]))
        }
    })
    
    # Create a new group when clicking on the createGroup button
    observeEvent(input$createGroup, {
        # Get groups for the data table that is visible and active
        groups <- getGroupsFrom(getActiveDataset())
        
        # Create new group(s) from user input
        new <- createGroupFromInput(input, output, session)
        
        if (!is.null(new)) {
            # Rename duplicated group names
            new <- renameGroups(new, groups)
            
            # Append the new group(s) to the groups already created
            groups <- rbind(new, groups)
            setGroupsFrom(getActiveDataset(), groups)
        }
    })
    
    # Render groups list and show interface to manage groups
    output$groupsTable <- renderDataTable({
        groups <- getGroupsFrom(getActiveDataset())
        
        # Show groups only if there is at least one group
        if (!is.null(groups) && nrow(groups) > 0) {
            # Show number of rows for each group
            rows <- lapply(groups[ , 4], length)
            groups[ , 4] <- unlist(rows)
            
            # Ordering the groups (plus safety net for cases with one row)
            ord <- c(1, 4, 2, 3)
            ordered <- groups[ , ord]
            if (!is.matrix(ordered)) {
                ordered <- matrix(ordered, ncol = 4)
                colnames(ordered) <- colnames(groups)[ord]
            }
            
            # Add checkboxes
            pick <- paste("<center><input number=", 1:nrow(ordered),
                          "name='checkGroups' type='checkbox'/></center>")
            res <- cbind(pick, ordered)
            colnames(res)[1] <- "<input name='checkAllGroups' type='checkbox'/>"
            return(res)
        }
    }, options = list(pageLength = 10, lengthChange = FALSE, scrollX = TRUE,
                      #filter = FALSE, info = FALSE, paginationType = "simple",
                      ordering = FALSE,
                      # Stack DataTable elements so they fit in the container
                      dom = paste0(
                          '<"row view-filter"<"col-sm-12"<"pull-left"l>',
                          '<"pull-right"f><"clearfix">>>',
                          'rt<"row view-pager"<"col-sm-12"<"text-center"ip>>>')),
    escape = FALSE)
    
    # Remove selected groups
    removeId <- "removeGroups"
    operateOnGroups(input, session, sharedData, FUN="remove", buttonId=removeId)
    
    # Merge selected groups
    mergeId <- "mergeGroups"
    operateOnGroups(input, session, sharedData, FUN=union, buttonId=mergeId,
                    symbol=" \u222A ")
    
    # Intersect selected groups
    intersectId <- "intersectGroups"
    operateOnGroups(input, session, sharedData, FUN=intersect,
                    buttonId=intersectId, symbol=" \u2229 ")
    
    
    observe({
        if (!is.null(sharedData$selectedGroups) &&
            all(sharedData$selectedGroups > 0) &&
            isTRUE(sharedData$javascriptSent) &&
            isTRUE(sharedData$javascriptRead)) {
            
            FUN <- sharedData$groupsFUN
            symbol <- sharedData$groupSymbol
            
            # Set operation groups as 0 and flag to FALSE
            session$sendCustomMessage(type = "setZero", "selectedGroups")
            sharedData$javascriptSent <- FALSE
            sharedData$javascriptRead <- FALSE
            
            # Get groups for the data table that is visible and active
            groups <- getGroupsFrom(getActiveDataset())
            
            # Create new set
            new <- NULL
            selected <- as.numeric(sharedData$selectedGroups)
            if (!identical(FUN, "remove")) {
                mergedFields <- lapply(1:3, function(i) {
                    names <- paste(groups[selected, i], collapse = symbol)
                    # Add parenthesis around new expression
                    names <- paste0("(", names, ")")
                    return(names)
                })
                rowNumbers <- sort(as.numeric(Reduce(FUN, groups[selected, 4])))
                new <- matrix(c(mergedFields, list(rowNumbers)), ncol = 4)
            }
            
            # Remove selected groups
            if (identical(FUN, "remove") || input$removeSetsUsed)
                groups <- groups[-selected, , drop=FALSE]
            
            # Add new groups to top (if there are any)
            if (!is.null(new)) {
                new <- renameGroups(new, groups)
                groups <- rbind(new, groups)
            }
            setGroupsFrom(getActiveDataset(), groups)
        }
    })
    
    
    
    
    # Render groups interface only if any group exists
    output$groupsList <- renderUI({
        groups <- getGroupsFrom(getActiveDataset())
        
        operationButton <- function(operation, operationId, ...) {
            actionButton(paste(operationId, "button", sep="-"), operation, ...)
        }
        
        # Don't show anything when there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            tagList(
                hr(),
                dataTableOutput(ns("groupsTable")),
                div(class="btn-group",
                    operationButton("Merge", ns(mergeId)),
                    operationButton("Intersect", ns(intersectId)),
                    # actionButton("complementGroups", ns("Complement")),
                    # actionButton("subtractGroups", ns("Subtract")),
                    operationButton("Remove", ns(removeId), icon = icon("times"))),
                checkboxInput(ns("removeSetsUsed"), "Remove original groups",
                              value = TRUE)
            )
        }
    })
}

attr(groupsUI, "loader") <- "data"
attr(groupsServer, "loader") <- "data"