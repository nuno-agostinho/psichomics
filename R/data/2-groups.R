## TODO(NunoA): help the user by suggesting column names; try this with JavaScript

## TODO(NunoA): Allow to select the checkboxes in groups from different "pages"
## in the table; maybe by memorising the ones selected in each page? Then always 
## show the ones selected as checked... try to do it by memorising when clicking 
## to go on next page or something and clean checkboxes when 
## merging/intersect/removing groups

## TODO(NunoA): Easily identify multiple operations by isolating each with
## parenthesis; e.g. using A + B ∩ C + D, we don't know if the original
## operations were either A + (B ∩ C) + D or (A + B) ∩ (C + D) or ...

## TODO(NunoA): build GREP expressions
## TODO(NunoA): isolate GREP expression to give a warning on failure

name <- "Groups"

groupByColumn <- function() { list(
    selectizeInput(id("groupColumn"), "Select column", choices = NULL,
                   options = list(
                       placeholder = "Start typing to search for columns")),
    helpText("Groups will be created automatically depending on the",
             "given column.")
)}

groupByRow <- function() { list(
    selectizeInput(
        id("groupRows"), "Select rows", choices = NULL,
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

groupByExpression <- function() { list (
    textInput(id("groupExpression"), "Subset expression"),
    helpText('Type ', tags$kbd('X > 8 & Y == "alive"'), ' to select rows with',
             'values higher than 8 for column X and "alive" for column Y.'),
    uiOutput(id("groupExpressionAutocomplete"))
)}

groupByGrep <- function() { list (
    textInput(id("grepExpression"), "GREP expression"),
    selectizeInput(id("grepColumn"),
                   "Select column to GREP",
                   choices = NULL,
                   options = list(
                       placeholder = "Start typing to search for columns"))
)}

#' Creates UI elements for the grouping feature
ui <- function() {
    checkId <- function (sign, what)
        sprintf("input[id='%s'] %s '%s'", id("subsetBy"), sign, what)
    
    list(
        radioButtons(id("subsetBy"), "Subset by", inline = TRUE,
                     c("Column", "Rows", "Expression", "Grep")),
        conditionalPanel(checkId("==", "Column"), groupByColumn()),
        conditionalPanel(checkId("==", "Rows"), groupByRow()),
        conditionalPanel(checkId("==", "Expression"), groupByExpression()),
        conditionalPanel(checkId("==", "Grep"), groupByGrep()),
        conditionalPanel(checkId("!=", "Column"),
                         textInput(id("groupName"), "Group name",
                                   placeholder = "Unnamed")),
        actionButton(id("createGroup"), "Create group", class ="btn-primary"),
        uiOutput(id("groupsList"))
    )
}

#' Set new groups according to the user input
#' 
#' The groups are inserted in a matrix
createGroupFromInput <- function (input, output, session) {
    active <- input[[id("datasetTab")]]
    if (is.null(active)) {
        errorModal(session, "Data missing", "Load some data first.")
        return(NULL)
    }
    
    type <- input[[id("subsetBy")]]
    
    columnInput <- input[[id("groupColumn")]]
    data <- getCategoryData()[[active]]
    
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
        rows <- input[[id("groupRows")]]
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
        group <- cbind(input[[id("groupName")]], type, strRows, list(rows))
    } else if (type == "Expression") {
        # Subset data using the given expression
        expr <- input[[id("groupExpression")]]
        
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
        group <- cbind(input[[id("groupName")]], type, expr, list(rows))
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

operateOnGroups <- function(input, session, sharedData, FUN, name, 
                            symbol = " ") {
    # Operate on selected groups when pressing the corresponding button
    observeEvent(input[[paste(name, "Button", sep = "_")]], {
        session$sendCustomMessage(type = "getCheckedBoxes", name)
        sharedData[[name]] <- TRUE
    })
    
    observe({
        if (!is.null(input[[name]]) && all(input[[name]] > 0) &&
            isTRUE(sharedData[[name]])) {
            # Set operation groups as 0 and flag to FALSE
            session$sendCustomMessage(type = "setZero", name)
            sharedData[[name]] <- FALSE
            
            # Get groups for the data table that is visible and active
            active <- input[[id("datasetTab")]]
            groups <- getGroupsFrom(active)
            
            # Create new set
            new <- NULL
            selected <- as.numeric(input[[name]])
            if (!identical(FUN, "remove")) {
                mergedFields <- lapply(1:3, function(i)
                    paste(groups[selected, i], collapse = symbol))
                rowNumbers <- sort(as.numeric(Reduce(FUN, groups[selected, 4])))
                new <- matrix(c(mergedFields, list(rowNumbers)), ncol = 4)
            }
            
            # Remove selected groups
            if (identical(FUN, "remove") || input[[id("removeSetsUsed")]])
                groups <- groups[-selected, , drop=FALSE]
            
            # Add new groups to top (if there are any)
            if (!is.null(new)) {
                new <- renameGroups(new, groups)
                groups <- rbind(new, groups)
            }
            setGroupsFrom(active, groups)
        }
    })
}

server <- function(input, output, session) {
    # Update available attributes to suggest in the group expression
    output[[id("groupExpressionAutocomplete")]] <- renderUI({
        active <- input[[id("datasetTab")]]
        attributes <- names(getCategoryData()[[active]])
        textComplete(id("groupExpression"), attributes)
    })
    
    # Update columns available for creating groups when there's loaded data
    observeEvent(input[[id("datasetTab")]], {
        active <- input[[id("datasetTab")]]
        for (i in id(c("groupColumn", "grepColumn"))) {
            updateSelectizeInput(session, i, selected = NULL,
                                 choices = names(getCategoryData()[[active]]))
        }
    })
    
    # Create a new group when clicking on the createGroup button
    observeEvent(input[[id("createGroup")]], {
        # Get groups for the data table that is visible and active
        active <- input[[id("datasetTab")]]
        groups <- getGroupsFrom(active)
        
        # Create new group(s) from user input
        new <- createGroupFromInput(input, output, session)
        
        if (!is.null(new)) {
            # Rename duplicated group names
            new <- renameGroups(new, groups)
            
            # Append the new group(s) to the groups already created
            groups <- rbind(new, groups)
            setGroupsFrom(active, groups)
        }
    })
    
    # Render groups list and show interface to manage groups
    output[[id("groupsTable")]] <- renderDataTable({
        ## TODO(NunoA): Allow to remove and merge selected rows from the groups
        ## This could be done using checkboxes; how to retrieve which checkboxes
        ## were checked? Possible with data table or javascript?
        active <- input[[id("datasetTab")]]
        groups <- getGroupsFrom(active)
        
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
                      ordering = FALSE),
    escape = FALSE)
    
    # Remove selected groups
    removeId <- id("removeGroups")
    operateOnGroups(input, session, sharedData, "remove", removeId)
    
    # Merge selected groups
    mergeId <- id("mergeGroups")
    operateOnGroups(input, session, sharedData, union, mergeId, " + ")
    
    # Intersect selected groups
    intersectId <- id("intersectGroups")
    operateOnGroups(input, session, sharedData, intersect, intersectId, " ∩ ")
    
    # Render groups interface only if any group exists
    output[[id("groupsList")]] <- renderUI({
        active <- input[[id("datasetTab")]]
        groups <- getGroupsFrom(active)
        
        operationButton <- function(operation, operationId, ...)
            actionButton(paste(operationId, "Button", sep = "_"),
                         operation, ...)
        
        # Don't show anything when there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            list(
                hr(),
                dataTableOutput(id("groupsTable")),
                div(class="btn-group",
                    operationButton("Merge", mergeId),
                    operationButton("Intersect", intersectId),
                    # actionButton(id("complementGroups"), "Complement"),
                    # actionButton(id("subtractGroups"), "Subtract"),
                    operationButton("Remove", removeId, icon = icon("times"))),
                checkboxInput(id("removeSetsUsed"), "Remove original groups", 
                              value = TRUE)
            )
        }
    })
}