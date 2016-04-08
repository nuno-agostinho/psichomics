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
             'values higher than 8 for column X and "alive" for column Y.')
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
        selectizeInput(id("subsetBy"), "Subset by",
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
    ns <- c("Names", "Subset", "Input", "Rows")
    if (is.matrix(group))
        colnames(group) <- ns
    else
        names(group) <- ns
    rownames(groups) <- NULL
    return(group)
}

server <- function(input, output, session) {
    # Update columns available for creating groups when there's loaded data
    observeEvent(input[[id("dataTypeTab")]], {
        active <- input[[id("dataTypeTab")]]
        for (i in id(c("groupColumn", "grepColumn"))) {
            updateSelectizeInput(session, i,
                                 selected = NULL,
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
        old <- groups
        groups <- rbind(new, groups)
        
        #Rename duplicated group names
        newNames <- unlist(new[ , "Names"])
        oldNames <- unlist(old[ , "Names"])
        groups[ , "Names"] <- renameDuplicated(newNames, oldNames)
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
            
            # Remove selected groups
            if (input[[id("removeSetsUsed")]])
                groups <- groups[-selected, , drop=FALSE]
            
            # Add new merged group
            groups <- rbind(new, groups)
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
            
            # Remove selected groups
            if (input[[id("removeSetsUsed")]])
                groups <- groups[-selected, , drop=FALSE]
            
            # Add new merged group
            groups <- rbind(new, groups)
            setGroupsFrom(active, groups)
        }
    })
    
    # Render groups interface only if any group exists
    output[[id("groupsList")]] <- renderUI({
        active <- input[[id("dataTypeTab")]]
        groups <- getGroupsFrom(active)
        
        # Don't show anything when there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            list(
                hr(),
                dataTableOutput(id("groupsTable")),
                div(class="btn-group",
                    actionButton(id("mergeGroups"), "Merge"),
                    actionButton(id("intersectGroups"), "Intersect"),
                    # actionButton(id("complementGroups"), "Complement"),
                    # actionButton(id("subtractGroups"), "Subtract"),
                    actionButton(id("removeGroups"), "Remove",
                                 icon = icon("times"))),
                checkboxInput(id("removeSetsUsed"), "Remove original groups", 
                              value = TRUE)
            )
        }
    })
}