#' Group selection interface
#' 
#' @param id Character: identifier of the group selection
#' @param label Character: selectize label
#' @param placeholder Character: selectize placeholder
#' @param noGroupsLabel Character: label to show when no groups may be selected
#' (if NULL, the option to show no groups will not be shown)
#' @param groupsLabel Character: label to show to the option of using groups
#' when no groups may be selected
#' 
#' @importFrom shiny fluidRow column uiOutput selectizeInput actionButton
#' radioButtons actionLink
#' 
#' @note To allow the user to (explicitly) select no groups, pass the 
#' \code{noGroupsLabel} and \code{groupsLabel} arguments.
#' 
#' @seealso selectGroupsServer getSelectedGroups
#' 
#' @return Interface for group selection
selectGroupsUI <- function (
    id, label, placeholder="Click on 'Groups' to create or edit groups",
    noGroupsLabel=NULL, groupsLabel=NULL) {
    editId <- paste0(id, "Edit")
    modalId <- paste0(id, "Modal")
    groupSelect <- selectizeInput(id, label, choices=NULL, multiple=TRUE, 
                                  width="auto",
                                  options=list(placeholder=placeholder,
                                               plugins=list('remove_button',
                                                            'drag_drop')))
    
    if ( !is.null(label) ) {
        if ( is.null(noGroupsLabel) ) {
            label <- column(12, groupSelect$children[[1]])
        } else {
            # Use label in radio buttons instead
            radioLabel <- groupSelect$children[[1]]
            label <- NULL
        }
        groupSelect$children[[1]] <- NULL
    }
    
    select <- fluidRow(
        uiOutput(modalId), label,
        column(10, groupSelect),
        column(2, actionButton(
            editId, "Groups", class="pull-right", class="btn-info",
            style="z-index: 1; position: relative;")))
    
    if ( !is.null(noGroupsLabel) ) {
        noGroupsId <- paste0(id, "Selection")
        
        choices <- c("noGroups", "groups")
        names(choices) <- c(noGroupsLabel, groupsLabel)
        
        select <- tagList(
            radioButtons(noGroupsId, radioLabel, choices=choices),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", noGroupsId, "groups"),
                select))
    }
    return(select)
}

#' Group selection logic
#' 
#' @param session Shiny session
#' @param id Character: identifier of the group selection
#' 
#' @importFrom shinyjs enable disable onclick toggleClass
#' 
#' @return Server logic for group selection
selectGroupsServer <- function(session, id) {
    datasetName <- "Clinical data"
    ns <- session$ns
    input <- session$input
    output <- session$output
    
    editId <- paste0(id, "Edit")
    modalId <- paste0(id, "Modal")
    showId <- uId <- paste0(id, "Show")
    uId <- paste0(id, "Call")
    
    output[[modalId]] <- renderUI({
        bsModal2(ns(showId), style="info", trigger=NULL, size=NULL,
                 div(icon("object-group"), "Groups"), groupsUI(ns(uId)))
    })
    
    # Toggle group selection interface when clicking the "Edit" button
    observeEvent(input[[editId]],
                 toggleModal(session, showId, toggle="open"))
    
    callModule(groupsServer, uId, datasetName)
    
    # Update groups shown in the interface
    observe({
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        groups <- rownames(groups)
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

#' Creates UI elements for the grouping feature
#' 
#' @param id Character: identifier
#' 
#' @return HTML elements
groupsUI <- function(id) {
    ns <- NS(id)
    
    groupOptions <- function(id, loaded) {
        if (id == "Patients") {
            title <- "Group by patients"
            choices <- getPatientId()
            dataset <- getClinicalData()
            example <- tagList(
                "For instance, to create groups by tumour stage, type",
                tags$b("tumor_stage"), "and select the first suggestion.")
        } else if (id == "Samples") {
            title <- "Group by samples"
            choices <- getSampleId()
            dataset <- getSampleInfo()
            example <- NULL
        }
        
        modalId <- gsub("(.*)Call\\-.*", "\\1Show", ns(id))
        missingData <- function(dataset, message, linkText)
            div(class="alert alert-danger", role="alert",
                icon("exclamation-circle"), message,
                tags$a(linkText, onclick=loadRequiredData( modalId )))
        
        if (loaded) {
            suggestedCols <- attr(dataset, "show")
            if (!is.null(suggestedCols)) {
                cols <- colnames(dataset)
                suggestedIndex <- match(suggestedCols, cols)
                suggestedIndex <- suggestedIndex[!is.na(suggestedIndex)]
                cols <- list("Start typing to search for attributes"="",
                             "Suggested attributes"=cols[suggestedIndex],
                             "Other attributes"=cols[-suggestedIndex])
            } else {
                cols <- colnames(dataset)
                cols <- c("Start typing to search for columns"="", cols)
            }
            
            navbarMenu(
                title,
                tabPanel("Attribute",
                         groupByAttribute(ns, cols, id, example)),
                tabPanel("Index/Identifier", groupById(ns, id, choices)),
                "----", "Advanced options",
                tabPanel("Subset expression", groupByExpression(ns, id)),
                tabPanel("Regular expression", groupByGrep(ns, cols, id)))
        } else if ( id == "Patients" ) {
            tabPanel(title, missingData(
                "Clinical data",
                "No clinical data loaded to group by patients.",
                "Please, load clinical data."))
        } else if ( id == "Samples" ) {
            tabPanel(title, missingData(
                "Inclusion levels",
                "No sample information loaded to group by samples.",
                "Please, load sample information."))
        }
    }
    
    tagList(
        uiOutput(ns("alert")),
        tabsetPanel(id=ns("groupBy"), type="pills",
                    groupOptions("Patients", !is.null(getPatientId()) ),
                    groupOptions("Samples", !is.null(getSampleId()) )),
        uiOutput(ns("groupsList")))
}

#' User interface to group by attribute
#' 
#' @param ns Namespace function
#' @param cols Character or list: name of columns to show
#' @param id Character: identifier
#' @param example Character: text to show as an example
#' 
#' @return HTML elements
groupByAttribute <- function(ns, cols, id, example) {
    if (!is.null(example)) example <- tagList(" ", example)
    
    tagList(
        helpText("Automatically create groups according to the unique values",
                 "for the selected attribute.", example),
        selectizeInput(ns(paste0("groupAttribute", id)), "Select attribute",
                       width="auto", choices=cols),
        actionButton(ns(paste0("createGroupAttribute", id)), "Create group",
                     class ="btn-primary")
    )}

#' User interface to group by row
#' 
#' @inheritParams groupByAttribute
#' @param choices Character: identifier suggestions
#' 
#' @importFrom shiny textInput
#' 
#' @return HTML elements
groupById <- function(ns, id, choices) {
    tagList(
        selectizeInput(
            ns(paste0("groupRows", id)), paste(id, "indexes or identifiers"),
            choices=choices, multiple=TRUE,
            # Allow to add new items
            width="auto", options=list(
                create=TRUE, createOnBlur=TRUE,
                # Hide discarded user-created items in the dropdown
                persist=FALSE)),
        helpText("Type ", tags$kbd("1:6, 8, 10:19"), "to create a group with",
                 "rows 1 to 6, 8 and 10 to 19."),
        textInput(ns(paste0("groupNameRows", id)), "Group name", width="auto",
                  placeholder="Unnamed"),
        actionButton(ns(paste0("createGroupRows", id)), "Create group", 
                     class="btn-primary")
    )
}

#' User interface to group by subset expression
#' 
#' @inheritParams groupByAttribute
#' 
#' @importFrom shiny textInput
#' 
#' @return HTML elements
groupByExpression <- function(ns, id) {
    tagList (
        textInput(ns(paste0("groupExpression", id)), "Subset expression",
                  width="auto"),
        helpText('Type ', tags$kbd('X > 8 & Y == "alive"'), ' to select rows',
                 'with values higher than 8 for column X and "alive" for',
                 'column Y.'),
        uiOutput(ns(paste0("groupExpressionSuggestions", id))),
        textInput(ns(paste0("groupNameSubset", id)), "Group name", width="auto",
                  placeholder="Unnamed"),
        actionButton(ns(paste0("createGroupSubset", id)), "Create group",
                     class="btn-primary")
    )
}

#' User interface to group by grep expression
#' 
#' @inheritParams groupByAttribute
#' 
#' @importFrom shiny textInput
#' 
#' @return HTML elements
groupByGrep <- function(ns, cols, id) {
    tagList (
        textInput(ns(paste0("grepExpression", id)), "Regular expression",
                  width="auto"),
        selectizeInput(ns(paste0("grepColumn", id)), "Select column to GREP",
                       choices=cols, width="auto"),
        textInput(ns(paste0("groupNameRegex", id)), "Group name", width="auto",
                  placeholder="Unnamed"),
        actionButton(ns(paste0("createGroupRegex", id)), "Create group", 
                     class="btn-primary")
    )}

#' Prepare to create group according to specific details
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' @param id Character: identifier of the group selection
#' @param type Character: type of group to create
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
createGroup <- function(session, input, output, id, type) {
    removeAlert(output)
    
    if (id == "Patients")
        dataset <- getClinicalData()
    else if (id == "Samples")
        dataset <- getSampleInfo()
    new <- createGroupFromInput(session, input, output, dataset, id, type)
    
    if (!is.null(new)) {
        # Rename duplicated group names
        groups <- getGroupsFrom("Clinical data", complete=TRUE)
        new <- renameGroups(new, groups)
        
        # Append the new group(s) to the groups already created
        groups <- rbind(new, groups)
        setGroupsFrom("Clinical data", groups)
    }
    updateSelectizeInput(session, paste0("groupAttribute", id),
                         selected=character())
}

#' Set new groups according to the user input
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' @param dataset Data frame or matrix: dataset of interest
#' @param id Character: identifier of the group selection
#' @param type Character: type of group to create
#' 
#' @return Matrix with the group names and respective elements
createGroupFromInput <- function (session, input, output, dataset, id, type) {
    if (type == "Attribute") {
        col <- input[[paste0("groupAttribute", id)]]
        if (col == "") return(NULL)
        group <- createGroupByAttribute(col, dataset)
        group <- cbind(names(group), type, col, group)
    } else if (type == "Index/Identifier") {
        rows <- input[[paste0("groupRows", id)]]
        strRows <- paste(rows, collapse=", ")
        
        identifiers <- switch(id, "Patients"=getPatientId(),
                              "Samples"=getSampleId())
        
        allRows <- createGroupById(session, rows, dataset, identifiers)
        group <- cbind(input[[paste0("groupNameRows", id)]], type, strRows,
                       list(allRows))
    } else if (type == "Subset") {
        # Subset dataset using the given expression
        expr <- input[[paste0("groupExpression", id)]]
        # Test expression before running
        set <- tryCatch(subset(dataset, eval(parse(text=expr))), error=return)
        
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorAlert(session, "Error in the subset expression.",
                       "Check if column names are correct.", br(),
                       "The following error was raised:",
                       tags$code(set$message))
            return(NULL)
        }
        
        rows <- match(rownames(set), rownames(dataset))
        rows <- rownames(dataset)[rows]
        group <- cbind(input[[paste0("groupNameSubset", id)]], type, expr,
                       list(rows))
    } else if (type == "Regex") {
        # Subset dataset column using given regular expression
        col <- input[[paste0("grepColumn", id)]]
        colData <- as.character(dataset[[col]])
        expr <- input[[paste0("grepExpression", id)]]
        
        # Test expression before running
        set <- tryCatch(grep(expr, colData), error=return)
        
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorAlert(session, "GREP expression error",
                       "The following error was raised:", br(),
                       tags$code(set$message))
            return(NULL)
        }
        
        set <- rownames(dataset)[set]
        strRows <- sprintf('"%s" in %s', expr, col)
        group <- cbind(input[[paste0("groupNameRegex", id)]], "GREP", strRows, 
                       list(set))
    }
    
    # Name group if empty
    if (group[[1]] == "") group[[1]] <- "Unnamed"
    
    # Standardise rows
    ns <- c("Names", "Subset", "Input", id)
    if (is.matrix(group))
        colnames(group) <- ns
    else
        names(group) <- ns
    rownames(group) <- NULL
    
    clinical <- getClinicalData()
    samples <- getSampleId()
    match <- getClinicalMatchFrom("Inclusion levels")
    
    # Replace sample indexes with sample names
    if (id == "Samples") {
        group[ , "Samples"] <- lapply(group[ , "Samples"],
                                      function(i) samples[i])
    }
    
    # Match patients with samples (or vice-versa)
    if (!is.null(clinical) && !is.null(samples) && !is.null(match)) {
        if (id == "Patients") {
            patients <- group[ , "Patients"]
            samples <- getMatchingSamples(patients, samples, clinical,
                                          match=match)
            group <- cbind(group, "Samples"=samples)
        } else if (id == "Samples") {
            samples <- group[ , "Samples"]
            patients <- lapply(group[ , "Samples"], function(i) {
                m <- match[i]
                return(unique(m[!is.na(m)]))
            })
            group <- cbind(group, "Patients"=patients)
            group <- group[ , c(1:3, 5, 4), drop=FALSE]
        }
    }
    return(group)
}

#' @inherit createGroupByAttribute
#' @export
createGroupByColumn <- function(col, dataset) {
    .Deprecated("createGroupByAttribute")
    createGroupByAttribute(col, dataset)
}

#' Create groups based on the unique values of a given column
#' 
#' @param col Character: column name
#' @param dataset Matrix or data frame: dataset
#' 
#' @return Named list with each unique value from a given column and respective
#' elements
#' @export
#' 
#' @examples 
#' df <- data.frame(gender=c("male", "female"),
#'                  stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))
#' rownames(df) <- paste0("patient-", LETTERS[1:8])
#' createGroupByAttribute(col="stage", dataset=df)
createGroupByAttribute <- function(col, dataset) {
    colData <- as.character(dataset[[col]])
    names(colData) <- rownames(dataset)
    
    # Replace missing values for "NA" so they are found using the `which` function
    colData[is.na(colData)] <- "NA"
    
    # Create groups according to the chosen column
    groupNames <- sort(unique(colData))
    group <- lapply(lapply(groupNames, `==`, colData), which)
    names(group) <- groupNames
    group <- sapply(group, names)
    return(group)
}

#' Create groups based on given row indexes or identifiers
#' 
#' @param session Shiny session
#' @param rows Character: comma-separated row indexes or identifiers
#' @param dataset Matrix or data frame: dataset
#' @param identifiers Character: available identifiers
#' 
#' @importFrom shiny tags
#' @return Character: values based on given row indexes or identifiers
createGroupById <- function(session, rows, dataset, identifiers) {
    # Check which strings match available identifiers
    matched <- rows %in% identifiers
    match <- match(rows[matched], identifiers)
    
    # Convert remaining parsable strings to a sequence of numbers
    parsable <- grepl("^[0-9]*:?[0-9]*$", rows[!matched])
    parsed <- unlist(lapply(rows[!matched][parsable],
                            function(row) eval(parse(text=row))))
    parsed <- sort(unique(parsed))
    # Check indexes higher than the number of patients available
    valid <- parsed <= length(identifiers)
    
    # Warn about invalid input
    invalid <- union(rows[!matched][!parsable], parsed[!valid])
    if (length(invalid) > 0) {
        discarded <- paste(invalid, collapse=", ")
        warningAlert(
            session, "The following ", length(invalid),
            " indexes or identifiers were discarded:", tags$code(discarded))
    }
    rows <- rownames(dataset)[unique(union(match, parsed[valid]))]
    return(rows)
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

#' Perform set operations on selected groups
#'
#' @param operation Character: set operation
#' @param groups Matrix: groups
#' @param selected Integer: index of rows regarding selected groups
#' @param symbol Character: Unicode symbol to visually indicate the operation
#' performed (" " by default)
#' @param groupName Character: group name (automatically created if NULL or
#' \code{""})
#' @param patients Character: all patients (required when performing the
#' \code{complement} operation)
#' @param samples Character: all samples (required when performing the 
#' \code{complement} operation)
#' @param matches Character: match between samples (as names) and patients (as
#' values)
#' 
#' @return Matrix containing groups (new group is in the first row)
setOperation <- function(operation, groups, selected, symbol=" ", 
                         groupName=NULL, patients=NULL, samples=NULL,
                         matches=NULL) {
    # Create new set
    new <- NULL
    if (!identical(operation, "remove") && !identical(operation, "rename")) {
        mergedFields <- lapply(1:3, function(i) {
            start <- NULL
            if (operation == "complement" && symbol != " ") {
                start  <- symbol
                symbol <- " \u222A "
            }
            
            if (length(selected) > 1) {
                names <- paste(groups[selected, i], collapse=symbol)
                # Add parenthesis around new expression
                names <- paste0("(", names, ")")
            } else {
                names <- as.character(groups[selected, i])
            }
            
            if (!is.null(start)) {
                names <- paste0("(", start, names, ")")
            }
            return(names)
        })
        if ( !is.null(groupName) && !identical(groupName, "") )
            mergedFields[[1]] <- groupName
        ncol <- 3
        
        operate <- function(col) {
            if (operation == "union") {
                Reduce("union", groups[selected, col])
            } else if (operation == "intersect") {
                Reduce("intersect", groups[selected, col])
            } else if (operation == "subtract") {
                if (length(selected) != 2)
                    stop("set subtract requires 2 groups")
                setdiff(groups[[selected[1], col]], groups[[selected[2], col]])
            } else if (operation == "symDiff") {
                setdiff(Reduce("union", groups[selected, col]),
                        Reduce("intersect", groups[selected, col]))
            } else if (operation == "complement") {
                if (col == "Samples")
                    universe <- samples
                else
                    universe <- patients
                setdiff(universe, Reduce(union, groups[selected, col]))
            }
        }
        
        if ("Samples" %in% colnames(groups)) {
            samples <- operate("Samples")
            ncol <- ncol + 1
        } else {
            samples <- NULL
        }
        
        if ("Patients" %in% colnames(groups)) {
            patients <- operate("Patients")
            if (!is.null(samples) && !is.null(matches)) {
                # Include patients if their samples were included in same group
                matched  <- unname(matches[samples])
                matched  <- matched[!is.na(matched)]
                if (length(matched) > 0)
                    patients <- unique(c(patients, matched))
            }
            ncol <- ncol + 1
        } else {
            patients <- NULL
        }
        
        if (!is.null(patients)) patients <- list(patients)
        if (!is.null(samples)) samples <- list(samples)
        new <- matrix(c(mergedFields, patients, samples), ncol=ncol)
    }
    
    # Add new groups to the top
    if (!is.null(new)) {
        new <- renameGroups(new, groups)
        groups <- rbind(new, groups)
    }
    
    if (identical(operation, "remove")) { # || input$removeSetsUsed)
        # Remove selected groups
        groups <- groups[-selected, , drop=FALSE]
    } else if (identical(operation, "rename")) {
        # Rename selected group
        renamed <- groups[selected, , drop=FALSE]
        renamed[, "Names"] <- groupName
        
        # Avoid groups with the same name
        renamed <- renameGroups(renamed, groups[-selected, , drop=FALSE])
        groups[selected, ] <- renamed[ , ]
        rownames(groups)[selected] <- groups[selected, "Names"]
    }
    return(groups)
}

#' Set operations on groups
#' 
#' This function can be used on groups to merge, intersect, subtract, etc.
#' 
#' @param input Shiny input
#' @param session Shiny session
#' @param buttonId Character: ID of the button to trigger operation
#' @param datasetName Character: name of dataset
#' @param sharedData Shiny app's global variable
#' 
#' @inheritParams setOperation
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
operateOnGroups <- function(input, session, operation, buttonId, symbol=" ",
                            datasetName, sharedData=sharedData) {
    ns <- session$ns
    # Operate on selected groups when pressing the corresponding button
    observeEvent(input[[paste(buttonId, "button", sep="-")]], {
        # Get groups from the dataset
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        selected <- input$groupsTable_rows_selected
        
        if (operation == "subtract") {
            selected  <- sort(selected)
        } else if (operation == "subtract2") {
            selected  <- rev(sort(selected))
            operation <- "subtract"
        }
        
        groupName <- input$groupName
        patients  <- seq(getPatientId())
        samples   <- seq(getSampleId())
        matches   <- getClinicalMatchFrom("Inclusion levels")
        groups    <- setOperation(operation, groups, selected, symbol, 
                                  groupName, patients, samples, matches)
        setGroupsFrom(datasetName, groups)
    })
}

#' Present groups table
#' 
#' @param datasetName Character: name of dataset
#' 
#' @return Matrix with groups ordered (or NULL if no groups exist)
showGroupsTable <- function(datasetName) {
    groups <- getGroupsFrom(datasetName, complete=TRUE)
    
    # Show groups only if there is at least one group
    if (!is.null(groups) && nrow(groups) > 0) {
        show <- NULL
        
        # Show number of patients for each group (if available)
        showPatients <- "Patients" %in% colnames(groups)
        if (showPatients) {
            patients <- lapply(groups[ , "Patients"], length)
            groups[ , "Patients"] <- unlist(patients)
            show <- 4
        }
        
        # Show number of samples for each group (if available)
        showSamples <- "Samples" %in% colnames(groups)
        if (showSamples) {
            samples <- lapply(groups[ , "Samples"], length)
            groups[ , "Samples"] <- unlist(samples)
            if (showPatients)
                show <- c(4, 5)
            else
                show <- 4
        }
        
        # Ordering the groups (plus safety net for cases with one row)
        ord <- c(1, show, 2, 3)
        ordered <- groups[ , ord, drop=FALSE]
        return(ordered)
    } else {
        return(NULL)
    }
}

#' Server function for data grouping
#'
#' @inheritParams operateOnGroups
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom shinyjs disabled enable disable hidden show hide
#' @importFrom shiny textInput
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
groupsServer <- function(input, output, session, datasetName) {
    ns <- session$ns
    
    # Create new group(s)
    createGroupOptions <- function(id) {
        observeEvent(input[[paste0("createGroupAttribute", id)]], {
            createGroup(session, input, output, id, type="Attribute")
        })
        
        observeEvent(input[[paste0("createGroupRows", id)]], {
            createGroup(session, input, output, id, type="Index/Identifier")
        })
        
        observeEvent(input[[paste0("createGroupSubset", id)]], {
            createGroup(session, input, output, id, type="Subset")
        })
        # Update available attributes to suggest in the subset expression
        output[[paste0("groupExpressionSuggestions", id)]] <- renderUI({
            if (id == "Patients")
                dataset <- getClinicalData()
            else if (id == "Samples")
                dataset <- getSampleInfo()
            textSuggestions(ns(paste0("groupExpression", id)), names(dataset))
        })
        
        observeEvent(input[[paste0("createGroupRegex", id)]], {
            createGroup(session, input, output, id, type="Regex")
        })
    }
    
    observe( createGroupOptions("Patients") )
    observe( createGroupOptions("Samples") )
    
    # Render groups list and show interface to manage groups
    output$groupsTable <- renderDataTable(
        showGroupsTable(datasetName), style="bootstrap", escape=FALSE, 
        server=TRUE, rownames=FALSE,
        options=list(pageLength=10, lengthChange=FALSE, scrollX=TRUE,
                     ordering=FALSE,
                     # Stack DataTable elements so they fit in the container
                     dom=paste0(
                         '<"row view-filter"<"col-sm-12"<"pull-left"l>',
                         '<"pull-right"f><"clearfix">>>',
                         'rt<"row view-pager"<"col-sm-12"<"text-center"ip>>>')))
    
    # Remove selected groups
    removeId <- "removeGroups"
    operateOnGroups(input, session, operation="remove", buttonId=removeId,
                    datasetName=datasetName)
    
    # Rename selected groups
    renameId <- "renameGroupName"
    operateOnGroups(input, session, operation="rename", buttonId=renameId,
                    datasetName=datasetName)
    
    # Merge selected groups
    mergeId <- "mergeGroups"
    operateOnGroups(input, session, operation="union", buttonId=mergeId, 
                    symbol=" \u222A ", datasetName=datasetName)
    
    # Intersect selected groups
    intersectId <- "intersectGroups"
    operateOnGroups(input, session, operation="intersect",
                    datasetName=datasetName, buttonId=intersectId, 
                    symbol=" \u2229 ")
    
    # The following set operations are organised inside a "More" button
    moreId <- "moreSetOperations"
    
    # Complement of selected group(s)
    complementId <- "complementGroups"
    operateOnGroups(input, session, operation="complement",
                    datasetName=datasetName, buttonId=complementId, 
                    symbol="U \u005c ")
    
    # Subtract elements from selected group
    subtractId  <- "subtractGroups"
    operateOnGroups(input, session, operation="subtract",
                    datasetName=datasetName, buttonId=subtractId, 
                    symbol=" \u005c ")
    
    subtract2Id <- "subtract2Groups"
    operateOnGroups(input, session, operation="subtract2",
                    datasetName=datasetName, buttonId=subtract2Id, 
                    symbol=" \u005c ")
    
    # Symmetric difference of selected group(s)
    symDiffId <- "symDiffGroups"
    operateOnGroups(input, session, operation="symDiff",
                    datasetName=datasetName, buttonId=symDiffId,
                    symbol=" \u2A01 ")
    
    # Disable set operations according to the selected rows
    observe({
        getButtonId <- function(id) paste(id, "button", sep="-")
        mergeButton      <- getButtonId(mergeId)
        intersectButton  <- getButtonId(intersectId)
        removeButton     <- getButtonId(removeId)
        
        getListId <- function(id) paste(id, "list", sep="-")
        moreButton       <- moreId
        complementButton <- getListId(complementId)
        subtractButton   <- getListId(subtractId)
        subtractButton2  <- getListId(subtract2Id)
        symDiffButton    <- getListId(symDiffId)
        
        selectedRows <- length(input$groupsTable_rows_selected)
        # Number of groups selected to enable each set operation
        # - complement: >= 1
        # - remove:     >= 1
        # - subtract:   2
        # - merge:      >= 2
        # - intersect:  >= 2
        # - sym diff:   >= 2
        
        if (selectedRows >= 1) {
            enable(moreButton)
            enable(complementButton)
            enable(removeButton)
        } else {
            disable(moreButton)
            disable(complementButton)
            disable(removeButton)
        }
        
        if (selectedRows == 2) {
            enable(subtractButton)
            enable(subtractButton2)
        } else {
            disable(subtractButton)
            disable(subtractButton2)
        }
        
        if (selectedRows >= 2) {
            enable(mergeButton)
            enable(intersectButton)
            enable(symDiffButton)
        } else {
            disable(mergeButton)
            disable(intersectButton)
            disable(symDiffButton)
        }
    })
    
    # Show group rename if only one group is selected
    observe({
        if (length(input$groupsTable_rows_selected) == 1)
            show("renameAlert", anim=TRUE, time=0.2)
        else
            hide("renameAlert", anim=TRUE, time=0.2)
    })
    
    # Disable rename button if no new name was given
    observe({
        renameButton <- paste(renameId, "button", sep="-")
        if (is.null(input$groupName) || input$groupName == "")
            disable(renameButton)
        else
            enable(renameButton)
    })
    
    # Render groups interface only if at least one group exists
    output$groupsList <- renderUI({
        groups <- getGroupsFrom(datasetName, complete=TRUE)
        
        operationElement <- function(operation, ..., id=NULL, icon=NULL,
                                     FUN=actionButton, width=NULL) {
            buttonId <- paste(id, "button", sep="-")
            button <- FUN(buttonId, operation, icon=icon, width=width, ...)
            
            if (identical(FUN, actionLink)) {
                itemId <-paste(id, "list", sep="-")
                button <- tags$li(id=itemId, button)
            }
            return(button)
        }
        
        operationButton <- function(...) operationElement(..., FUN=actionButton)
        operationLink   <- function(...) operationElement(..., FUN=actionLink)
        
        # Don't show anything when there are no groups
        if (!is.null(groups) && nrow(groups) > 0) {
            operations <- div(
                id=ns("setOperations"), class="btn-group",
                operationButton("Merge", id=ns(mergeId),
                                icon=setOperationIcon("union")),
                operationButton("Intersect", id=ns(intersectId),
                                icon=setOperationIcon("intersect")),
                tags$div(
                    class="btn-group", role="group",
                    tags$button("More", id=ns(moreId), tags$span(class="caret"),
                                class="btn btn-default dropdown-toggle",
                                "data-toggle"="dropdown",
                                "aria-haspopup"="true", "aria-expanded"="true"),
                    tags$ul(
                        class="dropdown-menu",
                        operationLink(
                            "Complement", id=ns(complementId),
                            helpText("Create a group with elements outside the",
                                     "selected group(s)"),
                            icon=setOperationIcon("complement-AB")),
                        operationLink(
                            "Subtract elements from upper-selected group",
                            helpText("Create a group with the exclusive",
                                     "elementes from the upper-selected group"),
                            id=ns(subtractId),
                            icon=setOperationIcon("difference-AB")),
                        operationLink(
                            "Subtract elements from lower-selected group",
                            helpText("Create a group with the exclusive",
                                     "elementes from the lower-selected group"),
                            id=ns(subtract2Id),
                            icon=setOperationIcon("difference-BA")),
                        operationLink(
                            "Symmetric difference",
                            helpText("Create a group with the non-intersecting",
                                     "elements of selected groups"),
                            id=ns(symDiffId),
                            icon=setOperationIcon("symmetric-difference")))),
                operationButton("Remove", id=ns(removeId), class="btn-danger",
                                icon=icon("times")))
            
            removeAllButton <- actionButton(ns("removeAll"), class="btn-danger",
                                            class="pull-right",
                                            "Remove all groups",
                                            icon=icon("trash"))
            
            renameButton <- operationButton("Rename", id=ns(renameId),
                                            class="pull-right",
                                            icon=icon("pencil"))
            nameField <- textInput(ns("groupName"), label=NULL, 
                                   placeholder="Rename selected group")
            nameField$attribs$style <- "margin: 0"
            renameInterface <- div(
                id=ns("renameAlert"), class="alert", role="alert",
                class="alert-info",
                style=" margin-top: 10px; margin-bottom: 0px;",
                div(class="input-group", nameField,
                    div(class="input-group-btn", renameButton)))
            
            tagList(
                hr(),
                dataTableOutput(ns("groupsTable")),
                helpText("Select groups by clicking on each one to perform the",
                         "following actions on them."),
                disabled(operations),
                removeAllButton,
                #checkboxInput(ns("removeSetsUsed"), "Remove original groups",
                #              value=TRUE)
                hidden(renameInterface)
            )
        }
    })
    
    observeEvent(input$removeAll, setGroupsFrom(datasetName, NULL))
}

#' Server function for data grouping (one call)
#' 
#' These functions only run once instead of running for every instance of groups
#' 
#' @inheritParams groupsServer
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
groupsServerOnce <- function(input, output, session) {
    # Update groups according to the availability of sample identifiers
    observe({
        group <- getGroupsFrom("Clinical data", complete=TRUE)
        # Ignore if there are no groups
        if (is.null(group)) return(NULL)
        
        samples <- getSampleId()
        clinical <- getClinicalData()
        match <- getClinicalMatchFrom("Inclusion levels")
        showSamples <- "Samples" %in% colnames(group)
        
        if ( is.null(samples) && showSamples ) {
            # Remove sample identifiers from groups
            if ( "Patients" %in% colnames(group) ) {
                # Groups are made with patients and samples; only remove samples
                group <- group[ , -match("Samples", colnames(group))]
            } else {
                # Groups are only made of samples; remove all groups
                group <- NULL
            }
            setGroupsFrom("Clinical data", group)
        } else if ( !is.null(samples) && !showSamples && !is.null(clinical) &&
                    !is.null(match) ) {
            # Update groups if previously made with patients only
            patients <- group[ , "Patients"]
            samples <- getMatchingSamples(patients, samples, clinical,
                                          match=match)
            group <- cbind(group, "Samples"=samples)
            setGroupsFrom("Clinical data", group)
        }
    })
}

#' Get selected groups for a given group selection element
#' 
#' @param input Shiny input
#' @param id Character: identifier of the group selection element
#' @inheritParams getGroupsFrom
#' @param filter Character: only get groups passed
#' 
#' @return List with selected groups (or NULL if no groups were selected)
getSelectedGroups <- function(input, id, samples=FALSE, dataset="Clinical data",
                              filter=NULL) {
    selection <- input[[paste0(id, "Selection")]]
    noGroups  <- !is.null(selection) && selection == "noGroups"
    selected  <- input[[id]]
    
    if ( noGroups || is.null(selected) ) {
        # User selects no groups (either explicitly or not)
        groups <- NULL
    } else {
        groups <- getGroupsFrom(dataset, samples=samples)[selected]
        
        if (!is.null(filter))
            groups <- lapply(groups, function(i) i[i %in% filter])
    }
    return(groups)
}

# attr(groupsUI, "loader") <- "data"
# attr(groupsServer, "loader") <- "data"