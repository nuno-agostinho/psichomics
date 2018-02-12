#' Group selection
#' 
#' Group selection interface and logic
#' 
#' @param id Character: identifier
#' @param label Character: selectize label
#' @param placeholder Character: selectize placeholder
#' @param noGroupsLabel Character: label to show when no groups may be selected
#' (if NULL, the option to show no groups will not be shown)
#' @param groupsLabel Character: label to show to the option of using groups
#' when no groups may be selected
#' @param maxItems Numeric: maximum number of selected items
#' 
#' @importFrom shiny fluidRow column uiOutput selectizeInput actionButton
#' radioButtons actionLink downloadLink
#' 
#' @note To allow the user to (explicitly) select no groups, pass the 
#' \code{noGroupsLabel} and \code{groupsLabel} arguments.
#' 
#' @return \code{selectGroupsUI}: Interface for group selection
selectGroupsUI <- function (
    id, label, placeholder="Click on 'Groups' to create or edit groups",
    noGroupsLabel=NULL, groupsLabel=NULL, maxItems=NULL) {
    editId <- paste0(id, "Edit")
    groupSelect <- selectizeInput(
        id, label, choices=NULL, multiple=TRUE, width="auto", options=list(
            plugins=list('remove_button', 'drag_drop'), maxItems=maxItems, 
            searchField=list("value", "label"),
            placeholder=placeholder, render = I(
                '{ option: function(item, escape) {
                       return "<div><b>" + escape(item.value) + "</b><small> " + 
                           escape(item.label) + "</small></div>"; },
                   item: function(item, escape) {
                       return "<div><b>" + escape(item.value) + "</b><small> " + 
                           escape(item.label) + "</small></div>";
                   } }')))
    
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
        label, column(10, groupSelect), column(2, actionButton(
            editId, "Groups", class="pull-right", class="btn-info",
            style="z-index: 1; position: relative;")))
    
    if ( !is.null(noGroupsLabel) ) {
        noGroupsId <- paste0(id, "Selection")
        
        choices <- c("noGroups", "groups")
        names(choices) <- c(noGroupsLabel, groupsLabel)
        
        select <- tagList(
            radioButtons(noGroupsId, radioLabel, choices=choices, width="100%"),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", noGroupsId, "groups"), 
                select))
    }
    return(select)
}

#' @rdname selectGroupsUI
#' 
#' @inheritParams getGroups
#' @param session Shiny session
#' @param preference Character: name of groups to pre-select, when available
#' (if NULL, all groups will be pre-selected)
#' 
#' @importFrom shinyjs enable disable onclick toggleClass runjs
#' 
#' @return \code{selectGroupsServer}: Server logic for group selection
selectGroupsServer <- function(session, id, type, preference=NULL) {
    ns     <- session$ns
    input  <- session$input
    output <- session$output
    
    editId  <- paste0(id, "Edit")
    modalId <- paste0(id, "Modal")
    
    # Open correct group interface
    observeEvent(input[[editId]], runjs(sprintf("showGroups('%s')", type)))
    
    # Update groups shown in the selectize element
    observe({
        groupTable <- getGroups(type, complete=TRUE)
        groups <- rownames(groupTable)
        if (is.null(groups)) {
            # Disable selection and animate button when clicking disabled input
            groups <- list()
            disable(id)
            onclick(id, runjs(paste0("$('#", ns(editId), 
                                     "').animateCss('rubberBand');")))
        } else {
            enable(id)
            onclick(id, NULL)
            
            if (type %in% c("Patients", "Samples")) {
                elem1 <- "Patients"
                elem2 <- "Samples"
            } else if (type %in% c("ASevents", "Genes")) {
                elem1 <- "ASevents"
                elem2 <- "Genes"
            }
            
            elem1Number <- NULL
            if (elem1 %in% colnames(groupTable))
                elem1Number <- sapply(groupTable[ , elem1], length)
            
            elem2Number <- NULL
            if (elem2 %in% colnames(groupTable))
                elem2Number <- sapply(groupTable[ , elem2], length)
            
            elem1 <- gsub("ASevents", "AS events", elem1)
            if (!is.null(elem1Number) && !is.null(elem2Number))
                ns <- sprintf("%s %s, %s %s", elem1Number, elem1, elem2Number,
                              elem2)
            else if (!is.null(elem1Number))
                ns <- sprintf("%s %s", elem1Number, elem1)
            else if (!is.null(elem2Number))
                ns <- sprintf("%s %s", elem2Number, elem2)
            names(groups) <- ns
        }
        
        currentSelection <- isolate(input[[id]])
        if (is.null(currentSelection)) {
            if (is.null(preference))
                selected <- groups
            else
                selected <- groups[groups %in% preference]
        } else {
            selected <- currentSelection[currentSelection %in% groups]
            if (length(selected) == 0) selected <- groups
        }
        updateSelectizeInput(session, id, choices=groups, selected=selected)
    })
}

#' Interface to manipulate data grouping
#' 
#' @param id Character: identifier
#' @param type Character: type of data for each the interface is intended
#' 
#' @return HTML elements
groupManipulationInput <- function(id, type) {
    ns <- NS(id)
    
    multiFisherTests <- TRUE
    if (type == "ASevents") multiFisherTests <- FALSE
    
    modalId <- gsub("(.*)Call\\-.*", "\\1Show", ns(id))
    missingData <- function(dataset, message, linkText) {
        div(class="alert alert-danger", role="alert",
            icon("exclamation-circle"), message,
            tags$a(linkText, onclick=loadRequiredData( modalId )))
    }
    
    groupOptions <- function(id, title) {
        if (id == "Patients") {
            example <- tagList(
                "For instance, to create groups by tumour stage, type",
                tags$b("tumor_stage"), "and select the first suggestion.")
        } else if (id == "Samples") {
            example <- NULL
        }
        
        cols <- c("No attributes available"="")
        indexIdentifiersUI <- groupById(ns, id)
        if (id %in% c("Patients", "Samples")) {
            navbarMenu(
                title,
                tabPanel("Attribute", groupByAttribute(ns, cols, id, example)),
                tabPanel("Index/Identifier", indexIdentifiersUI),
                "----", "Advanced options",
                tabPanel("Subset expression", groupByExpression(ns, id)),
                tabPanel("Regular expression", groupByGrep(ns, cols, id)))
        } else {
            tabPanel(title, indexIdentifiersUI)
        }
    }
    
    if (type == "Samples") {
        title  <- "By patients"
        first  <- groupOptions("Patients", title)
        # firstAlert  <- tabPanel(title, value="NoPatients", missingData(
        #     "Clinical data", "No clinical data loaded to group by patients.",
        #     "Please, load clinical data."))
        
        title  <- "By samples"
        second <- groupOptions("Samples", title)
        # secondAlert <- tabPanel(title, value="NoSamples", missingData(
        #     "Inclusion levels",
        #     "No sample information loaded to group by samples.",
        #     "Please, load sample information."))
    } else if (type == "ASevents") {
        title  <- "By splicing events"
        first  <- groupOptions("ASevents", title)
        # firstAlert <- tabPanel(title, value="ASevents", missingData(
        #     "Alternative splicing events",
        #     "No alternative splicing events are available.",
        #     "Please, load or quantify splicing events."))
        
        title  <- "By genes"
        second <- groupOptions("Genes", title)
        # secondAlert <- tabPanel(title, value="NoGenes", missingData(
        #     "Genes", "No genes or alternative splicing events are available.",
        #     paste("Please, load genes and/or load and quantify",
        #           "splicing events.")))
    }
    
    sidebarLayout(
        sidebarPanel( 
            tabsetPanel(id=ns("groupBy"), type="pills", first, second) ),
        mainPanel( renderGroupInterface(ns, multiFisherTests) ))
}

#' @rdname appUI
#' @importFrom shinyBS bsCollapse bsCollapsePanel
groupsUI <- function(id, tab) {
    ns <- NS(id)
    
    tab(icon="object-group", title="Groups",
        uiOutput(ns("alert")),
        tabsetPanel(
            id="groupsTypeTab",
            tabPanel("Patient and sample groups",
                     groupManipulationInput(ns("sampleGroupModule"), 
                                            "Samples")),
            tabPanel("Splicing event and gene groups",
                     groupManipulationInput(ns("ASeventGroupModule"),
                                            "ASevents"))))
}

#' Render group interface
#' 
#' @param ns Namespace function
#' @param multiFisherTests Boolean: allow to perform multiple Fisher exact test
#' between groups
#' 
#' @importFrom shiny actionButton actionLink downloadLink tags
#' @importFrom shinyjs disabled
#' @importFrom colourpicker colourInput
#' 
#' @return HTML elements
renderGroupInterface <- function(ns, multiFisherTests=TRUE) {
    renameId    <- "renameGroupName"
    setColourId <- "setGroupColour"
    removeId    <- "removeGroups"
    
    # Set operations
    mergeId      <- "mergeGroups"
    intersectId  <- "intersectGroups"
    moreId       <- "moreSetOperations"
    complementId <- "complementGroups"
    subtractId   <- "subtractGroups"
    subtract2Id  <- "subtract2Groups"
    symDiffId    <- "symDiffGroups"
    
    # Save and load groups
    saveLoadId           <- "saveLoadGroups"
    groupIndependenceId  <- "groupTesting"
    saveSelectedGroupsId <- "saveSelectedGroups"
    saveAllGroupsId      <- "saveAllGroups"
    loadGroupsId         <- "loadGroups"
    testGroupsId         <- "testGroups"
    
    operationElement <- function(operation, ..., id=NULL, icon=NULL,
                                 FUN=actionButton, width=NULL, disable=TRUE) {
        buttonId <- paste(id, "button", sep="-")
        button <- FUN(buttonId, operation, icon=icon, width=width, ...)
        
        if (identical(FUN, actionLink) || identical(FUN, downloadLink)) {
            itemId <- paste(id, "list", sep="-")
            button <- tags$li(id=itemId, button)
        }
        
        if (disable) button <- disabled(button)
        return(button)
    }
    
    operationButton <- function(...) operationElement(..., FUN=actionButton)
    operationLink   <- function(...) operationElement(..., FUN=actionLink)
    downloadContent <- function(...) operationElement(..., FUN=downloadLink)
    
    # Set operations
    complementLink <- operationLink(
        "Complement", id=ns(complementId),
        helpText("Create a group with the elements outside the",
                 "selected group(s)"),
        icon=setOperationIcon("complement-AB"),
        disable=FALSE)
    
    subtractLink <- operationLink(
        "Subtract elements from upper-selected group",
        helpText("Create a group with the exclusive",
                 "elements from the upper-selected group"),
        id=ns(subtractId),
        icon=setOperationIcon("difference-AB"))
    
    subtract2Link <- operationLink(
        "Subtract elements from lower-selected group",
        helpText("Create a group with the exclusive",
                 "elements from the lower-selected group"),
        id=ns(subtract2Id),
        icon=setOperationIcon("difference-BA"))
    
    symDiffLink <- operationLink(
        "Symmetric difference",
        helpText("Create a group with the non-intersecting",
                 "elements of selected groups"),
        id=ns(symDiffId),
        icon=setOperationIcon("symmetric-difference"))
    
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
            tags$ul(class="dropdown-menu dropdown-menu-right",
                    complementLink, subtractLink, subtract2Link, symDiffLink)))
    
    # Save and load groups
    saveSelectedGroupsLink <- downloadContent(
        icon("user"), class=NULL, 
        "Save elements from selected group(s)",
        helpText("Export a file containing identifiers by select groups"),
        id=ns(saveSelectedGroupsId))
    
    saveAllGroupsLink <- downloadContent(
        icon("users"), class=NULL, "Save elements from all groups",
        helpText("Export a file containing identifiers by every group"),
        id=ns(saveAllGroupsId),
        disable=FALSE)
    
    loadGroupsLink <- operationLink(
        "Load groups", 
        helpText("Import a file containing identifiers by group"),
        id=ns(loadGroupsId), icon=icon("plus-circle"), disable=FALSE)
    
    saveLoadGroups <- tags$div(
        class="btn-group", role="group",
        tags$button(icon("folder-open"), id=ns(saveLoadId),
                    tags$span(class="caret"),
                    class="btn btn-default dropdown-toggle",
                    "data-toggle"="dropdown",
                    "aria-haspopup"="true", "aria-expanded"="true"),
        tags$ul(class="dropdown-menu dropdown-menu-right",
                saveSelectedGroupsLink, saveAllGroupsLink,
                tags$li(role="separator", class="divider"), loadGroupsLink))
    
    if (multiFisherTests) {
        testGroupIndependenceLink <- operationLink(
            "Test group independence against subject and sample categories",
            helpText(
                "Perform multiple independent Fisher's exact tests between",
                "each of the selected groups and their complementary group",
                "against the values from categorical variables"),
            id=ns(testGroupsId), icon=icon("plus-circle"), 
            style="white-space: normal;")
        
        groupIndependenceSet <- tags$div(
            class="btn-group", role="group",
            tags$button(icon("link"), id=ns(groupIndependenceId),
                        tags$span(class="caret"),
                        class="btn btn-default dropdown-toggle",
                        "data-toggle"="dropdown",
                        "aria-haspopup"="true", "aria-expanded"="true"),
            tags$ul(class="dropdown-menu dropdown-menu-right", 
                    style="width: 300px;", testGroupIndependenceLink))
    }
    
    removeButton <- operationButton("Remove", id=ns(removeId),
                                    class="btn-danger",
                                    icon=icon("times"))
    
    removeAllLink <- operationLink("Remove all groups",
                                   id=ns("removeAll"),
                                   class="li-danger",
                                   icon=icon("trash"),
                                   disable=FALSE)
    
    # Rename interface
    renameButton <- operationButton("Rename", id=ns(renameId),
                                    class="pull-right", icon=icon("pencil"))
    nameField <- textInput(ns("groupName"), label=NULL,
                           placeholder="Rename selected group")
    nameField$attribs$style <- "margin: 0; width: auto;"
    
    # Colour selection interface
    colourSelector <- colourInput(ns("groupColour"), label=NULL)
    colourSelector[[2]][["class"]] <- paste(colourSelector[[2]][["class"]],
                                            "groups-colourpicker")
    colourSelector[[2]][["style"]] <- paste("margin-bottom: 0px !important;",
                                            "width: auto;")
    setColourButton <- operationButton("Set colour", id=ns(setColourId),
                                       class="pull-right", disable=FALSE,
                                       icon=icon("paint-brush"))
    
    singleGroupSelectedInterface <- div(
        id=ns("singleGroupSelected"), class="alert", role="alert",
        class="alert-info", style="margin-top: 10px; margin-bottom: 0px;",
        fluidRow(
            column(7, div(class="input-group", nameField,
                          div(class="input-group-btn", renameButton))),
            column(5, div(class="input-group", colourSelector,
                          div(class="input-group-btn", setColourButton)))))
    
    tagList(
        dataTableOutput(ns("groupsTable")),
        helpText("Select groups by clicking on each one to perform the",
                 "following actions on them."),
        operations, saveLoadGroups, 
        if (multiFisherTests) groupIndependenceSet,
        tags$div(class="btn-group pull-right",
                 removeButton,
                 tags$button(type="button", tags$span(class="caret"),
                             class="btn btn-danger dropdown-toggle",
                             "data-toggle"="dropdown",
                             "aria-haspopup"="true",
                             "aria-expanded"="false"),
                 tags$ul(class="dropdown-menu",
                         style="background-color: #d9534f;",
                         style="border-color: #d43f3a;", removeAllLink)),
        # checkboxInput(ns("removeSetsUsed"), "Remove original groups",
        #               value=TRUE)
        hidden(singleGroupSelectedInterface),
        plotOutput(ns(groupTestId <- "groupIndependenceTestingPlot"), 
                   height="200px", 
                   hover=hoverOpts(ns(paste0(groupTestId, "-hover")), 
                                   delay=50, delayType="throttle")),
        uiOutput(ns(paste0(groupTestId, "-tooltip"))))
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
                       width="auto", choices=cols,
                       options=list(lockOptgroupOrder=TRUE)),
        actionButton(ns(paste0("createGroupAttribute", id)), "Create group",
                     class ="btn-primary")
    )
}

#' User interface to group by row
#' 
#' @inheritParams groupByAttribute
#' 
#' @importFrom shiny textInput
#' 
#' @return HTML elements
groupById <- function(ns, id) {
    sid <- gsub("s$", "", id)
    sid <- gsub("ASevent", "Splicing event", sid)
    
    tagList(
        selectizeInput(
            ns(paste0("groupRows", id)), paste(sid, "indexes or identifiers"),
            choices=NULL, multiple=TRUE, width="auto", options=list(
                create=TRUE, createOnBlur=TRUE, # Allow to add new items
                plugins=list('remove_button'), persist=FALSE)),
        helpText("Example: ", tags$kbd("1:6, 8, 10:19"), "creates a group with",
                 "rows 1 to 6, 8 and 10 to 19. You can also input identifiers",
                 "instead of indexes."),
        textInput(ns(paste0("groupNameRows", id)), "Group label", width="auto",
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
        helpText('Example: ', tags$kbd('X > 8 & Y == "alive"'), ' selects rows',
                 'with values higher than 8 for column X and "alive" for',
                 'column Y.'),
        uiOutput(ns(paste0("groupExpressionSuggestions", id))),
        textInput(ns(paste0("groupNameSubset", id)), "Group label", 
                  width="auto", placeholder="Unnamed"),
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
        textInput(ns(paste0("groupNameRegex", id)), "Group label", width="auto",
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
    else
        dataset <- NULL
    
    new <- createGroupFromInput(session, input, output, dataset, id, type)
    if (!is.null(new)) appendNewGroups(id, new)
    
    updateSelectizeInput(session, paste0("groupAttribute", id),
                         selected=character())
}

#' Assign colours to groups
#' 
#' @param new Matrix: groups to which colours will be assigned
#' @param groups Matrix: groups to check which colours are already assigned
#' 
#' @return Groups with an added column to state the colour
assignColours <- function(new, groups=NULL) {
    strong <- c("#08419E", "#EF9636", "#D33E6A", "#00C652", 
                "#4C71DB", "#8F033B", "#F89CD1", "#05CFC0")
    medium <- c("#7D87B6", "#EFB893", "#E17A90", "#8CD59A", 
                "#8696DC", "#9E646F", "#F6C4DF", "#9CDDD5")
    light  <- c("#BEC1D2", "#EAD2C7", "#E7AFBA", "#C6DECA", 
                "#B6BBE0", "#D6BBC0", "#F2E1EA", "#D3E7E5")
    colours <- c(strong, medium, light)
    
    # Avoid setting colours previously assigned
    priority <- NULL
    if (!is.null(groups) && "Colour" %in% colnames(groups))
        priority <- colours[!colours %in% groups[ , "Colour"]]
    
    # Repeat default colours when there are many groups
    reps    <- ceiling((nrow(new) - length(priority)) / length(colours))
    colours <- rep(colours, reps)
    colours <- c(priority, colours)[nrow(new):1]
    
    hasColnames <- !is.null(colnames(new))
    new <- cbind(new, colours)
    if (hasColnames)
        colnames(new)[ncol(new)] <- "Colour"
    else
        colnames(new) <- NULL
    return(new)
}

#' Append new groups to already existing groups
#' 
#' Retrieve previous groups, rename duplicated group names in the new groups
#' and append new groups to the previous ones
#' 
#' @inheritParams getGroups
#' @param new Rows of groups to be added
#' @param clearOld Boolean: clear old groups?
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
appendNewGroups <- function(type, new, clearOld=FALSE) {
    # Rename duplicated group names
    if (clearOld) 
        groups <- NULL
    else
        groups <- getGroups(type, complete=TRUE)
    
    new <- renameGroups(new, groups)
    new <- assignColours(new, groups)
    
    if (clearOld) {
        groups <- new
    } else {
        # Append the new group(s) to the groups already created
        groups <- rbind(new, groups)
    }
    setGroups(type, groups)
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
                              "Samples"=getSampleId(),
                              "ASevents"=getASevents(),
                              "Genes"=getGenes())
        allRows <- createGroupById(session, rows, identifiers)
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
    
    if (id %in% c("Patients", "Samples")) {
        patients <- getPatientId()
        samples  <- getSampleId()
        match    <- getClinicalMatchFrom("Inclusion levels")
        
        # Match patients with samples (or vice-versa)
        if (!is.null(patients) && !is.null(samples) && !is.null(match)) {
            if (id == "Patients") {
                patients <- group[ , "Patients"]
                samples <- getMatchingSamples(patients, samples, patients,
                                              match=match)
                group <- cbind(group, "Samples"=samples)
            } else if (id == "Samples") {
                samples2patients <- function(i, match) {
                    m <- match[i]
                    return(unique(m[!is.na(m)]))
                }
                patients <- lapply(group[ , "Samples"], samples2patients, match)
                group <- cbind(group, "Patients"=patients)
                group <- group[ , c(1:3, 5, 4), drop=FALSE]
            }
        }
    } else if (id %in% c("ASevents", "Genes")) {
        # Match AS events with genes (or vice-versa)
        ASevents <- getASevents()
        if (!is.null(ASevents)) {
            if (id == "ASevents") {
                ASevents <- group[ , "ASevents"]
                genes <- lapply(group[ , "ASevents"], parseSplicingEvent)
                genes <- lapply(genes, "[[", "gene")
                uniqueUnlist <- function(i) unique(unlist(i))
                group <- cbind(group, "Genes"=lapply(genes, uniqueUnlist))
            } else if (id == "Genes") {
                genes <- group[ , "Genes"]
                ASeventGenes <- matchSplicingEventsWithGenes(ASevents)
                filterBasedOnGenes <- function(gene, ASeventGenes)
                    ASeventGenes[names(ASeventGenes) %in% gene]
                
                # Process TCGA gene ID
                genes <- lapply(genes, function(gene) gsub("\\|.*$", "", gene))
                ASevents <- lapply(genes, filterBasedOnGenes, ASeventGenes)
                group <- cbind(group, ASevents)
                group <- group[ , c(1:3, 5, 4), drop=FALSE]
            }
        } else if (id == "Genes") {
            ASevents <- lapply(seq(nrow(group)), function(i) character(0))
            group <- cbind(group, ASevents)
            group <- group[ , c(1:3, 5, 4), drop=FALSE]
        }
    }
    return(group)
}

#' @rdname createGroupByAttribute
#' @export
createGroupByColumn <- function(col, dataset) {
    .Deprecated("createGroupByAttribute")
    createGroupByAttribute(col, dataset)
}

#' Split elements into groups based on a given column of a dataset
#' 
#' Elements are identified by their respective row name.
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
    if (!col %in% colnames(dataset))
        stop("The given attribute was not found in the column names of the",
             "dataset")
    
    colData <- as.character(dataset[[col]])
    names(colData) <- rownames(dataset)
    
    # Replace missing values for "NA" so they are included by the "which" function
    colData[is.na(colData)] <- "NA"
    
    # Create groups according to the chosen column
    groupNames <- sort(unique(colData))
    group <- lapply(lapply(groupNames, `==`, colData), which)
    names(group) <- groupNames
    group <- lapply(group, names)
    return(group)
}

#' Create groups based on given row indexes or identifiers
#' 
#' @param session Shiny session
#' @param rows Character: comma-separated row indexes or identifiers
#' @param identifiers Character: available identifiers
#' 
#' @importFrom shiny tags
#' @return Character: values based on given row indexes or identifiers
createGroupById <- function(session, rows, identifiers) {
    # Check which strings match available identifiers
    matched <- rows %in% identifiers
    match <- match(rows[matched], identifiers)
    
    # Convert remaining parsable strings to a sequence of numbers
    parsable <- grepl("^[0-9]*:?[0-9]*$", rows[!matched])
    if (length(parsable) > 0) {
        parsed <- unlist(lapply(rows[!matched][parsable],
                                function(row) eval(parse(text=row))))
        parsed <- sort(unique(parsed))
    } else {
        parsed <- NULL
    }
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
    rows <- identifiers[unique(union(match, parsed[valid]))]
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
#' @param first Character: identifiers of the first element (required when 
#' performing the \code{complement} operation)
#' @param second Character: identifiers of the second element (required when 
#' performing the \code{complement} operation)
#' @param matches Character: match between samples (as names) and patients (as
#' values)
#' @param type Character: type of group where set operations are to be performed
#' @param assignColoursToGroups Boolean: assign colours to new groups? FALSE by
#' default
#' 
#' @return Matrix containing groups (new group is in the first row)
setOperation <- function(operation, groups, selected, symbol=" ", 
                         groupName=NULL, first=NULL, second=NULL, matches=NULL,
                         type="Samples",  assignColoursToGroups=FALSE) {
    # Create new set
    new <- NULL
    if (!identical(operation, "remove") && !identical(operation, "rename")) {
        setOperationString <- function(i) {
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
                if (!is.null(selected))
                    names <- paste0("(", start, names, ")")
                else
                    names <- "Universe"
            }
            return(names)
        }
        setOperated <- lapply(1:3, setOperationString)
        if ( !is.null(groupName) && !identical(groupName, "") && 
             operation == "rename")
            setOperated[[1]] <- groupName
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
                if (col == "Samples" || col == "Genes")
                    universe <- second
                else
                    universe <- first
                setdiff(universe, Reduce(union, groups[selected, col]))
            }
        }
        
        if (type == "Samples") {
            firstStr  <- "Patients"
            secondStr <- "Samples"
        } else if (type == "ASevents") {
            firstStr  <- "ASevents"
            secondStr <- "Genes"
        }
        
        if (secondStr %in% colnames(groups) || !is.null(second)) {
            second <- operate(secondStr)
            ncol <- ncol + 1
        } else {
            second <- NULL
        }
        
        if (firstStr %in% colnames(groups) || !is.null(first)) {
            first <- operate(firstStr)
            if (!is.null(second) && !is.null(matches)) {
                # Include the first element if their respective second elements
                # were included in same group
                matched  <- unname(matches[second])
                matched  <- matched[!is.na(matched)]
                if (length(matched) > 0)
                    first <- unique(c(first, matched))
            }
            ncol <- ncol + 1
        } else {
            first <- NULL
        }
        
        if (!is.null(first)) first <- list(first)
        if (!is.null(second)) second <- list(second)
        new <- matrix(c(setOperated, first, second), ncol=ncol)
    }
    
    if (!is.null(new)) {
        # Add new groups to the top
        new <- renameGroups(new, groups)
        if (assignColoursToGroups) new <- assignColours(new, groups)
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
    
    # Correctly fill column names if empty
    if (is.null(colnames(groups))) {
        names <- c("Names", "Subset", "Input")
        if (!is.null(first))              names <- c(names, firstStr)
        if (!is.null(second))             names <- c(names, secondStr)
        if (length(names) < ncol(groups)) names <- c(names, "Colour")
        colnames(groups) <- names
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
#' @param sharedData Shiny app's global variable
#' 
#' @inheritParams setOperation
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
operateOnGroups <- function(input, session, operation, buttonId, symbol=" ",
                            type, sharedData=sharedData) {
    # Operate on selected groups when pressing the corresponding button
    observeEvent(input[[paste(buttonId, "button", sep="-")]], {
        # Get groups from the dataset
        groups <- getGroups(type, complete=TRUE)
        selected <- input$groupsTable_rows_selected
        
        if (operation == "subtract") {
            selected  <- sort(selected)
        } else if (operation == "subtract2") {
            selected  <- rev(sort(selected))
            operation <- "subtract"
        }
        
        groupName <- input$groupName
        if (!is.null(groupName) && groupName != "" && operation == "rename")
            updateTextInput(session, "groupName", value="")
        
        if (type == "Samples") {
            first   <- getPatientId()
            second  <- getSampleId()
            matches <- getClinicalMatchFrom("Inclusion levels")
        } else if (type == "ASevents") {
            first   <- getASevents()
            second  <- getGenes()
            matches <- getGenesFromSplicingEvents(first)
        }
        groups <- setOperation(operation, groups, selected, symbol, groupName, 
                               first, second, matches, type,
                               assignColoursToGroups=TRUE)
        setGroups(type, groups)
    })
}

#' Present groups table
#' 
#' @inheritParams getGroups
#' @importFrom shiny tags
#' 
#' @return Matrix with groups ordered (or NULL if no groups exist)
showGroupsTable <- function(type) {
    groups <- getGroups(type, complete=TRUE)
    
    # Show groups only if there is at least one group
    if (!is.null(groups) && nrow(groups) > 0) {
        show <- NULL
        
        if (type == "Samples") {
            elem1 <- "Patients"
            elem2 <- "Samples"
        } else if (type == "ASevents") {
            elem1 <- "ASevents"
            elem2 <- "Genes"
        }
        
        # Show number of patients or AS events for each group (if available)
        showElem1 <- elem1 %in% colnames(groups)
        if (showElem1) {
            groups[ , elem1] <- unlist( lapply(groups[ , elem1], length) )
            show <- 4
        }
        
        # Show number of samples or genes for each group (if available)
        showElem2 <- elem2 %in% colnames(groups)
        if (showElem2) {
            groups[ , elem2] <- unlist( lapply(groups[ , elem2], length) )
            
            if (showElem1)
                show <- c(4, 5)
            else
                show <- 4
        }
        
        # Show colours
        if ("Colour" %in% colnames(groups)) {
            colour <- match("Colour", colnames(groups))
            groups[ , colour] <- lapply(groups[ , colour], function(col) {
                div <- tags$div(style=paste0(
                    "background-color: ", col, ";",
                    "height: 20px;", "border-radius: 5px;"))
                return(as.character(div))
            })
            show <- c(show, colour)
        }
        
        # Ordering the groups (plus safety net for cases with one row)
        ord <- c(1, show, 2, 3)
        ordered <- groups[ , ord, drop=FALSE]
        colnames(ordered)[1] <- "Group"
        colnames(ordered) <- gsub("ASevents", "AS events", colnames(ordered))
        
        # Unlist items as required for newer versions of DT
        if (nrow(ordered) == 1)
            ordered <- t(apply(ordered, 1, unlist))
        else
            ordered <- apply(ordered, 2, unlist)
        return(data.frame(ordered, stringsAsFactors=FALSE))
    } else {
        return(NULL)
    }
}

#' Logic server to manipulate data grouping
#' 
#' @inheritParams groupsServer
#' @param type Character: type of data for each the interface is intended
#' 
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom shinyjs toggleState enable disable hidden show hide
#' @importFrom shiny textInput
#' @importFrom shinyBS updateCollapse
#' @importFrom colourpicker updateColourInput
#' 
#' @return HTML elements
groupManipulation <- function(input, output, session, type) {
    ns <- session$ns
    
    # Update attributes to select in "group by attributes" panel
    if (type == "Samples") {
        updateAttributes <- function(id) {
            if (id == "Patients") {
                attrs <- getPatientAttributes()
            } else if (id == "Samples") {
                attrs <- getSampleAttributes()
            }
            
            suggestedCols <- attr(attrs, "default")
            if (!is.null(suggestedCols)) {
                suggestedIndex <- match(suggestedCols, attrs)
                suggestedIndex <- suggestedIndex[!is.na(suggestedIndex)]
                cols <- list("Start typing to search for attributes"="",
                             "Suggested attributes"=attrs[suggestedIndex],
                             "Other attributes"=attrs[-suggestedIndex])
            } else {
                cols <- c("Start typing to search for attributes"="", attrs)
            }
            
            updateSelectizeInput(session, paste0("groupAttribute", id),
                                 choices=cols)
        }
        
        observe(updateAttributes("Samples"))
        observe(updateAttributes("Patients"))
    }
    
    # Update identifiers to suggest in index/identifiers panel
    observe({
        if (type == "Samples") {
            updateSelectizeInput(session, paste0("groupRows", "Patients"),
                                 choices=getPatientId(), server=TRUE)
            updateSelectizeInput(session, paste0("groupRows", "Samples"),
                                 choices=getSampleId(), server=TRUE)
        } else if (type == "ASevents") {
            updateSelectizeInput(session, paste0("groupRows", "ASevents"),
                                 choices=getASevents(), server=TRUE)
            updateSelectizeInput(session, paste0("groupRows", "Genes"),
                                 choices=getGenes(), server=TRUE)
        }
    })
    
    # Create new group(s)
    createGroupOptions <- function(id, hasAttributes=TRUE) {
        collapse <- "groupCollapse"
        panel    <- "Data groups"
        
        observeEvent(input[[paste0("createGroupRows", id)]], {
            createGroup(session, input, output, id, type="Index/Identifier")
            updateCollapse(session, collapse, open=panel)
        })
        
        if (!hasAttributes) return(NULL)
        
        observeEvent(input[[paste0("createGroupAttribute", id)]], {
            createGroup(session, input, output, id, type="Attribute")
            updateCollapse(session, collapse, open=panel)
        })
        
        observeEvent(input[[paste0("createGroupSubset", id)]], {
            createGroup(session, input, output, id, type="Subset")
            updateCollapse(session, collapse, open=panel)
        })
        
        # Update available attributes to suggest in the subset expression
        output[[paste0("groupExpressionSuggestions", id)]] <- renderUI({
            if (id == "Patients")
                attrs <- getPatientAttributes()
            else if (id == "Samples")
                attrs <- getSampleAttributes()
            textSuggestions(ns(paste0("groupExpression", id)), attrs)
        })
        
        observeEvent(input[[paste0("createGroupRegex", id)]], {
            createGroup(session, input, output, id, type="Regex")
            updateCollapse(session, collapse, open=panel)
        })
    }
    
    if (type == "Samples") {
        observe( createGroupOptions("Patients") )
        observe( createGroupOptions("Samples") )
    } else if (type == "ASevents") {
        observe( createGroupOptions("ASevents") )
        observe( createGroupOptions("Genes") )
    }
    
    # Render data table with groups
    output$groupsTable <- renderDataTable({
        dataTable <- showGroupsTable(type)
        if (!is.null(dataTable)) {
            plus <- '<i class="fa fa-plus-circle" aria-hidden="true"></i>'
            return(cbind(' ' = plus, dataTable))
        }
        
        if (type == "Samples") {
            elems <- c("Patients", "Samples")
        } else if (type == "ASevents") {
            elems <- c("ASevents", "Genes")
        }
        
        cols <- c("Names", elems, "Subset", "Input")
        mf <- matrix(ncol=5, dimnames=list(NA, cols))
        return(mf[-1, ])
    }, style="bootstrap", escape=FALSE, server=TRUE, rownames=FALSE,
    options=list(
        pageLength=10, lengthChange=FALSE, scrollX=TRUE, ordering=FALSE,
        columnDefs = list(
            list(orderable=FALSE, className='details-control', targets=0)),
        language=list(zeroRecords="No groups available to display")),
    callback = JS(
        "var plus  = '<i class=\"fa fa-plus-circle\" aria-hidden=\"true\"></i>';",
        "var minus = '<i class=\"fa fa-minus-circle\" aria-hidden=\"true\"></i>';",
        "var cols = table.columns()[0].slice(-2);",
        "table.columns(cols).visible(false, false);",
        "table.columns.adjust().draw(false);",
        "var format = function(d) {
        return '<table class=\"table table-details\" border=\"0\">'+
        '<tr>'+
        '<td>Subset:</td>'+
        '<td>'+ d[d.length - 2] +'</td>'+
        '</tr>'+
        '<tr>'+
        '<td>Input:</td>'+
        '<td>' + d[d.length - 1] + '</td>'+
        '</tr>'+
        '</table>';
};",
        "table.on('click', 'td.details-control', function() {
        var td = $(this),
        row = table.row(td.closest('tr'));
        if (row.child.isShown()) {
        row.child.hide();
        td.html(plus);
        } else {
        row.child( format(row.data()), 'no-padding' ).show();
        td.html(minus);
        }
        });"))
    
    # Remove selected groups
    removeId <- "removeGroups"
    operateOnGroups(input, session, operation="remove", buttonId=removeId,
                    type=type)
    
    # Rename selected groups
    renameId <- "renameGroupName"
    operateOnGroups(input, session, operation="rename", buttonId=renameId,
                    type=type)
    
    # Merge selected groups
    mergeId <- "mergeGroups"
    operateOnGroups(input, session, operation="union", buttonId=mergeId, 
                    symbol=" \u222A ", type=type)
    
    # Intersect selected groups
    intersectId <- "intersectGroups"
    operateOnGroups(input, session, operation="intersect", type=type, 
                    buttonId=intersectId, symbol=" \u2229 ")
    
    # The following set operations are organised inside a "More" button
    moreId <- "moreSetOperations"
    
    # Complement of selected group(s)
    complementId <- "complementGroups"
    operateOnGroups(input, session, operation="complement", type=type, 
                    buttonId=complementId, symbol="U \u005c ")
    
    # Subtract elements from selected group
    subtractId  <- "subtractGroups"
    operateOnGroups(input, session, operation="subtract", type=type, 
                    buttonId=subtractId, symbol=" \u005c ")
    
    subtract2Id <- "subtract2Groups"
    operateOnGroups(input, session, operation="subtract2", type=type, 
                    buttonId=subtract2Id, symbol=" \u005c ")
    
    # Symmetric difference of selected group(s)
    symDiffId <- "symDiffGroups"
    operateOnGroups(input, session, operation="symDiff", type=type, 
                    buttonId=symDiffId, symbol=" \u2A01 ")
    
    saveSelectedGroupsId <- "saveSelectedGroups"
    saveAllGroupsId      <- "saveAllGroups"
    loadGroupsId         <- "loadGroups"
    testGroupsId         <- "testGroups"
    
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
        selGroupsButton  <- getListId(saveSelectedGroupsId)
        allGroupsButton  <- getListId(saveAllGroupsId)
        testGroupsButton <- getListId(testGroupsId)
        
        # Disable data group export when no groups were yet created
        groups <- getGroups(type, complete=TRUE)
        toggleState(allGroupsButton, !is.null(groups) && nrow(groups) > 0)
        
        selectedRows <- length(input$groupsTable_rows_selected)
        # Number of groups selected to enable each set operation
        # - complement: >= 0
        # - remove:     >= 1
        # - subtract:   2
        # - merge:      >= 2
        # - intersect:  >= 2
        # - sym diff:   >= 2
        # Also, group independence testing: >= 1
        
        if (selectedRows >= 1) {
            enable(removeButton)
            enable(selGroupsButton)
            enable(testGroupsButton)
        } else {
            disable(removeButton)
            disable(selGroupsButton)
            disable(testGroupsButton)
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
    
    # Prepare interface for when only one group is selected
    observe({
        selected <- input$groupsTable_rows_selected
        
        # Check table selection only after the table is updated
        allRows <- isolate(input$groupsTable_rows_all)
        groups  <- isolate(getGroups(type, complete=TRUE))
        updated <- length(allRows) == nrow(groups)
        
        if (length(selected) == 1 && updated) {
            show("singleGroupSelected", anim=TRUE, time=0.2)
            
            colour <- groups[[selected, "Colour"]]
            updateColourInput(session, "groupColour", value=colour)
        } else {
            hide("singleGroupSelected", anim=TRUE, time=0.2)
        }
    })
    
    # Set colour of selected group
    observeEvent(input[["setGroupColour-button"]], {
        selected <- input$groupsTable_rows_selected
        groups   <- getGroups(type, complete=TRUE)
        groups[[selected, "Colour"]] <- input$groupColour
        setGroups(type, groups)
    })
    
    # Disable rename button if no new name was given
    observe({
        renameButton <- paste(renameId, "button", sep="-")
        toggleState(renameButton, 
                    !is.null(input$groupName) && input$groupName != "")
    })
    
    # Remove all groups
    observeEvent(input[["removeAll-button"]], setGroups(type, NULL))
    
    # Export groups to file
    exportGroupsToFile <- function(groups, file, type, match) {
        if (type == "Samples") {
            first  <- "Patients"
            second <- "Samples"
        } else if (type == "ASevents") {
            first  <- "Genes"
            second <- "ASevents"
        }
        
        data <- NULL
        for (g in seq(nrow(groups))) {
            group       <- groups[g, ]
            name        <- group$Names
            firstGroup  <- group[[first]]
            secondGroup <- group[[second]]
            
            if (!is.null(secondGroup) && length(secondGroup) > 0) {
                if (!is.null(firstGroup) && length(firstGroup) > 0) {
                    # Match first and second elements available in the group
                    matchingElems <- match[names(match) %in% secondGroup]
                    matchingElems[!matchingElems %in% firstGroup] <- NA
                } else {
                    matchingElems <- NA
                }
                data <- rbind(data, 
                              cbind(name, names(matchingElems), matchingElems))
                
                # Get remaining first elements
                firstGroup <- firstGroup[!firstGroup %in% matchingElems]
                if (length(firstGroup) > 0)
                    data <- rbind(data, cbind(name, NA, firstGroup))
            } else {
                data <- rbind(data, cbind(name, NA, firstGroup))
            }
        }
        second <- gsub("ASevent", "AS event", second)
        strs <- paste(gsub("s$", "", c(second, first)), "ID")
        colnames(data) <- c("Group", strs)
        rownames(data) <- NULL
        write.table(data, file, quote=FALSE, row.names=FALSE, sep="\t")
    }
    
    output[["saveSelectedGroups-button"]] <- downloadHandler(
        filename = function() {
            paste0("psichomics groups ", getCategory(), ".txt")
        }, content = function(file) {
            # Get selected groups
            selected <- sort(input$groupsTable_rows_selected)
            groups   <- getGroups(type, complete=TRUE)
            groups   <- groups[selected, , drop=FALSE]
            
            if (type == "Samples")
                match <- getClinicalMatchFrom("Inclusion levels")
            else if (type == "ASevents" && !is.null(getASevents()) )
                match <- getGenesFromSplicingEvents(getASevents())
            else
                match <- NULL
            exportGroupsToFile(groups, file, type, match)
        }
    )
    
    output[["saveAllGroups-button"]] <- downloadHandler(
        filename = function() {
            paste0("psichomics groups ", getCategory(), ".txt")
        }, content = function(file) {
            groups <- getGroups(type, complete=TRUE)
            if (type == "Samples")
                match <- getClinicalMatchFrom("Inclusion levels")
            else if (type == "ASevents" && !is.null(getASevents()) )
                match <- getGenesFromSplicingEvents(getASevents())
            else
                match <- NULL
            exportGroupsToFile(groups, file, type, match)
        }
    )
    
    # Import groups from file
    observeEvent(input[["loadGroups-button"]], {
        importGroupsFrom <- function (file, first, second, type, match) {
            # Check if file contains group information
            format <- NULL
            if (type == "Samples")
                strs <- paste(c("Sample", "Patient"), "ID")
            else if (type == "ASevents")
                strs <- paste(c("Gene", "AS event"), "ID")
            format$check <- c("Group", strs)
            format$colNames   <- 1
            format$ignoreRows <- 1
            
            head <- fread(file, header=FALSE, nrows=6, stringsAsFactors=FALSE, 
                          data.table=FALSE)
            
            isGroupFile <- checkFileFormat(format, head)
            if (isGroupFile) {
                groups     <- loadFile(format, file)
                firstElem  <- as.character(groups[[strs[2]]])
                secondElem <- as.character(groups[[strs[1]]])
                groups     <- as.character(groups$Group)
                groups[is.na(groups)] <- "NA" # Retrieve groups named NA
                
                if (length(firstElem) == 0 && length(secondElem) == 0) 
                    return(NULL)
                
                # Prepare patients or splicing events per group
                patientGroupNames <- NULL
                if (length(firstElem) > 0) {
                    matched  <- firstElem %in% first
                    firstElem <- split(firstElem[matched], groups[matched])
                    firstElem <- lapply(firstElem, unique)
                    patientGroupNames <- names(firstElem)
                }
                
                # Prepare samples or genes per group
                sampleGroupNames <- NULL
                if (length(secondElem) > 0) {
                    matched <- secondElem %in% second
                    secondElem <- split(secondElem[matched], groups[matched])
                    sampleGroupNames <- names(secondElem)
                }
                
                # Prepare matrix with imported groups
                groupNames <- unique(c(patientGroupNames, sampleGroupNames))
                imported   <- matrix(nrow=length(groupNames), ncol=3)
                colnames(imported) <- c("Names", "Subset", "Input")
                imported[ , "Names"]   <- groupNames
                rownames(imported)     <- groupNames
                
                info <- paste("Groups loaded from file", file)
                imported[ , "Subset"] <- info
                imported[ , "Input"]  <- info
                
                fillEmptyGroups <- function(g, vec) {
                    group <- vec[[g]]
                    if (is.null(group))
                        return(character(0))
                    else
                        return(group)
                }
                
                if (!is.null(first)) {
                    if (length(firstElem) > 0) {
                        firstElem <- lapply(groupNames, fillEmptyGroups, 
                                            firstElem)
                        firstElem[sapply(firstElem, is.null)] <- character(0)
                        names(firstElem) <- groupNames
                    } else if (!is.null(match)) {
                        # Provide patients/AS events based on the samples/genes
                        secondTofirst <- function(i, match) {
                            m <- match[i]
                            return(unique(m[!is.na(m)]))
                        }
                        firstElem <- lapply(secondElem, secondTofirst, match)
                    } else {
                        firstElem <- list(character(0))
                    }
                    firstElem <- lapply(firstElem, unique)
                    imported <- cbind(imported, "YYY"=firstElem)
                    
                    if (type == "Samples")       firstStr <- "Patients"
                    else if (type == "ASevents") firstStr <- "ASevents"
                    colnames(imported)[colnames(imported) == "YYY"] <- firstStr
                }
                
                if (!is.null(second)) {
                    if (length(secondElem) > 0) {
                        secondElem <- lapply(groupNames, fillEmptyGroups, 
                                             secondElem)
                        secondElem <- lapply(secondElem, unique)
                        names(secondElem) <- groupNames
                    } else {
                        secondElem <- list(character(0))
                    }
                    imported <- cbind(imported, "ZZZ"=secondElem)
                    
                    if (type == "Samples")       secondStr <- "Samples"
                    else if (type == "ASevents") secondStr <- "Genes"
                    colnames(imported)[colnames(imported) == "ZZZ"] <- secondStr
                }
                
                # Keep original order if possible
                if ( all(unique(groups) %in% rownames(imported)) )
                    imported <- imported[unique(groups), ]    
                
                return(imported)
            } else {
                return(NULL)
            }
        }
        
        groupsFile <- fileBrowser()
        
        isolate({
            if (type == "Samples") {
                first  <- getPatientId()
                second <- getSampleId()
                match  <- getClinicalMatchFrom("Inclusion levels")
            } else if (type == "ASevents") {
                first  <- getASevents() 
                second <- getGenes()
                match  <- getGenesFromSplicingEvents(first)
            }
        })
        
        imported <- tryCatch(
            importGroupsFrom(groupsFile, first, second, type, match), 
            error=return)
        
        if (!is.null(imported) && !is(imported, "error"))
            appendNewGroups(type, imported)
    })
    
    if (type == "Samples") {
        # Test group indepedence for patient/sample groups only
        observeEvent(input[["testGroups-button"]], {
            isolate({
                # Get selected groups
                selected <- sort(input$groupsTable_rows_selected)
                groups   <- getGroups("Samples", complete=TRUE)
                groups   <- groups[selected, , drop=FALSE]
                
                # Get subject and sample attributes
                clinical   <- getClinicalData()
                sampleAttr <- getSampleInfo()
            })
            
            res <- NULL
            # Test against categorical variables from clinical attributes
            if (!is.null(clinical)) {
                reference  <- setNames(groups[ , "Patients"], rownames(groups))
                categories <- parseCategoricalGroups(clinical)
                res        <- rbind(res, testGroupIndependence(
                    reference, categories, clinical))
            }
            closeProgress()
            
            # Test against categorical variables from sample attributes
            if (!is.null(sampleAttr)) {
                reference  <- setNames(groups[ , "Samples"], rownames(groups))
                categories <- parseCategoricalGroups(sampleAttr)
                res        <- rbind(res, testGroupIndependence(
                    reference, categories, sampleAttr))
            }
            closeProgress()
            setGroupIndependenceTesting(res)
        })
        
        output$groupIndependenceTestingPlot <- renderPlot({
            tests <- getGroupIndependenceTesting()
            if (!is.null(tests)) plotGroupIndependence(tests, top=25, 
                                                       textSize=15)
        })
        
        output[["groupIndependenceTestingPlot-tooltip"]] <- renderUI({
            tests <- getGroupIndependenceTesting()
            if (!is.null(tests)) {
                df    <- tests
                top   <- 25
                hover <- input[["groupIndependenceTestingPlot-hover"]]
                x     <- "Attributes"
                y     <- "Reference"
                
                # Sort and select attributes with the lowest p-values
                ord <- order(df[["Adjusted p-value"]])
                attrs <- unique(df[["Attributes"]][ord])[seq(top)]
                df <- df[df[["Attributes"]] %in% attrs, ]
                # Avoid log of zeroes
                df$pvalue <- df[["Adjusted p-value"]] + .Machine$double.xmin
                
                # Modify data frame according to data representation
                row <- match(df$Reference, unique(df$Reference))
                df2 <- data.frame(row=row, col=seq(length(row) / max(row)))
                rownames(df2) <- rownames(df)
                df2 <- cbind(df2, df)
                
                point <- nearPoints(df2, hover, threshold=10, maxpoints=1, 
                                    addDist=TRUE, xvar="col", yvar="row")
                if (nrow(point) == 0) return(NULL)
                
                # Calculate point position inside the image as percent of total 
                # dimensions from left (horizontal) and from top (vertical)
                xDomain    <- hover$domain$right - hover$domain$left
                left_pct   <- (hover$x - hover$domain$left) / xDomain
                yDomain    <- hover$domain$top - hover$domain$bottom
                bottom_pct <- (hover$y - hover$domain$bottom) / yDomain
                
                # Calculate distance from left and bottom in pixels
                xRange    <- hover$range$right - hover$range$left
                left_px   <- hover$range$left + left_pct * xRange + 40
                yRange    <- hover$range$bottom - hover$range$top
                bottom_px <- hover$range$bottom + bottom_pct * yRange - 140
                
                trItem <- function(key, value) 
                    tags$tr(tags$td(tags$b(key)), tags$td(value))
                
                # Prepare contigency table
                ref <- as.character(point[["Reference"]])
                ref <- gsub("vs others", "", ref, fixed=TRUE)
                contigencyTable <- point[["Contigency table"]][[1]]
                rownames(contigencyTable) <- c(ref, "Others")
                contigencyTable <- table2html(contigencyTable)
                
                wellPanel(
                    class="well-sm",
                    style=paste0(
                        "position: absolute; z-index: 100;",
                        "background-color: rgba(245, 245, 245, 0.85); ",
                        "left:", left_px, "px; bottom:", bottom_px, "px;"),
                    tags$table(
                        class="table table-condensed", 
                        style="margin-bottom: 0; width: auto;",
                        tags$thead( trItem("Attribute", point$Attributes) ),
                        tags$tbody(
                            trItem("p-value", signifDigits(point$pvalue)),
                            trItem("p-value (BH)",
                                   signifDigits(point$`Adjusted p-value`)))),
                    tags$div(style="padding-left: 5px;",
                             tags$b("Contigency table"), contigencyTable))
            }
        })
    }
}

#' @rdname appServer
#'
#' @inheritParams operateOnGroups
groupsServer <- function(input, output, session) {
    callModule(groupManipulation, "sampleGroupModule",  "Samples")
    callModule(groupManipulation, "ASeventGroupModule", "ASevents")
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
        group <- getGroups("Samples", complete=TRUE)
        # Ignore if there are no groups
        if (is.null(group)) return(NULL)
        
        samples     <- getSampleId()
        patients    <- getPatientId()
        match       <- getClinicalMatchFrom("Inclusion levels")
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
            setGroups("Samples", group)
        } else if ( !is.null(samples) && !showSamples && !is.null(patients) &&
                    !is.null(match) && "Patients" %in% colnames(group)) {
            # Update groups if previously made with patients only
            samples <- getMatchingSamples(group[ , "Patients"], samples, 
                                          patients, match=match)
            group <- cbind(group, "Samples"=samples)
            setGroups("Samples", group)
        }
    })
    
    # Create groups by sample types when loading TCGA data
    observe({
        sampleInfo <- getSampleInfo()
        if (!is.null(sampleInfo) && any(grepl("^TCGA", rownames(sampleInfo)))) {
            new    <- createGroupByAttribute("Sample types", sampleInfo)
            groups <- cbind("Names"=names(new), "Subset"="Attribute",
                            "Input"="Sample types", "Samples"=new)

            # Match samples with patients (if loaded)
            patients <- isolate(getPatientId())
            if (!is.null(patients)) {
                indiv <- lapply(new, function(i)
                    unname(getPatientFromSample(i, patientId=patients)))
                groups <- cbind(groups[ , 1:3, drop=FALSE], "Patients"=indiv,
                                groups[ ,   4, drop=FALSE])
            }
            
            if (!is.null(groups)) 
                isolate( appendNewGroups("Samples", groups, clearOld=TRUE) )
        }
    })
}

#' @rdname selectGroupsUI
#'  
#' @inheritParams getGroups
#' @param input Shiny input
#' @param filter Character: get groups only if they are present in this argument
#' (if TCGA-styled gene symbols, they will be "converted" to gene symbols alone)
#' 
#' @return \code{getSelectedGroups}: List with selected groups (or NULL if no
#' groups were selected)
getSelectedGroups <- function(input, id, type, filter=NULL) {
    selection <- input[[paste0(id, "Selection")]]
    noGroups  <- !is.null(selection) && selection == "noGroups"
    selected  <- input[[id]]
    
    if ( noGroups || is.null(selected) ) {
        # User selects no groups (either explicitly or not)
        groups <- NULL
    } else {
        groups <- getGroups(type)
        colour <- attr(groups, "Colour")
        groups <- groups[selected]
        
        if (!is.null(filter)) {
            names(filter) <- filter
            if (type == "Genes") {
                # Compare TCGA-styled genes
                original <- filter
                filter   <- gsub("\\|.*", "", filter)
                filter[filter == "?"] <- original[filter == "?"]
                names(filter) <- original
            }
            groups <- lapply(groups, function(i) names(filter)[filter %in% i])
        }
        attr(groups, "Colour") <- colour[selected]
    }
    return(groups)
}


# Multiple group independence testing -------------------------------------
# Inspiration from https://rud.is/projects/facetedheatmaps.html

#' Parse categorical columns in a data frame
#' 
#' Retrieve elements grouped by their unique group based on each categorical 
#' column
#' 
#' @param df Data frame
#' 
#' @seealso \code{\link{testGroupIndependence}} and 
#' \code{\link{plotGroupIndependence}}
#' 
#' @return List of lists containing values based on rownames of \code{df}
#' @export
#' @examples 
#' df <- data.frame("race"=c("caucasian", "caucasian", "asian"),
#'                  "gender"=c("male", "female", "male"))
#' rownames(df) <- paste("patient", 1:3)
#' parseCategoricalGroups(df)
parseCategoricalGroups <- function(df) {
    isCategorical <- sapply(seq(ncol(df)), function(i) is.factor(df[ , i]))
    categories <- df[ , isCategorical, drop=FALSE]
    
    groups <- lapply(colnames(categories), function(col) {
        divisions <- as.character(df[ , col])
        divisions[is.na(divisions)] <- "NA" # Do not dismiss missing values
        split(rownames(df), divisions)
    })
    names(groups) <- colnames(categories)
    return(groups)
}

#' Multiple independence tests between a reference group and list of groups
#' 
#' Uses Fisher's exact test.
#' 
#' @param ref Character: identifier of elements in reference group
#' @param groups List of characters: list of groups where each element contains
#' the identifiers of respective elements
#' @param elements Character: all patient identifiers
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#' 
#' @importFrom stats p.adjust fisher.test
#' 
#' @details The following methods for p-value adjustment are supported by using 
#' the respective string in the \code{pvalueAdjust} argument:
#' \itemize{
#'     \item{\code{none}: Do not adjust p-values}
#'     \item{\code{BH}: Benjamini-Hochberg's method (false discovery rate)}
#'     \item{\code{BY}: Benjamini-Yekutieli's method (false discovery rate)}
#'     \item{\code{bonferroni}: Bonferroni correction (family-wise error rate)}
#'     \item{\code{holm}: Holm's method (family-wise error rate)}
#'     \item{\code{hochberg}: Hochberg's method (family-wise error rate)}
#'     \item{\code{hommel}: Hommel's method (family-wise error rate)}
#' }
#' 
#' @return Returns a \code{groupIndependenceTest} object: a list where each 
#' element is a list containing:
#' \item{attribute}{Name of the original groups compared against the reference
#' groups}
#' \item{table}{Contigency table used for testing}
#' \item{pvalue}{Fisher's exact test's p-value}
testSingleIndependence <- function(ref, groups, elements, pvalueAdjust="BH") {
    # Number of intersections between reference and groups of interest
    updateProgress("Calculating intersections", 
                   detail="Reference versus categorical groups")
    getIntersectionsNumber <- function(reference, attribute)
        length(intersect(reference, attribute))
    intersections1 <- lapply(groups, sapply, getIntersectionsNumber, ref)
    
    # Number of intersections between complement and groups of interest
    updateProgress("Calculating intersections", 
                   detail="Complement versus categorical groups")
    getComplementOf <- function(ref, elements) setdiff(elements, ref)
    complement <- getComplementOf(ref, elements)
    intersections2 <- lapply(groups, sapply, getIntersectionsNumber, complement)
    
    # Fisher's exact method for group independence testing
    updateProgress("Performing multiple Fisher's exact tests")
    groupIndependenceTesting <- function(i, intersections1, intersections2) {
        mat    <- rbind(intersections1[[i]], intersections2[[i]])
        updateProgress("Performing multiple Fisher's exact tests", 
                       console=FALSE)
        if ( all(sort(unique(as.numeric(mat))) %in% 0:1) ) {
            pvalue <- 1
        } else {
            pvalue <- tryCatch(fisher.test(mat)$p.value, error=return)
            if (is(pvalue, "error")) pvalue <- 1
        }
        col <- names(intersections1)[[i]]
        return(list(attribute=col, table=mat, pvalue=pvalue))
    }
    res <- lapply(seq(intersections1), groupIndependenceTesting, 
                  intersections1, intersections2)
    names(res) <- sapply(res, "[[", "attribute")
    
    adjusted <- p.adjust(sapply(res, "[[", "pvalue"), pvalueAdjust)
    for (i in seq(res)) res[[i]]$pvalueAdjusted <- adjusted[[i]]
    class(res) <- c(class(res), "groupIndependenceTest")
    return(res)
}

#' Multiple independence tests between reference groups and list of groups
#' 
#' Test multiple contigency tables comprised by two groups (one reference group
#' and another containing remaing elements) and provided groups.
#' 
#' @param ref List of character: list of groups where each element contains the
#' identifiers of respective elements
#' @param groups List of characters: list of groups where each element contains
#' the identifiers of respective elements
#' @param elements Character: all available elements (if a data frame is given,
#' its rownames will be used)
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#' 
#' @details The following methods for p-value adjustment are supported by using 
#' the respective string in the \code{pvalueAdjust} argument:
#' \itemize{
#'     \item{\code{none}: Do not adjust p-values}
#'     \item{\code{BH}: Benjamini-Hochberg's method (false discovery rate)}
#'     \item{\code{BY}: Benjamini-Yekutieli's method (false discovery rate)}
#'     \item{\code{bonferroni}: Bonferroni correction (family-wise error rate)}
#'     \item{\code{holm}: Holm's method (family-wise error rate)}
#'     \item{\code{hochberg}: Hochberg's method (family-wise error rate)}
#'     \item{\code{hommel}: Hommel's method (family-wise error rate)}
#' }
#' 
#' @return \code{multiGroupIndependenceTest} object, a data frame containing:
#' \item{attribute}{Name of the original groups compared against the reference
#' groups}
#' \item{table}{Contigency table used for testing}
#' \item{pvalue}{Fisher's exact test's p-value}
#' 
#' @seealso \code{\link{parseCategoricalGroups}} and 
#' \code{\link{plotGroupIndependence}}
#' 
#' @export
#' @examples 
#' elements <- paste("patients", 1:10)
#' ref      <- elements[5:10]
#' groups   <- list(race=list(asian=elements[1:3],
#'                            white=elements[4:7],
#'                            black=elements[8:10]),
#'                  region=list(european=elements[c(4, 5, 9)],
#'                              african=elements[c(6:8, 10)]))
#' groupTesting <- testGroupIndependence(ref, groups, elements)
#' # View(groupTesting)
testGroupIndependence <- function(ref, groups, elements, pvalueAdjust="BH") {
    if (is.data.frame(elements)) elements <- rownames(elements)
    
    if (is.list(ref)) {
        updateProgress("Group independence testing", 
                       divisions=3 + length(ref) * length(groups))
        obj <- lapply(seq(ref), function(i) {
            updateProgress("Group independence testing", detail=names(ref)[[i]])
            testSingleIndependence(ref[[i]], groups, elements)
        })
    } else {
        updateProgress("Group independence testing", 
                       divisions=3 + length(groups))
        obj <- list("Reference"=testSingleIndependence(ref, groups, elements))
    }
    names(obj) <- names(ref)
    
    df <- NULL
    for (g in seq(obj)) {
        name <- names(obj)[[g]]
        if (is.null(name) || identical(name, "")) name <- "Reference"
        
        pvalue   <- sapply(obj[[g]], "[[", "pvalue")
        adjusted <- sapply(obj[[g]], "[[", "pvalueAdjusted")
        mat <- lapply(obj[[g]], "[[", "table")
        new <- data.frame("Reference"=paste(name, "vs others"), 
                          "Attributes"=names(pvalue),
                          "p-value"=pvalue,
                          "Adjusted p-value"=adjusted,
                          "Contigency table"=I(mat))
        df <- rbind(df, new)
    }
    colnames(df) <- c("Reference", "Attributes", "p-value", "Adjusted p-value",
                      "Contigency table")
    class(df) <- c(class(df), "multiGroupIndependenceTest")
    return(df)
}

#' Plot -log10(p-values) of the results obtained after multiple group 
#' independence testing
#' 
#' @param groups \code{multiGroupIndependenceTest} object (obtained after 
#' running \code{\link{testGroupIndependence}})
#' @param top Integer: number of attributes to render
#' @param textSize Integer: size of the text
#' @param colourLow Character: name or HEX code of colour for lower values
#' @param colourMid Character: name or HEX code of colour for middle values
#' @param colourHigh Character: name or HEX code of colour for higher values
#' @param colourMidpoint Numeric: midpoint to identify middle values
#' 
#' @importFrom ggplot2 ggplot aes geom_tile theme coord_equal unit element_blank
#' element_text scale_fill_gradient2
#' 
#' @seealso \code{\link{parseCategoricalGroups}} and 
#' \code{\link{testGroupIndependence}}
#' 
#' @return \code{ggplot} object
#' @export
#' 
#' @examples 
#' elements <- paste("patients", 1:50)
#' ref      <- elements[10:50]
#' groups   <- list(race=list(asian=elements[1:3],
#'                            white=elements[4:7],
#'                            black=elements[8:10]),
#'                  region=list(european=elements[c(4, 5, 9)],
#'                              african=elements[c(6:8, 10:50)]))
#' groupTesting <- testGroupIndependence(ref, groups, elements)
#' plotGroupIndependence(groupTesting)
plotGroupIndependence <- function(groups, top=50, textSize=10, 
                                  colourLow="lightgrey", colourMid="blue", 
                                  colourHigh="orange", colourMidpoint=150) {
    # Sort and select attributes with the lowest p-values
    df <- groups
    ord <- order(df[["Adjusted p-value"]])
    attrs <- head(unique(df[["Attributes"]][ord]), top)
    df <- df[df[["Attributes"]] %in% attrs, ]
    # Avoid log of zeroes
    df$pvalue <- df[["Adjusted p-value"]] + .Machine$double.xmin
    
    plot <- ggplot(df, aes_string(x="Attributes", y="Reference")) + 
        geom_tile(aes_string(fill="-log10(pvalue)"), color="white", size=0.1) +
        coord_equal() +
        theme(axis.title.x=element_text(size=textSize),
              axis.text.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_text(size=textSize),
              axis.ticks=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              legend.position="bottom",
              legend.title=element_text(size=textSize),
              legend.text=element_text(size=textSize - 4),
              legend.key.width=unit(0.5, "cm"),
              legend.key.height=unit(0.2, "cm")) + 
        scale_fill_gradient2(low=colourLow, mid=colourMid, high=colourHigh, 
                             midpoint=colourMidpoint)
    return(plot)
}

attr(groupsUI, "loader") <- "app"
attr(groupsServer, "loader") <- "app"