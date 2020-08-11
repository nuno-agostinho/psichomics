#' Group selection
#'
#' Group selection interface and logic
#'
#' @param id Character: identifier
#' @param label Character: \code{selectize} label
#' @inheritParams getGroups
#' @param placeholder Character: \code{selectize} placeholder
#' @param noGroupsLabel Character: label to explicitly allow to select no groups
#' (if \code{NULL}, this option is not displayed to the user)
#' @param groupsLabel Character: label to explicitly allow to select groups
#' (only required if \code{noGroupsLabel} is not \code{NULL})
#' @param maxItems Numeric: maximum number of groups to select
#' @param returnAllDataLabel Character: label to allow to return data outside
#' selected groups as belonging to an outside group (if \code{NULL}, this option
#' is not displayed to the user)
#' @param returnAllDataValue Boolean: default value to whether return all data
#' or not (only required if \code{returnAllDataLabel} is not \code{NULL})
#'
#' @importFrom shiny fluidRow column uiOutput selectizeInput actionButton
#' radioButtons actionLink downloadLink
#'
#' @note To allow the user to (explicitly) select no groups, pass the
#' \code{noGroupsLabel} and \code{groupsLabel} arguments.
#'
#' @return \code{selectGroupsUI}: Interface for group selection
#' @keywords internal
selectGroupsUI <- function (
    id, label, type, placeholder="Type to search groups",
    noGroupsLabel=NULL, groupsLabel=NULL, maxItems=NULL,
    returnAllDataLabel=NULL, returnAllDataValue=FALSE) {

    editId <- paste0(id, "Edit")
    onItemAdd <- sprintf(
        "function(value, $item) {
            var editLabel = 'Create or edit groups...';
            if (value === editLabel) {
                showGroups('%s');
                this.removeItem(editLabel);
                this.blur();
            }
        }", type)
    groupSelect <- selectizeInput(
        id, label, choices=NULL, multiple=TRUE, width="auto", options=list(
            plugins=list('remove_button', 'drag_drop'), maxItems=maxItems,
            searchField=list("value", "label"), placeholder=placeholder,
            onItemAdd=I(onItemAdd), render=I(
                '{option: renderGroupSelection, item: renderGroupSelection}')))

    if ( !is.null(label) ) {
        if ( is.null(noGroupsLabel) ) {
            label <- groupSelect$children[[1]]
        } else {
            # Use label in radio buttons instead
            radioLabel <- groupSelect$children[[1]]
            label <- NULL
        }
        groupSelect$children[[1]] <- NULL
    }
    select <- tagList(label, groupSelect)

    if ( !is.null(returnAllDataLabel) ) {
        noGroupsId <- paste0(id, "ShowAllData")
        select <- tagList(select, checkboxInput(
            noGroupsId, returnAllDataLabel, value=returnAllDataValue))
    }
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
#' (if \code{NULL}, all groups will be pre-selected)
#'
#' @importFrom shinyjs enable disable onclick toggleClass runjs
#'
#' @return \code{selectGroupsServer}: Server logic for group selection
selectGroupsServer <- function(session, id, type, preference=NULL) {
    ns     <- session$ns
    input  <- session$input
    output <- session$output

    modalId <- paste0(id, "Modal")
    # Update groups shown in the selectize element
    observe({
        groupTable <- getGroups(type, complete=TRUE)
        groups <- rownames(groupTable)
        if (is.null(groups)) {
            groups <- list()
        } else {
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

            ns1 <- ns2 <- NULL
            if (!is.null(elem1Number)) {
                elem1 <- gsub("ASevents", "AS events", elem1)
                elem1 <- ifelse(elem1Number == 1, gsub(".$", "", elem1), elem1)
                ns1 <- sprintf("%s %s", elem1Number, elem1)
            }
            if (!is.null(elem2Number)) {
                elem2 <- ifelse(elem2Number == 1, gsub(".$", "", elem2), elem2)
                ns2 <- sprintf("%s %s", elem2Number, elem2)
            }

            if (!is.null(ns1) && !is.null(ns2)) {
                ns <- paste(ns1, ns2, sep=", ")
            } else if (!is.null(ns1)) {
                ns <- ns1
            } else if (!is.null(ns2)) {
                ns <- ns2
            }

            names(groups) <- paste(ns, unlist(groupTable[ , "Colour"]))
        }

        selected <- isolate(input[[id]])
        selected <- selected[selected %in% groups]
        if ( is.null(selected) ) selected <- character()
        groups <- c(groups, "Edit"="Create or edit groups...")
        updateSelectizeInput(session, id, choices=groups, selected=selected)
    })
}

#' Interface to manipulate data grouping
#'
#' @param id Character: identifier
#' @param type Character: type of data for each the interface is intended
#'
#' @return HTML elements
#' @keywords internal
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
        identifierUI <- groupById(ns, id)
        if (id %in% c("Patients", "Samples")) {
            navbarMenu(
                title,
                tabPanel("Attribute", groupByAttribute(ns, cols, id, example)),
                tabPanel("Index/Identifier", identifierUI),
                "----", "Advanced options",
                tabPanel("Subset expression", groupByExpression(ns, id)),
                tabPanel("Regular expression", groupByGrep(ns, cols, id)))
        } else if (id == "Genes") {
            navbarMenu(
                title,
                tabPanel("Gene names", identifierUI),
                tabPanel("Pre-made gene lists",
                         groupByPreMadeList(ns, getGeneList(), id)))
        } else {
            tabPanel(title, identifierUI)
        }
    }

    if (type == "Samples") {
        title       <- "By subjects"
        first       <- groupOptions("Patients", title)
        firstAlert  <- tabPanel(
            title, value="NoSubjects", missingData(
                "Subject information",
                "No subject information available.",
                "Please load subject information."))

        title       <- "By samples"
        second      <- groupOptions("Samples", title)
        secondAlert <- tabPanel(
            title, value="NoSamples", missingData(
                "Sample information",
                "No sample information available.",
                "Please load sample information."))
    } else if (type == "ASevents") {
        title      <- "By splicing events"
        first      <- groupOptions("ASevents", title)
        firstAlert <- tabPanel(
            title, value="NoASevents", missingData(
                "Alternative splicing events",
                "No alternative splicing events are available.",
                "Please load or quantify splicing events."))

        title       <- "By genes"
        second      <- groupOptions("Genes", title)
        secondAlert <- tabPanel(
            title, value="NoGenes", missingData(
                "Genes", "No genes are available.",
                "Please load gene expression or load/quantify splicing events.")
        )
    }

    sidebarLayout(
        sidebarPanel(
            uiOutput(ns("alert-side")),
            tabsetPanel(id=ns("groupBy"), type="pills",
                        first, firstAlert, second, secondAlert) ),
        mainPanel(
            uiOutput(ns("alert-main")),
            renderGroupInterface(ns, multiFisherTests) ))
}

#' @rdname appUI
#' @importFrom shinyBS bsCollapse bsCollapsePanel
groupsUI <- function(id, tab) {
    ns <- NS(id)

    tab(icon="object-ungroup", title="Groups", tabsetPanel(
        id="groupsTypeTab",
        tabPanel(
            "Subject and sample groups",
            groupManipulationInput(ns("sampleGroupModule"), "Samples")),
        tabPanel(
            "Splicing event and gene groups",
            groupManipulationInput(ns("ASeventGroupModule"), "ASevents"))))
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
#' @keywords internal
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
        icon=setOperationIcon("complement-AB"), disable=FALSE)

    subtractLink <- operationLink(
        "Subtract elements from upper-selected group",
        helpText("Select two groups for subtraction operations"),
        id=ns(subtractId), icon=setOperationIcon("difference-AB"))

    subtract2Link <- operationLink(
        "Subtract elements from lower-selected group",
        helpText("Select two groups for subtraction operations"),
        id=ns(subtract2Id), icon=setOperationIcon("difference-BA"))

    symDiffLink <- operationLink(
        "Symmetric difference",
        helpText("Select two or more groups for symmetric difference"),
        id=ns(symDiffId), icon=setOperationIcon("symmetric-difference"))

    operations <- div(
        id=ns("setOperations"), class="btn-group",
        style="margin-top: 4px; margin-bottom: 4px;",
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
        icon("user"), class=NULL, "Save selected groups",
        id=ns(saveSelectedGroupsId))

    saveAllGroupsLink <- downloadContent(
        icon("users"), class=NULL, "Save all groups", id=ns(saveAllGroupsId),
        disable=FALSE)

    loadGroupsLink <- operationLink(
        "Load groups", id=ns(loadGroupsId), icon=icon("plus-circle"),
        disable=FALSE)

    saveLoadGroups <- tags$div(
        class="btn-group", role="group",
        style="margin-top: 4px; margin-bottom: 4px;",
        tags$button("Save and load", icon("folder-open"), id=ns(saveLoadId),
                    tags$span(class="caret"),
                    class="btn btn-default dropdown-toggle",
                    "data-toggle"="dropdown", "aria-haspopup"="true",
                    "aria-expanded"="true"),
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
            style="margin-top: 4px; margin-bottom: 4px;",
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
    nameField$children[[1]] <- NULL
    nameField$attribs$style <- "margin: 0; width: auto;"
    nameField$children[[1]][[2]]$style <- "border-radius: 4px 0 0 4px;"

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
        class="alert-info",
        style="margin-top: 10px; margin-bottom: 0px;",
        style="padding-top: 11px; padding-bottom: 11px;",
        fluidRow(
            column(7, style="margin-top: 4px; margin-bottom: 4px;",
                   div(class="input-group", nameField,
                       div(class="input-group-btn", renameButton))),
            column(5, style="margin-top: 4px; margin-bottom: 4px;",
                   div(class="input-group", colourSelector,
                       div(class="input-group-btn", setColourButton)))))

    removeGroups <- tags$div(
        class="btn-group pull-right",
        style="margin-top: 4px; margin-bottom: 4px;",
        removeButton,
        tags$button(type="button", tags$span(class="caret"),
                    class="btn btn-danger dropdown-toggle",
                    "data-toggle"="dropdown",
                    "aria-haspopup"="true",
                    "aria-expanded"="false"),
        tags$ul(class="dropdown-menu",
                style="background-color: #d9534f;",
                style="border-color: #d43f3a;", removeAllLink))

    groupTestId <- "groupIndependenceTestingPlot"
    tagList(
        dataTableOutput(ns("groupsTable")),
        helpText("Click on groups in the table above to select them and",
                 "perform the following actions."),
        removeGroups, operations, saveLoadGroups,
        if (multiFisherTests) groupIndependenceSet,
        # checkboxInput(ns("removeSetsUsed"), "Remove original groups",
        #               value=TRUE)
        hidden(singleGroupSelectedInterface),
        plotOutput(ns(groupTestId), height="200px",
                   hover=hoverOpts(ns(paste0(groupTestId, "-hover")),
                                   delay=50, delayType="throttle")),
        uiOutput(ns(paste0(groupTestId, "-tooltip"))))
}

#' Data grouping interface
#'
#' @param ns Namespace function
#' @param cols Character or list: name of columns to show
#' @param id Character: identifier
#' @param example Character: text to show as an example
#'
#' @return HTML elements
#' @keywords internal
groupByAttribute <- function(ns, cols, id, example) {
    if (!is.null(example)) example <- tagList(" ", example)

    tagList(
        helpText("Automatically create groups according to the unique values",
                 "for the selected attribute.", example),
        selectizeInput(ns(paste0("groupAttribute", id)), "Select attribute",
                       width="auto", choices=cols, options=list(
                           lockOptgroupOrder=TRUE,
                           placeholder="Type to search attributes")),
        uiOutput(ns(paste0("previewGroups", id))),
        actionButton(ns(paste0("createGroupAttribute", id)), "Create groups",
                     class ="btn-primary"))
}

#' @rdname groupByAttribute
#'
#' @param data List: list of groups with elements
#'
#' @importFrom shiny helpText tags
groupByPreMadeList <- function(ns, data, id) {
    cols <- preparePreMadeGroupForSelection(data)

    tagList(
        helpText("Load pre-made, literature-based lists of genes."),
        selectizeInput(ns(paste0("groupAttribute", id)), "Gene list",
                       width="auto", choices=cols, options=list(
                           lockOptgroupOrder=TRUE, placeholder="Gene list")),
        tags$small(uiOutput(ns(paste0("geneListNumber", id)))),
        tags$small(tags$b("Source:"),
                   helpText(textOutput(ns(paste0("geneListSource", id))))),
        actionButton(ns(paste0("loadPreMadeGroup", id)), "Create group",
                     class="btn-primary")
    )
}

#' Prepare list of pre-made groups for a \code{selectize} element
#'
#' @param groups List of list of characters
#'
#' @return List
#' @keywords internal
preparePreMadeGroupForSelection <- function(groups) {
    res <- lapply(groups, names)
    for (ns in names(res)) {
        tmp <- res[[ns]]
        res[[ns]] <- paste(ns, res[[ns]], sep="|||")
        names(res[[ns]]) <- tmp
    }
    return(res)
}

#' Select pre-made groups from a selected item
#'
#' @param groups List of list of characters
#' @param selected Character: selected item
#'
#' @return Elements of selected item
#' @keywords internal
selectPreMadeGroup <- function(groups, selected, genes=NULL) {
    selected <- strsplit(selected, "\\|\\|\\|| ~ ")[[1]]
    first    <- selected[[1]]
    second   <- selected[[2]]
    res      <- groups[[first]][[second]]

    # Filter genes based on those currently available
    if (!is.null(genes)) {
        valid     <- res %in% genes
        total     <- length(res)
        discarded <- sum(!valid)
        perc      <- round(discarded / total * 100)

        res <- res[valid]
        if (discarded > 0) {
            msg <- tagList(tags$b(sprintf("%s out of %s genes (%s%%)",
                                          discarded, total, perc)),
                           "do not match loaded data and were discarded")
            attr(res, "discarded") <- msg
        }
    }
    attr(res, "title") <- sprintf("%s (%s)", second, first)
    attr(res, "citation") <- attr(groups[[first]], "citation")
    return(res)
}

#' @rdname groupByAttribute
#' @importFrom shiny textInput
groupById <- function(ns, id) {
    sid <- gsub("s$", "", id)
    sid <- gsub("ASevent", "Splicing event", sid)

    tagList(
        selectizeInput(
            ns(paste0("groupRows", id)), paste(sid, "indexes or identifiers"),
            choices=NULL, multiple=TRUE, width="auto", options=list(
                create=TRUE, createOnBlur=TRUE, # Allow to add new items
                plugins=list('remove_button'), persist=FALSE,
                placeholder="Type to search identifiers")),
        helpText("Example: ", tags$kbd("1:6, 8, 10:19"), "creates a group with",
                 "items 1 to 6, 8 and 10 to 19. You can also input identifiers",
                 "instead of indexes."),
        textInput(ns(paste0("groupNameRows", id)), "Group label", width="auto",
                  placeholder="Unlabelled group"),
        actionButton(ns(paste0("createGroupRows", id)), "Create group",
                     class="btn-primary")
    )
}

#' @rdname groupByAttribute
#' @importFrom shiny textInput
groupByExpression <- function(ns, id) {
    tagList (
        textInput(ns(paste0("groupExpression", id)), "Subset expression",
                  width="auto", placeholder="Insert subset expression"),
        helpText(
            'Examples: ', tags$ul(
                tags$li(
                    tags$kbd('`X` > 8 & `Y` == "alive"'),
                    ' selects items whose values are higher than 8 for X',
                    ' and equal to ', tags$code("alive"), ' for Y.'),
                tags$li(
                    tags$kbd('grepl("Tumour", `Z`)'),
                    ' selects items whose values contain the word ',
                    tags$code("Tumour"), ' (case sensitive) in Z.'),
                tags$li(
                    tags$kbd('grepl("Tumour", `Z`, ignore.case=TRUE)'),
                    ' selects items whose values contain the word ',
                    tags$code("Tumour"), ' (ignores case) in Z.'),
                tags$li(
                    tags$kbd('!grepl("Tumour", `Z`, ignore.case=TRUE)'),
                    ' selects items whose values ', tags$b('do not'),
                    ' contain the word ', tags$code("Tumour"),
                    ' (ignores case) in Z.'))),
        uiOutput(ns(paste0("groupExpressionSuggestions", id))),
        textInput(ns(paste0("groupNameSubset", id)), "Group label",
                  width="auto", placeholder="Unlabelled group"),
        actionButton(ns(paste0("createGroupSubset", id)), "Create group",
                     class="btn-primary")
    )
}

#' @rdname groupByAttribute
#' @importFrom shiny textInput
groupByGrep <- function(ns, cols, id) {
    tagList (
        textInput(ns(paste0("grepExpression", id)), "Regular expression",
                  width="auto", placeholder="Insert regular expression"),
        selectizeInput(ns(paste0("grepColumn", id)), "Select attribute to GREP",
                       choices=cols, width="auto", options=list(
                           placeholder="Type to search attributes")),
        textInput(ns(paste0("groupNameRegex", id)), "Group label", width="auto",
                  placeholder="Unlabelled group"),
        actionButton(ns(paste0("createGroupRegex", id)), "Create group",
                     class="btn-primary"))
}

#' Prepare to create group according to specific details
#'
#' @inheritParams createGroupFromInput
#' @param output Shiny output
#'
#' @inherit psichomics return
#' @keywords internal
createGroup <- function(session, input, output, id, type, selected=NULL,
                        expr=NULL, groupNames=NULL) {
    removeAlert(output, alertId="alert-side")
    if (id == "Patients") {
        dataset <- getClinicalData()
    } else if (id == "Samples") {
        dataset <- getSampleInfo()
    } else {
        dataset <- NULL
    }
    new <- createGroupFromInput(session, input, output, dataset, id, type,
                                selected, expr, groupNames)
    if (!is.null(new)) appendNewGroups(id, new)
}

#' Assign colours to groups
#'
#' @param new Matrix: groups to which colours will be assigned
#' @param groups Matrix: groups to check which colours are already assigned
#'
#' @return Groups with an added column to state the colour
#' @keywords internal
assignColours <- function(new, groups=NULL) {
    strong <- c("#08419E", "#EF9636", "#D33E6A", "#00C652",
                "#4C71DB", "#8F033B", "#F89CD1", "#05CFC0")
    medium <- c("#7D87B6", "#EFB893", "#E17A90", "#8CD59A",
                "#8696DC", "#9E646F", "#F6C4DF", "#9CDDD5")
    light  <- c("#BEC1D2", "#EAD2C7", "#E7AFBA", "#C6DECA",
                "#B6BBE0", "#D6BBC0", "#F2E1EA", "#D3E7E5")
    colours <- c(strong, medium, light)

    # Avoid setting colours that were previously assigned
    priority <- NULL
    if (!is.null(groups) && "Colour" %in% colnames(groups)) {
        priority <- colours[!colours %in% groups[ , "Colour"]]
    }

    # Repeat default colours when there are many groups
    reps    <- ceiling((nrow(new) - length(priority)) / length(colours))
    colours <- rep(colours, reps)
    colours <- c(priority, colours)[nrow(new):1]

    hasColnames <- !is.null(colnames(new))
    new <- cbind(new, colours)
    if (hasColnames) {
        colnames(new)[ncol(new)] <- "Colour"
    } else {
        colnames(new) <- NULL
    }
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
#' @inherit psichomics return
#' @keywords internal
appendNewGroups <- function(type, new, clearOld=FALSE) {
    # Rename duplicated group names
    if (clearOld)  {
        groups <- NULL
    } else {
        groups <- getGroups(type, complete=TRUE)
    }
    new <- renameGroups(new, groups)

    # Assign colours if not previously assigned
    hasAssignedColours <- "Colour" %in% colnames(new)
    if (!hasAssignedColours) new <- assignColours(new, groups)

    if (clearOld) {
        groups <- new
    } else {
        # Append the new group(s) to the groups already created
        groups <- rbind(new, groups)
    }
    setGroups(type, groups)
}

#' Match subjects and samples in a group
#'
#' @param id Character: identifier
#' @param group Data frame: group
#'
#' @return Data frame with groups containing matching elements
#' @keywords internal
matchGroupSubjectsAndSamples <- function(id, group) {
    subjects <- getSubjectId()
    samples  <- getSampleId()
    match    <- getClinicalMatchFrom("Inclusion levels")

    # Match subjects with samples (or vice-versa)
    if (!is.null(subjects) && !is.null(samples) && !is.null(match)) {
        if (id == "Patients") {
            subjects <- group[ , "Patients"]
            samples <- getSampleFromSubject(subjects, samples, subjects,
                                            match=match)
            group <- cbind(group, "Samples"=samples)
        } else if (id == "Samples") {
            samples2subjects <- function(i, match) {
                m <- match[i]
                return(unique(m[!is.na(m)]))
            }
            subjects <- lapply(group[ , "Samples"], samples2subjects, match)
            group <- cbind(group, "Patients"=subjects)

            lastCol <- ncol(group)
            group   <- group[ , c(seq(lastCol - 2), lastCol, lastCol - 1),
                              drop=FALSE]
        }
    }
    return(group)
}

#' Match AS events and genes in a group
#'
#' @inheritParams matchGroupSubjectsAndSamples
#'
#' @return Data frame with groups containing matching elements
#' @keywords internal
matchGroupASeventsAndGenes <- function(id, group, ASevents) {
    # Match AS events with genes (or vice-versa)
    if (!is.null(ASevents)) {
        if (id == "ASevents") {
            ASevents <- group[ , "ASevents"]
            genes <- lapply(group[ , "ASevents"], parseSplicingEvent,
                            data=ASevents)
            genes <- lapply(genes, "[[", "gene")
            group <- cbind(group, "Genes"=lapply(
                genes, function(i) unique(unlist(i))))
        } else if (id == "Genes") {
            genes <- group[ , "Genes"]
            ASeventGenes <- matchSplicingEventsWithGenes(ASevents)
            filterBasedOnGenes <- function(gene, ASeventGenes)
                ASeventGenes[names(ASeventGenes) %in% gene]

            # Process TCGA gene ID
            genes <- lapply(genes, function(gene) gsub("\\|.*$", "", gene))
            ASevents <- lapply(genes, filterBasedOnGenes, ASeventGenes)
            group <- cbind(group, ASevents)
            group <- group[ , c(1, 2, 3, 5, 4), drop=FALSE]
        }
    } else if (id == "Genes") {
        ASevents <- lapply(seq(nrow(group)), function(i) character(0))
        group <- cbind(group, ASevents)
        group <- group[ , c(1, 2, 3, 5, 4), drop=FALSE]
    }
    return(group)
}

#' Set new groups according to the user input
#'
#' @inheritParams appServer
#' @param dataset Data frame or matrix: dataset of interest
#' @param id Character: identifier of the group selection
#' @param type Character: type of group to create
#' @param selected Character: selected item
#' @param expr Character: expression
#' @param groupNames Character: group names
#'
#' @return Matrix with the group names and respective elements
#' @keywords internal
createGroupFromInput <- function(session, input, output, dataset, id, type,
                                 selected=NULL, expr=NULL, groupNames=NULL) {
    if (type == "Attribute") {
        if (selected == "") return(NULL)
        group <- createGroupByAttribute(selected, dataset)
        group <- cbind(names(group), type, selected, group)
    } else if (type == "Index/Identifier") {
        strRows <- paste(selected, collapse=", ")
        identifiers <- switch(id, "Patients"=getSubjectId(),
                              "Samples"=getSampleId(),
                              "ASevents"=getASevents(),
                              "Genes"=getGenes())
        allRows <- createGroupById(session, selected, identifiers)
        group <- cbind(groupNames, type, strRows, list(allRows))
    } else if (type == "Subset") {
        # Test expression before running
        set <- tryCatch(subset(dataset, eval(parse(text=expr))), error=return)
        # Display error
        if ("simpleError" %in% class(set)) {
            errorAlert(session, title="Issue with subset expression",
                       "Check if column names are correct.", br(),
                       "The following error was raised:",
                       tags$code(set$message), alertId="alert-side",
                       caller="Data grouping")
            return(NULL)
        }
        rows  <- match(rownames(set), rownames(dataset))
        rows  <- rownames(dataset)[rows]
        group <- cbind(groupNames, type, expr, list(rows))
    } else if (type == "Regex") {
        # Subset dataset column using given regular expression
        colData    <- as.character(dataset[[selected]])
        # Test expression before running
        set <- tryCatch(grep(expr, colData), error=return)
        # Show error to the user
        if ("simpleError" %in% class(set)) {
            errorAlert(session, title="Issue with GREP expression",
                       "The following error was raised:", br(),
                       tags$code(set$message), alertId="alert-side",
                       caller="Data grouping")
            return(NULL)
        }
        set     <- rownames(dataset)[set]
        strRows <- sprintf('"%s" in %s', expr, selected)
        group   <- cbind(groupNames, "GREP", strRows, list(set))
    } else if (type == "PreMadeList") {
        group        <- selectPreMadeGroup(getGeneList(), selected, getGenes())
        groupNames   <- attr(group, "title")
        discardedMsg <- attr(group, "discarded")
        selected     <- gsub("|||", " ~ ", selected, fixed=TRUE)
        group        <- cbind(groupNames, type, selected, list(group))

        if (!is.null(discardedMsg)) {
            removeAlert(output, alertId="alert-main")
            warningAlert(session, title="Discarded values", discardedMsg,
                         alertId="alert-main", caller="Data grouping")
        }
    }
    # Name group if empty
    if (group[[1]] == "") group[[1]] <- "Unlabelled group"
    # Standardise rows
    ns <- c("Names", "Subset", "Input", id)
    if (is.matrix(group)) {
        colnames(group) <- ns
    } else {
        names(group) <- ns
    }
    rownames(group) <- NULL

    if (id %in% c("Patients", "Samples")) {
        group <- matchGroupSubjectsAndSamples(id, group)
    } else if (id %in% c("ASevents", "Genes")) {
        group <- matchGroupASeventsAndGenes(id, group, getASevents())
    }
    return(group)
}

#' Split elements into groups based on a given column of a dataset
#'
#' Elements are identified by their respective row name.
#'
#' @param col Character: column name
#' @param dataset Matrix or data frame: dataset
#'
#' @family functions for data grouping
#' @return Named list with each unique value from a given column and respective
#' elements
#' @export
#'
#' @examples
#' df <- data.frame(gender=c("male", "female"),
#'                  stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))
#' rownames(df) <- paste0("subject-", LETTERS[1:8])
#' createGroupByAttribute(col="stage", dataset=df)
createGroupByAttribute <- function(col, dataset) {
    if (!col %in% colnames(dataset))
        stop("The given attribute was not found in the column names of the",
             "dataset")

    colData <- as.character(dataset[[col]])
    names(colData) <- rownames(dataset)

    # Replace missing values for NA so they are included by which()
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
#'
#' @return Character: values based on given row indexes or identifiers
#' @keywords internal
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
    # Check indexes higher than the number of subjects available
    valid <- parsed <= length(identifiers)

    # Warn about invalid input
    invalid <- union(rows[!matched][!parsable], parsed[!valid])
    if (length(invalid) > 0) {
        discarded <- paste(invalid, collapse=", ")
        warningAlert(session, title="Discarded values", sprintf(
            "The following %s indexes or identifiers were discarded:",
            length(invalid)), tags$code(discarded),
            caller="Data grouping", alertId="alert-main")
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
#' @keywords internal
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
#' performed
#' @param groupName Character: group name (automatically created if \code{NULL}
#' or \code{""})
#' @param first Character: identifiers of the first element (required when
#' performing the \code{complement} operation)
#' @param second Character: identifiers of the second element (required when
#' performing the \code{complement} operation)
#' @param matches Character: match between samples (as names) and subjects (as
#' values)
#' @param type Character: type of group where set operations are to be performed
#' @param assignColoursToGroups Boolean: assign colours to new groups?
#'
#' @return Matrix containing groups (new group is in the first row)
#' @keywords internal
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
        setOperated <- lapply(seq(3), setOperationString)
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
#' @inheritParams appServer
#' @param buttonId Character: ID of the button to trigger operation
#' @param sharedData Shiny app's global variable
#'
#' @inheritParams setOperation
#'
#' @inherit psichomics return
#' @keywords internal
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
            first   <- getSubjectId()
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
#' @return Matrix with groups ordered (or \code{NULL} if there are no groups)
#' @keywords internal
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

        # Show number of subjects or AS events for each group (if available)
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

#' Check type of groups within file
#'
#' @param file Character: file path
#'
#' @return Type of group: \code{Samples}, \code{ASevents} or \code{NULL}
#' @keywords internal
checkGroupType <- function(file) {
    head <- fread(file, header=FALSE, nrows=1, stringsAsFactors=FALSE,
                  data.table=FALSE)

    prepareGroupFileFormat <- function(check) {
        formatSamples <- list(check=check, colNames=1, ignoreRows=1,
                              checkIndex=1, rowCheck=TRUE)
    }
    headerSamples   <- c("Group", paste(c("Sample", "Patient"), "ID"))
    formatSamples   <- prepareGroupFileFormat(headerSamples)
    isSamplesGroup  <- checkFileFormat(formatSamples, head)

    headerASevents  <- c("Group", paste(c("AS event", "Gene"), "ID"))
    formatASevents  <- prepareGroupFileFormat(headerASevents)
    isASeventsGroup <- checkFileFormat(formatASevents, head)

    if (isSamplesGroup) {
        type <- "Samples"
    } else if (isASeventsGroup) {
        type <- "ASevents"
    } else {
        type <- NA
    }
    return(type)
}

#' Export groups to a file
#'
#' @param groups Matrix with groups
#' @param file Character: path to output file
#' @param match Match between elements within groups
#'
#' @return Saves groups to file
#' @keywords internal
exportGroupsToFile <- function(groups, file, match=NULL) {
    if (nrow(groups) == 0) stop("Groups must have one or more rows.")

    if ("Samples" %in% colnames(groups)) {
        uniqueID  <- "Samples"
        matchedID <- "Patients"
    } else if ("Genes" %in% colnames(groups)) {
        uniqueID  <- "ASevents"
        matchedID <- "Genes"
    } else {
        stop("Unexpected group format.")
    }
    data <- NULL
    for (g in seq(nrow(groups))) {
        group        <- groups[g, ]
        name         <- group[["Names"]]
        uniqueGroup  <- group[[uniqueID]]
        matchedGroup <- group[[matchedID]]
        colour       <- group[["Colour"]]

        if (!is.null(uniqueGroup) && length(uniqueGroup) > 0) {
            if (!is.null(matchedGroup) && length(matchedGroup) > 0) {
                # Match unique and matched elements available in the group
                matchedElems <- match[names(match) %in% uniqueGroup]
                matchedElems[!matchedElems %in% matchedGroup] <- NA
                data <- rbind(
                    data,
                    cbind(name, names(matchedElems), matchedElems, colour))

                # Get remaining matched elements
                matchedGroup <- matchedGroup[!matchedGroup %in% matchedElems]
                if (length(matchedGroup) > 0) {
                    data <- rbind(data, cbind(name, NA, matchedGroup, colour))
                }

                # Get remaining unique elements
                uniqueGroup <- uniqueGroup[
                   !uniqueGroup %in% names(matchedElems)]
                if (length(uniqueGroup) > 0) {
                    data <- rbind(data, cbind(name, uniqueGroup, NA, colour))
                }
            } else {
                data <- rbind(data, cbind(name, uniqueGroup, NA, colour))
            }
        } else {
            if (length(matchedGroup) == 0) matchedGroup <- NA
            data <- rbind(data, cbind(name, NA, matchedGroup, colour))
        }
    }
    uniqueID <- gsub("ASevent", "AS event", uniqueID)
    strs <- paste(gsub("s$", "", c(uniqueID, matchedID)), "ID")

    if (ncol(data) == 4) {
        cols <- c("Group", strs, "Colour")
    } else {
        cols <- c("Group", strs)
    }
    colnames(data) <- cols
    rownames(data) <- NULL
    write.table(data, file, quote=FALSE, row.names=FALSE, sep="\t")
}

processCol <- function(col, group) {
    uniqueNonNA <- function(x) unique(x[!is.na(x)])
    res <- lapply(split(col, group), uniqueNonNA)
    res <- res[unique(group)] # Order by original group names
    return(res)
}

discardElemsInList <- function(data, keep=NULL, element=NULL) {
    if (is.null(keep)) return(data)
    filterByMatches <- function(x, keep) return(x[x %in% keep])
    res       <- lapply(data, filterByMatches, keep)
    kept      <- length(unique(unlist(res)))
    total     <- length(unique(unlist(data)))
    discarded <- total - kept
    perc      <- round(discarded / total * 100)
    if (discarded > 0) {
        element <- paste0(tolower(element), "s")
        attr(res, "discarded") <- sprintf("%s out of %s %s (%s%%)",
                                          discarded, total, element, perc)
    }
    return(res)
}

#' Import groups from a file
#'
#' @param file Character: path to file
#' @param uniqueElems Character: vector of unique elements (samples or
#' alternative splicing events)
#' @param matchingElems Character: vector of matching elements (subjects or
#' genes)
#' @param match Match between elements within groups
#'
#' @return Matrix with groups
#' @keywords internal
importGroupsFrom <- function(file, uniqueElems=NULL, matchingElems=NULL,
                             match=NULL, type=NULL) {
    if (is.null(type)) type <- checkGroupType(file)

    if (is.null(type)) {
        stop("File does not contain group data.")
    } else if (type == "Samples") {
        uniqueID  <- "Sample"
        matchedID <- "Patient"
    } else if (type == "ASevents") {
        uniqueID  <- "AS event"
        matchedID <- "Gene"
    }
    strs <- paste(c(uniqueID, matchedID), "ID")

    # Load groups file
    groups <- fread(file, colClasses="character")
    groups[["Group"]][is.na(groups[["Group"]])] <- "NA" # Read NA as group name
    groupNames <- unique(groups[["Group"]])

    data <- matrix(nrow=length(groupNames), ncol=3,
                   dimnames=list(groupNames, c("Names", "Subset", "Input")))
    data[ , "Names"]  <- groupNames
    data[ , "Subset"] <- "Groups loaded from file"
    data[ , "Input"]  <- file

    # Process columns containing samples/AS events
    uniqueCol   <- processCol(groups[[strs[[1]]]], groups[["Group"]])
    # Process columns containing subjects/genes
    matchingCol <- processCol(groups[[strs[[2]]]], groups[["Group"]])

    # Include available matching elements if uniqueCol or matchingCol are empty
    if (!is.null(match)) {
        anyUniqueIdLoaded  <- any(sapply(uniqueCol,   length))
        anyMatchedIdLoaded <- any(sapply(matchingCol, length))
        for (group in groupNames) {
            if (!anyUniqueIdLoaded) {
                uniqueCol[[group]] <- unique(
                    names(match[match %in% matchingCol[[group]]]))
            } else if (!anyMatchedIdLoaded) {
                matchingCol[[group]] <- unique(
                    unname(match[uniqueCol[[group]]]))
            }
        }
    }

    # Discard elements based on whether they are present or not
    uniqueCol   <- discardElemsInList(uniqueCol,   uniqueElems,   uniqueID)
    matchingCol <- discardElemsInList(matchingCol, matchingElems, matchedID)

    df <- cbind(data, matchingCol, uniqueCol)
    colnames(df)[4:5] <- c(paste0(matchedID, "s"), type)

    if (nrow(df) == 0) {
        return(NULL)
    } else if (type == "ASevents") {
        df <- df[ , c(1, 2, 3, 5, 4), drop=FALSE] # Fix column order
    }

    # Include colour if available
    colours <- groups[["Colour"]]
    if (!is.null(colours)) {
        colours <- setNames(colours, groups[["Group"]])
        df <- cbind(df, "Colour"=colours[groupNames])
    }

    # Attach message with number of discarded elements
    uniqueColDiscarded   <- attr(uniqueCol,   "discarded")
    matchingColDiscarded <- attr(matchingCol, "discarded")
    if (!is.null(uniqueColDiscarded) || !is.null(matchingColDiscarded)) {
        msg <- tags$ul(tags$li(
            tags$b(paste(c(matchingColDiscarded, uniqueColDiscarded),
                         collapse=" and ")),
            "do not match loaded data and were discarded"))
        attr(df, "discarded") <- msg
    }
    return(df)
}

#' Logic server to manipulate data grouping
#'
#' @inheritParams appServer
#' @param type Character: type of data for each the interface is intended
#'
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom shinyjs toggleState enable disable hidden show hide
#' @importFrom shiny textInput showTab hideTab
#' @importFrom shinyBS updateCollapse
#' @importFrom colourpicker updateColourInput
#'
#' @return HTML elements
#' @keywords internal
groupManipulation <- function(input, output, session, type) {
    ns <- session$ns

    # Update attributes available for data grouping
    if (type == "Samples") {
        updateAttributes <- function(id) {
            if (id == "Patients") {
                attrs <- getSubjectAttributes()
            } else if (id == "Samples") {
                attrs <- getSampleAttributes()
            }

            suggestedCols <- attr(attrs, "default")
            if (!is.null(suggestedCols)) {
                suggestedIndex <- match(suggestedCols, attrs)
                suggestedIndex <- suggestedIndex[!is.na(suggestedIndex)]
                cols <- list("Suggested attributes"=attrs[suggestedIndex],
                             "Other attributes"=attrs[-suggestedIndex])
            } else {
                cols <- attrs
            }

            if (is.null(cols)) cols <- list()
            updateSelectizeInput(session, paste0("groupAttribute", id),
                                 choices=cols, selected=character())
            updateSelectizeInput(session, paste0("grepColumn", id),
                                 choices=cols, selected=character())
        }

        observe(updateAttributes("Samples"))
        observe(updateAttributes("Patients"))
    }

    # Update identifiers to suggest in index/identifiers panel
    observe({
        convertNull2List <- function(x) {
            if (is.null(x)) x <- list()
            return(x)
        }
        toggleTabsBy <- function(data, dataTab, noDataTab, inputId="groupBy",
                                 select=FALSE) {
            if (is.null(data)) {
                showTab(inputId, noDataTab)
                hideTab(inputId, dataTab)
            } else {
                hideTab(inputId, noDataTab)
                showTab(inputId, dataTab, select=select)
            }
        }

        genes    <- getGenes()
        ASevents <- getASevents()
        if (type == "Samples") {
            samples <- getSampleId()
            toggleTabsBy(samples, "By samples", "NoSamples", select=TRUE)
            updateSelectizeInput(
                session, paste0("groupRows", "Samples"),
                choices=convertNull2List(samples), server=TRUE)

            subjects <- getSubjectId()
            toggleTabsBy(subjects, "By subjects", "NoSubjects", select=TRUE)
            updateSelectizeInput(
                session, paste0("groupRows", "Patients"),
                choices=convertNull2List(subjects), server=TRUE)
        } else if (type == "ASevents") {
            toggleTabsBy(genes, "By genes", "NoGenes", select=TRUE)
            updateSelectizeInput(
                session, paste0("groupRows", "Genes"),
                choices=convertNull2List(genes), server=TRUE)

            toggleTabsBy(ASevents, "By splicing events", "NoASevents",
                         select=TRUE)
            updateSelectizeInput(
                session, paste0("groupRows", "ASevents"),
                choices=convertNull2List(ASevents), server=TRUE)
        }
    })

    createGroupOptions <- function(id, hasAttributes=TRUE) {
        clearSelection <- function(session, element) {
            updateSelectizeInput(session, element, selected=character())
        }

        # Group based on index or identifiers ----------------------------------
        observeEvent(input[[paste0("createGroupRows", id)]], {
            isolate({
                selected   <- input[[paste0("groupRows", id)]]
                groupNames <- input[[paste0("groupNameRows", id)]]
            })
            createGroup(session, input, output, id, type="Index/Identifier",
                        selected=selected, groupNames=groupNames)
            clearSelection(session, paste0("groupRows", id))
            clearSelection(session, paste0("groupNameRows", id))
        })

        # Remaining code is unneeded unless attributes are available
        if (!hasAttributes) return(NULL)

        # Group by attribute ---------------------------------------------------
        observeEvent(input[[paste0("createGroupAttribute", id)]], {
            selected <- isolate(input[[paste0("groupAttribute", id)]])
            createGroup(session, input, output, id, type="Attribute",
                        selected=selected)
            clearSelection(session, paste0("groupAttribute", id))
        })

        # Pre-made list of genes -----------------------------------------------
        observeEvent(input[[paste0("loadPreMadeGroup", id)]], {
            selected <- isolate(input$groupAttributeGenes)
            createGroup(session, input, output, id, type="PreMadeList",
                        selected=selected)
        })

        output[[paste0("geneListNumber", id)]] <- renderUI({
            selected <- input$groupAttributeGenes
            group    <- selectPreMadeGroup(getGeneList(), selected)
            return(tagList(tags$b("Genes:"), length(group)))
        })

        output[[paste0("geneListSource", id)]] <- renderText({
            selected <- input$groupAttributeGenes
            group    <- selectPreMadeGroup(getGeneList(), selected)
            return(attr(group, "citation"))
        })

        # Group subset ---------------------------------------------------------
        observeEvent(input[[paste0("createGroupSubset", id)]], {
            isolate({
                expr       <- input[[paste0("groupExpression", id)]]
                groupNames <- input[[paste0("groupNameSubset", id)]]
            })
            createGroup(session, input, output, id, type="Subset",
                        expr=expr, groupNames=groupNames)
            clearSelection(session, paste0("groupExpression", id))
            clearSelection(session, paste0("groupNameSubset", id))
        })

        # Update available attributes to suggest in the subset expression
        output[[paste0("groupExpressionSuggestions", id)]] <- renderUI({
            if (id == "Patients")
                attrs <- getSubjectAttributes()
            else if (id == "Samples")
                attrs <- getSampleAttributes()
            attrs <- sprintf("`%s`", attrs)
            textSuggestions(ns(paste0("groupExpression", id)), attrs)
        })

        # Group based on regular expression ------------------------------------
        observeEvent(input[[paste0("createGroupRegex", id)]], {
            isolate({
                selected   <- input[[paste0("grepColumn", id)]]
                expr       <- input[[paste0("groupExpression", id)]]
                groupNames <- input[[paste0("groupNameSubset", id)]]
            })
            createGroup(session, input, output, id, type="Regex",
                        selected=selected, expr=expr, groupNames=groupNames)
            clearSelection(session, paste0("groupNameRegex", id))
            clearSelection(session, paste0("groupNameRegex", id))
            clearSelection(session, paste0("groupNameRegex", id))
        })

        # Preview groups to be created
        output[[paste0("previewGroups", id)]] <- renderUI({
            col <- input[[paste0("groupAttribute", id)]]
            if (is.null(col) || col == "") return(NULL)

            if (id == "Patients")
                dataset <- getClinicalData()
            else if (id == "Samples")
                dataset <- getSampleInfo()

            group <- createGroupByAttribute(col, dataset)

            # Only preview up to a given number of groups
            groupsToPreview <- 6
            sub             <- head(group, groupsToPreview)
            sub             <- cbind(names(sub), sub)
            colnames(sub)   <- c("Group", id)
            sub             <- matchGroupSubjectsAndSamples(id, sub)

            table <- cbind(sub[ , "Group"],
                           if ("Patients" %in% colnames(sub))
                               as.character(sapply(sub[ , "Patients"], length)),
                           if ("Samples" %in% colnames(sub))
                               as.character(sapply(sub[ , "Samples"], length)))
            colnames(table) <- gsub("Patient", "Subject", colnames(sub))

            totalGroups <- length(group)
            if (totalGroups > groupsToPreview) {
                table <- rbind(table, rep("(...)", ncol(table)))
                title <- sprintf("Preview %s out of %s groups",
                                 groupsToPreview, totalGroups)
            } else {
                title <- sprintf("Preview %s group%s", totalGroups,
                                 ifelse(totalGroups > 1, "s", ""))
            }

            htmlTable <- table2html(table, rownames=FALSE, thead=TRUE,
                                    class="table table-condensed table-striped")
            htmlTable <- gsub("<th>", '<th style="text-align: right">',
                              htmlTable)
            panel <- tags$div(class="panel panel-info",
                              style="z-index: 101; position: relative;",
                              tags$div(class="panel-heading", title), htmlTable)
            return(panel)
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
    callback=JS("renderGroupTable(table);"),
    options=list(pageLength=10, lengthChange=FALSE, scrollX=TRUE,
                 ordering=FALSE, columnDefs = list(list(
                     orderable=FALSE, className='details-control', targets=0)),
                 language=list(zeroRecords="No groups available to display")))

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

    output[["saveSelectedGroups-button"]] <- downloadHandler(
        filename = function() {
            paste0("psichomics groups ", getCategory(), ".txt")
        }, content = function(file) {
            # Get selected groups
            selected <- sort(input$groupsTable_rows_selected)
            groups   <- getGroups(type, complete=TRUE)
            groups   <- groups[selected, , drop=FALSE]

            if (type == "Samples") {
                match <- getClinicalMatchFrom("Inclusion levels")
            } else if (type == "ASevents" && !is.null(getASevents()) ) {
                match <- getGenesFromSplicingEvents(getASevents())
            } else {
                match <- NULL
            }
            exportGroupsToFile(groups, file, match)
        }
    )

    output[["saveAllGroups-button"]] <- downloadHandler(
        filename = function() {
            paste0("psichomics groups ", getCategory(), ".txt")
        }, content = function(file) {
            groups <- getGroups(type, complete=TRUE)
            if (type == "Samples") {
                match <- getClinicalMatchFrom("Inclusion levels")
            } else if (type == "ASevents" && !is.null(getASevents()) ) {
                match <- getGenesFromSplicingEvents(getASevents())
            } else {
                match <- NULL
            }
            exportGroupsToFile(groups, file, match)
        }
    )

    # Import groups from file
    observeEvent(input[["loadGroups-button"]], {
        groupsFile <- fileBrowser()
        if (is.na(groupsFile)) return(NULL) # Action cancelled by the user

        groupType <- checkGroupType(groupsFile)
        isolate({
            if (groupType == "Samples") {
                uniqueElems   <- getSampleId()
                matchingElems <- getSubjectId()
                match         <- getClinicalMatchFrom("Inclusion levels")
                groupTypeMsg  <- "samples/subjects"
            } else if (groupType == "ASevents") {
                uniqueElems   <- getASevents()
                matchingElems <- getGenes()
                match         <- getGenesFromSplicingEvents(uniqueElems)
                groupTypeMsg  <- "AS events/genes"
            }
        })

        removeAlert(output, alertId="alert-main")
        imported <- tryCatch(
            importGroupsFrom(groupsFile, uniqueElems, matchingElems, match,
                             type=groupType),
            error=return)
        if (!is.null(imported) && !is(imported, "error")) {
            appendNewGroups(groupType, imported)
            groups <- nrow(imported)
            title  <- paste(groups, ifelse(groups == 1, "group", "groups"),
                            sprintf("loaded (containing %s)", groupTypeMsg))

            msg <- attr(imported, "discarded")
            if (!is.null(msg)) {
                warningAlert(
                    session, title=title, msg, br(),
                    "Groups imported based on file:", tags$code(groupsFile),
                    alertId="alert-main", caller="Data grouping")
            } else {
                successAlert(
                    session, title=title,
                    "Groups imported based on file:", tags$code(groupsFile),
                    alertId="alert-main", caller="Data grouping")
            }
        } else {
            errorAlert(session, title="File not loaded", imported$message,
                       alertId="alert-main", caller="Data grouping")
        }
    })

    if (type == "Samples") {
        # Test group indepedence for subject/sample groups only
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

                # Prepare contingency table
                ref <- as.character(point[["Reference"]])
                ref <- gsub("vs others", "", ref, fixed=TRUE)
                contingencyTable <- point[["Contingency table"]][[1]]
                rownames(contingencyTable) <- c(ref, "Others")
                contingencyTable <- table2html(contingencyTable)

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
                             tags$b("Contingency table"), contingencyTable))
            }
        })
    }
}

#' @rdname appServer
groupsServer <- function(input, output, session) {
    callModule(groupManipulation, "sampleGroupModule",  "Samples")
    callModule(groupManipulation, "ASeventGroupModule", "ASevents")
}

#' Server function for data grouping (one call)
#'
#' These functions only run once instead of running for every instance of groups
#'
#' @inheritParams appServer
#'
#' @inherit psichomics return
#' @keywords internal
groupsServerOnce <- function(input, output, session) {
    # Update groups according to the availability of sample identifiers
    observe({
        group <- getGroups("Samples", complete=TRUE)
        # Ignore if there are no groups
        if (is.null(group)) return(NULL)

        samples     <- getSampleId()
        subjects    <- getSubjectId()
        match       <- getClinicalMatchFrom("Inclusion levels")
        showSamples <- "Samples" %in% colnames(group)

        if ( is.null(samples) && showSamples ) {
            # Remove sample identifiers from groups
            if ( "Patients" %in% colnames(group) ) {
                # Groups are made with subjects and samples; only remove samples
                group <- group[ , -match("Samples", colnames(group))]
            } else {
                # Groups are only made of samples; remove all groups
                group <- NULL
            }
            setGroups("Samples", group)
        } else if ( !is.null(samples) && !showSamples && !is.null(subjects) &&
                    !is.null(match) && "Patients" %in% colnames(group)) {
            # Update groups if previously made with subjects only
            samples <- getSampleFromSubject(group[ , "Patients"], samples,
                                            subjects, match=match)
            group <- cbind(group, "Samples"=samples)
            setGroups("Samples", group)
        }
    })

    # Create groups based on pre-made list of genes when loading gene or
    # splicing data
    observe({
        geneExp <- getGeneExpression()
        psi     <- getInclusionLevels()

        if (!is.null(geneExp) || !is.null(psi)) {
            # Intersect gene list with available genes
            genes    <- getGenes()
            geneList <- getGeneList(genes)
            if (is.null(geneList)) return(NULL)

            groups     <- unlist(geneList, recursive=FALSE)
            groupNames <- unlist(lapply(names(geneList), function(i)
                sprintf("%s (%s)", names(geneList[[i]]), i)))
            selected   <- unlist(lapply(names(geneList), function(i)
                paste(i, names(geneList[[i]]), sep=" ~ ")))
            groups     <- cbind(groupNames, "PreMadeList", selected, groups)

            # Standardise rows
            ns <- c("Names", "Subset", "Input", "Genes")
            colnames(groups) <- ns
            rownames(groups) <- NULL
            groups <- matchGroupASeventsAndGenes("Genes", groups, getASevents())

            if (!is.null(groups))
                isolate( appendNewGroups("Genes", groups, clearOld=TRUE) )
        }
    }, label="groups_premadeListOfGenes")

    # Create groups by sample types when loading TCGA data
    observe({
        sampleInfo <- getSampleInfo()
        if (!is.null(sampleInfo) && any(grepl("^TCGA", rownames(sampleInfo)))) {
            new    <- createGroupByAttribute("Sample types", sampleInfo)
            groups <- cbind("Names"=names(new), "Subset"="Attribute",
                            "Input"="Sample types", "Samples"=new)

            # Match samples with subjects (if loaded)
            subjects <- isolate(getSubjectId())
            if (!is.null(subjects)) {
                indiv <- lapply(new, function(i)
                    unname(getSubjectFromSample(i, patientId=subjects)))
                groups <- cbind(groups[ , seq(3), drop=FALSE], "Patients"=indiv,
                                groups[ ,      4, drop=FALSE])
            }

            if (!is.null(groups))
                isolate( appendNewGroups("Samples", groups, clearOld=TRUE) )
        }
    })
}

#' @rdname selectGroupsUI
#'
#' @inheritParams getGroups
#' @inheritParams appServer
#' @param filter Character: get groups only if they are present in this argument
#' (if TCGA-styled gene symbols, they will be "converted" to gene symbols alone)
#'
#' @return \code{getSelectedGroups}: List with selected groups (or \code{NULL}
#' when no groups are selected)
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
#' @seealso \code{\link{testGroupIndependence}()} and
#' \code{\link{plotGroupIndependence}()}
#'
#' @return List of lists containing values based on rownames of \code{df}
#' @export
#' @examples
#' df <- data.frame("race"=c("caucasian", "caucasian", "asian"),
#'                  "gender"=c("male", "female", "male"))
#' rownames(df) <- paste("subject", 1:3)
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
#' @param elements Character: all subject identifiers
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
#' \item{table}{Contingency table used for testing}
#' \item{pvalue}{Fisher's exact test's p-value}
#'
#' @keywords internal
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
#' Test multiple contingency tables comprised by two groups (one reference group
#' and another containing remaining elements) and provided groups.
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
#' @family functions for data grouping
#' @return \code{multiGroupIndependenceTest} object, a data frame containing:
#' \item{attribute}{Name of the original groups compared against the reference
#' groups}
#' \item{table}{Contingency table used for testing}
#' \item{pvalue}{Fisher's exact test's p-value}
#'
#' @seealso \code{\link{parseCategoricalGroups}()} and
#' \code{\link{plotGroupIndependence}()}
#'
#' @export
#' @examples
#' elements <- paste("subjects", 1:10)
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
                          "Contingency table"=I(mat))
        df <- rbind(df, new)
    }
    colnames(df) <- c("Reference", "Attributes", "p-value", "Adjusted p-value",
                      "Contingency table")
    class(df) <- c(class(df), "multiGroupIndependenceTest")
    return(df)
}

#' Discard grouped samples if not within a sample vector
#'
#' @param groups Named list of samples
#' @param samples Character: vector with all available samples
#' @param clean Boolean: clean results?
#'
#' @importFrom stats na.omit
#'
#' @return Groups without samples not found in \code{samples}
#' @keywords internal
discardOutsideSamplesFromGroups <- function(groups, samples, clean=FALSE) {
    filterSamplesInGroups <- function(i) ifelse(i %in% samples, i, NA)
    g <- lapply(groups, filterSamplesInGroups)
    g <- lapply(g, na.omit)
    if (clean) g <- lapply(g, as.character)
    attr(g, "Colour") <- attr(groups, "Colour")[names(g)]
    return(g)
}

#' Plot \code{-log10(p-values)} of the results obtained after multiple group
#' independence testing
#'
#' @param groups \code{multiGroupIndependenceTest} object (obtained after
#' running \code{\link{testGroupIndependence}()})
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
#' @seealso \code{\link{parseCategoricalGroups}()} and
#' \code{\link{testGroupIndependence}()}
#'
#' @family functions for data grouping
#' @return \code{ggplot} object
#' @export
#'
#' @examples
#' elements <- paste("subjects", 1:50)
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
