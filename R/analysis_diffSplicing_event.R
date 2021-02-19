#' @rdname appUI
#'
#' @importFrom highcharter highchartOutput
#' @importFrom shiny tagList uiOutput NS sidebarLayout numericInput h3 mainPanel
#' actionButton sidebarPanel
diffSplicingEventUI <- function(id) {
    ns <- NS(id)
    return(diffEventUI(id, ns, psi=TRUE))
}

#' @rdname appServer
#'
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide
diffSplicingEventServer <- function(input, output, session) {
    ns <- session$ns

    selectGroupsServer(session, "diffGroups", "Samples")

    observe({
        psi <- getInclusionLevels()
        event <- getEvent()
        if (is.null(psi) || is.null(event) || event == "") return(NULL)
        data <- as.numeric(psi[event, ])
        updateCheckboxInput(session, "rug", value=length(data) < 500)
    })

    observeEvent(input$analyse, {
        # Get splicing event's inclusion levels
        psi <- getInclusionLevels()
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        }

        # Get selected event
        event <- getEvent()
        if (is.null(event) || event == "") {
            errorModal(session, "No event selected",
                       "Please, select an alternative splicing event.",
                       caller="Differential splicing analysis")
            return(NULL)
        }

        # Prepare groups of samples to analyse
        groups <- getSelectedGroups(input, "diffGroups", "Samples",
                                    filter=colnames(psi))
        colour <- attr(groups, "Colour")
        if ( !is.null(groups) ) {
            attrGroups <- groups
            psi <- psi[ , unlist(groups), drop=FALSE]
            groups <- rep(names(groups), sapply(groups, length))
        } else {
            attrGroups <- "All samples"
            groups <- rep(attrGroups, ncol(psi))
        }

        # Check if analyses were already performed
        stats <- getDifferentialSplicing()
        if (!is.null(stats) && identical(attrGroups, attr(stats, "groups"))) {
            stat <- stats[event, ]
        } else {
            stat <- NULL
        }

        # Separate samples by their groups
        eventPSI <- as.numeric(psi[event, ])
        names(eventPSI) <- colnames(psi)
        eventPSI <- filterGroups(eventPSI, groups, 2)
        groups   <- attr(eventPSI, "Groups")
        attr(groups, "Colour") <- colour

        title <- parseSplicingEvent(event, char=TRUE, pretty=TRUE)
        plot  <- plotDistribution(eventPSI, groups, title=title,
                                  type=isolate(input$plotType), rug=input$rug)
        output$density <- renderHighchart(plot)

        output$eventDiagrams <- renderUI({
            parsed <- parseSplicingEvent(event)
            if (is.null(parsed$type)) return(NULL)
            isMXE <- parsed$type == "MXE"
            constitutive <- suppressWarnings(
                plotSplicingEvent(
                    style="position: absolute; top: 321px; left: 52px",
                    constitutiveWidth=40, alternativeWidth=40, intronWidth=0,
                    event, class=NULL, showPath=FALSE, showText=FALSE,
                    showAlternative1=FALSE, showAlternative2=TRUE)[[1]])
            alternative <- suppressWarnings(
                plotSplicingEvent(
                    style="position: absolute; top: 321px; right: 25px",
                    constitutiveWidth=40, alternativeWidth=40, intronWidth=0,
                    event, class=NULL, showPath=FALSE, showText=FALSE,
                    showAlternative1=TRUE, showAlternative2=!isMXE)[[1]])
            return(tagList(HTML(constitutive), HTML(alternative)))
        })

        output$basicStats <- renderUI(basicStats(eventPSI, groups))
        output$ttest      <- renderUI(ttest(eventPSI, groups, stat))
        output$wilcox     <- renderUI(wilcox(eventPSI, groups, stat))
        output$kruskal    <- renderUI(kruskal(eventPSI, groups, stat))
        output$levene     <- renderUI(levene(eventPSI, groups, stat))
        output$fligner    <- renderUI(fligner(eventPSI, groups, stat))
        # output$fisher   <- renderUI(fisher(eventPSI, groups))
        # output$spearman <- renderUI(spearman(eventPSI, groups))

        show("survivalButton")
        show("singleEventInfo")
    })

    observeEvent(input$missingInclusionLevels,
                 missingDataGuide("Inclusion levels"))
    observeEvent(input$missingIncLevelsButton,
                 missingDataGuide("Inclusion levels"))

    # Toggle options only if required data is available
    observe({
        # Get splicing event's inclusion levels
        psi <- getInclusionLevels()
        if (is.null(psi)) {
            show("missingIncLevels")
            hide("singleEventOptions")
            hide("survivalButton")
            hide("singleEventInfo")
        } else {
            hide("missingIncLevels")
            show("singleEventOptions")
        }
    })
}

attr(diffSplicingEventUI, "loader") <- "diffSplicing"
attr(diffSplicingEventUI, "name") <- "Individual alternative splicing event"
attr(diffSplicingEventServer, "loader") <- "diffSplicing"
