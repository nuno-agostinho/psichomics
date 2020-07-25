#' Create HTML table from data frame or matrix
#' 
#' @param data Data frame or matrix
#' @param rownames Boolean: print row names?
#' @param colnames Boolean: print column names?
#' @param class Character: table class
#' @param style Character: table style
#' @param thead Boolean: add a \code{thead} tag to the first row?
#' 
#' @importFrom xtable xtable print.xtable
#' @importFrom shiny HTML
#' 
#' @return HTML elements
#' @keywords internal
table2html <- function(data, rownames=TRUE, colnames=TRUE, class=NULL, 
                       style=NULL, thead=FALSE) {
    table <- xtable(data)
    table <- print(table, type="html", print.results=FALSE,
                   include.rownames=rownames, include.colnames=colnames)
    html <- HTML(table)
    
    if (thead)
        html <- gsub("(<tr>[[:space:]]*<th>.*</th>[[:space:]]*</tr>)",
                     "<thead>\\1</thead>", html)
    
    if (!is.null(class))
        class <- sprintf('class="%s"', paste(class, collapse=" "))
    
    if (!is.null(style))
        style <- sprintf('style="%s"', paste(style, collapse=" "))
    
    rep <- paste(class, style)
    if (length(rep) > 0) 
        html <- gsub("border=1", rep, html, fixed=TRUE)
    
    return(html)
}

#' Interface for interactive \code{\link[ggplot2]{ggplot}}
#' 
#' @param id Character: identifier
#' 
#' @importFrom shiny tagList plotOutput brushOpts hoverOpts uiOutput 
#' actionButton
#' @importFrom shinyjs hidden
#' 
#' @return HTML elements
#' @keywords internal
ggplotUI <- function(id) {
    idd <- function(str) paste(id, str, sep="-")
    plotId    <- idd("plot")
    tooltipId <- idd("tooltip")
    brushId   <- idd("brush")
    hoverId   <- idd("hover")
    resetId   <- idd("resetZoom")
    tagList(
        # Mimic Highcharts button to reset zoom level
        hidden(actionButton(
            resetId, "Reset zoom",
            style="font-size: 13px;",
            style="background-color: #f7f7f7;",
            style="border-color: #cccccc;",
            style="padding-bottom: 5px;", style="padding-top: 5px;",
            style="padding-left: 9px;", style="padding-right: 9px;",
            style="position: absolute;", style="z-index: 1;",
            style="top: 20px;", style="right: 35px;")),
        uiOutput(tooltipId),
        plotOutput(plotId,
                   brush=brushOpts(brushId, resetOnNew=TRUE),
                   hover=hoverOpts(hoverId, delay=50, delayType="throttle")))
}

#' Create the interface for the tooltip of a plot
#' 
#' @param df Data frame
#' @param hover Mouse hover information for a given plot as retrieved from
#' \code{\link[shiny]{hoverOpts}}
#' @param x Character: name of the variable used for the X axis
#' @param y Character: name of the variable used for the Y axis
#' @param eventData Alternative splicing event information (if available)
#' 
#' @importFrom shiny tags nearPoints wellPanel
#' 
#' @return HTML elements
#' @keywords internal
ggplotTooltip <- function(df, hover, x, y, eventData=NULL) {
    point <- nearPoints(df, hover, threshold=10, maxpoints=1, addDist=TRUE,
                        xvar=x, yvar=y)
    if (nrow(point) == 0) return(NULL)
    
    trItem <- function(key, value) tags$tr(tags$td(tags$b(key)), tags$td(value))
    
    thisPoint <- rownames(point)
    if ( areSplicingEvents(thisPoint, data=eventData) ) {
        res  <- prepareEventInfoTooltip(thisPoint, data=eventData)
        head <- tags$div(style="padding: 5px;", tags$small(tags$b(thisPoint)))
        if (!is.null(res)) {
            gene    <- trItem("Gene", res$gene[[1]])
            type    <- trItem("Event type", res$subtype[[1]])
            coord   <- trItem("Coordinates", res$coord[[1]])
            diagram <- plotSplicingEventHelper(thisPoint, data=eventData)
            if (!is.null(diagram)) diagram <- trItem("Diagram", diagram)
            params  <- tagList(gene, type, coord, diagram)
        } else {
            params  <- NULL
        }
    } else {
        head   <- tags$thead(trItem("Gene", thisPoint))
        params <- NULL
    }
    
    left <- hover$coords_css$x + 10
    top  <- hover$coords_css$y + 5
    tooltipStyle <- paste0("position: absolute; z-index: 100;",
                           "background-color: rgba(245, 245, 245, 0.85);",
                           "right: calc(100%% - %spx); top: %spx;")
    tooltipStyle <- sprintf(tooltipStyle, left, top)
    tooltip <- wellPanel(class="well-sm", style=tooltipStyle, tags$table(
        class="table table-condensed", style="margin-bottom: 0;", head,
        tags$tbody(params,
                   trItem(x, roundDigits(point[[x]])),
                   trItem(y, roundDigits(point[[y]])))))
    return(tooltip)
}

#' Logic set to create an interactive \code{\link[ggplot2]{ggplot}}
#' 
#' @inheritParams appServer
#' @param id Character: identifier
#' @param plot Character: plot expression (if \code{NULL}, no plots are
#' rendered)
#' @inheritParams ggplotTooltip
#' 
#' @importFrom shiny renderPlot renderUI
#' 
#' @inherit psichomics return
#' @keywords internal
ggplotServer <- function(input, output, id, plot=NULL, df=NULL, x=NULL, 
                         y=NULL, eventData=NULL) {
    idd <- function(str) paste(id, str, sep="-")
    output[[idd("plot")]] <- renderPlot(plot)
    
    if (is.null(plot)) {
        output[[idd("tooltip")]] <- renderUI(NULL)
    } else {
        output[[idd("tooltip")]] <- renderUI(
            ggplotTooltip(df, input[[idd("hover")]], x, y, eventData))
    }
}

#' @rdname ggplotServer
#' 
#' @note Insert \code{ggplotAuxSet} outside any observer (so it is only run 
#' once)
ggplotAuxServer <- function(input, output, id) {
    idd <- function(str) paste(id, str, sep="-")
    
    # Save zoom coordinates according to brushed area of the plot
    observe({
        brush <- input[[idd("brush")]]
        if (!is.null(brush)) {
            setZoom(id, brush)
            setSelectedPoints(id, NULL)
        }
    })
    
    # Toggle visibility of reset zoom button
    observe({
        zoom <- getZoom(id)
        if (is.null(zoom)) {
            hide(idd("resetZoom"))
        } else {
            show(idd("resetZoom"))
        }
    })
    
    # Reset zoom when clicking the respective button
    observeEvent(input[[idd("resetZoom")]], {
        setZoom(id, NULL)
        setSelectedPoints(id, NULL)
    })
}

inheritAttrs <- function(original, objectToCopyFrom,
                         avoid=c("names", "row.names", "class")) {
    notNames <- !names(attributes(objectToCopyFrom)) %in% c(
        names(attributes(original)), avoid)
    attributes(original) <- c(attributes(original),
                              attributes(objectToCopyFrom)[notNames])
    colnames(original) <- colnames(objectToCopyFrom)
    return(original)
}
