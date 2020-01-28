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
        plotOutput(plotId,
                   brush=brushOpts(brushId, resetOnNew=TRUE),
                   hover=hoverOpts(hoverId, delay=50, delayType="throttle")),
        uiOutput(tooltipId))
}

#' Create the interface for the tooltip of a plot
#' 
#' @param df Data frame
#' @param hover Mouse hover information for a given plot as retrieved from
#' \code{\link[shiny]{hoverOpts}}
#' @param x Character: name of the variable used for the X axis
#' @param y Character: name of the variable used for the Y axis
#' 
#' @importFrom shiny tags nearPoints wellPanel
#' 
#' @return HTML elements
#' @keywords internal
ggplotTooltip <- function(df, hover, x, y) {
    point <- nearPoints(df, hover, threshold=10, maxpoints=1, addDist=TRUE,
                        xvar=x, yvar=y)
    if (nrow(point) == 0) return(NULL)
    
    # Calculate point position inside the image as percent of total 
    # dimensions from left (horizontal) and from top (vertical)
    xDomain   <- hover$domain$right - hover$domain$left
    right_pct <- (hover$domain$right - hover$x) / xDomain
    yDomain   <- hover$domain$top - hover$domain$bottom
    top_pct   <- (hover$domain$top - hover$y) / yDomain
    
    # Calculate distance from left and bottom in pixels
    xRange   <- hover$range$right - hover$range$left
    right_px <- right_pct * xRange + 25
    yRange   <- hover$range$bottom - hover$range$top
    top_px   <- hover$range$top + top_pct * yRange + 2
    
    trItem <- function(key, value) tags$tr(tags$td(tags$b(key)), tags$td(value))
    
    thisPoint <- rownames(point)
    if ( areSplicingEvents(thisPoint) ) {
        event  <- parseSplicingEvent(thisPoint, pretty=TRUE)
        strand <- ifelse(event$strand == "+", "forward", "reverse")
        gene   <- paste(event$gene[[1]], collapse=" or ")
        type   <- trItem("Event type", event$type)
        coord  <- trItem(
            "Coordinates", 
            sprintf("chr %s: %s to %s (%s strand)", event$chrom,
                    event$pos[[1]][[1]], event$pos[[1]][[2]], strand))
    } else {
        gene  <- thisPoint
        type  <- NULL
        coord <- NULL
    }
    
    # Tooltip
    wellPanel(
        class="well-sm",
        style=paste0("position: absolute; z-index: 100;",
                     "background-color: rgba(245, 245, 245, 0.85); ",
                     "right:", right_px, "px; top:", top_px, "px;"),
        tags$table(class="table table-condensed", style="margin-bottom: 0;",
                   tags$thead( trItem("Gene", gene)),
                   tags$tbody( type, coord,
                               trItem(x, roundDigits(point[[x]])),
                               trItem(y, roundDigits(point[[y]])))))
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
                         y=NULL) {
    idd <- function(str) paste(id, str, sep="-")
    output[[idd("plot")]] <- renderPlot(plot)
    
    if (is.null(plot)) {
        output[[idd("tooltip")]] <- renderUI(NULL)
    } else {
        output[[idd("tooltip")]] <- renderUI(
            ggplotTooltip(df, input[[idd("hover")]], x, y))
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
        if (is.null(zoom))
            hide(idd("resetZoom"))
        else
            show(idd("resetZoom"))
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
