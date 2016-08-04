.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Start the visual interface by running the function ",
                          "psichomics()")
}

#' Modified version of shinyBS::bsModal
#' 
#' bsModal is used within the UI to create a modal window. This allows to use
#' the footer.
#' 
#' @inheritParams shinyBS::bsModal
#' @param footer UI set: List of elements to include in the footer
#' @param style Character: message style can be "warning", "error", "info" or 
#' NULL
#' @param size Character: Modal size ("small", "default" or "large")
#' 
#' @importFrom shiny tags HTML
#' @importFrom htmltools attachDependencies htmlDependency
#' @importFrom utils packageVersion
bsModal2 <- function (id, title, trigger, ..., size=NULL, footer=NULL, 
                      style = NULL)  {
    if (!is.null(style)) {
        style <- match.arg(style, c("info", "warning", "error"))
        modalHeader <- paste("modal-header", style)
    } else {
        modalHeader <- "modal-header"
    }
    if (!is.null(size)) {
        if (size == "large") size = "modal-lg"
        else if (size == "small") size = "modal-sm"
        else if (size == "default") size = "modal-sm"
        size <- paste("modal-dialog", size)
    }
    else size <- "modal-dialog"
    bsTag <- tags$div(
        class = "modal sbs-modal fade", 
        id = id, tabindex = "-1", `data-sbs-trigger` = trigger, 
        tags$div(
            class = size, tags$div(
                class = "modal-content", 
                tags$div(
                    class = modalHeader,
                    tags$button(
                        type = "button",
                        class = "close", `data-dismiss` = "modal"), 
                    tags$h4(class = "modal-title", title)), 
                tags$div(class = "modal-body", list(...)), 
                tags$div(
                    class = "modal-footer",
                    tags$button(type = "button", 
                                       class = "btn btn-default",
                                       `data-dismiss` = "modal", 
                                       "Close"),
                    footer))))
    shinyBSDep <- htmlDependency("shinyBS", packageVersion("shinyBS"),
                                 src=c(href="sbs"), script="shinyBS.js", 
                                 stylesheet="shinyBS.css")
    attachDependencies(bsTag, shinyBSDep)
}

#' Allows to add id to an image
#' 
#' @param name Character: name of the icon
#' @param class Additional classes to customise the icon style
#' @param lib Icon library to use (either "font-awesome" or "glyphicon")
#' @param ... Extra arguments to the icon tag
#' 
#' @importFrom htmltools htmlDependencies
icon2 <- function (name, class = NULL, lib = "font-awesome", ...) {
    prefixes <- list(`font-awesome` = "fa", glyphicon = "glyphicon")
    prefix <- prefixes[[lib]]
    if (is.null(prefix)) {
        stop("Unknown font library '", lib, "' specified. Must be one of ", 
             paste0("\"", names(prefixes), "\"", collapse = ", "))
    }
    iconClass <- ""
    if (!is.null(name)) 
        iconClass <- paste0(prefix, " ", prefix, "-", name)
    if (!is.null(class)) 
        iconClass <- paste(iconClass, class)
    iconTag <- tags$i(class = iconClass, ...)
    if (lib == "font-awesome") {
        htmlDependencies(iconTag) <- htmlDependency(
            "font-awesome", "4.5.0", c(href = "shared/font-awesome"), 
            stylesheet = "css/font-awesome.min.css")
    }
    iconTag
}

#' Disable a tab from the navbar
#' @importFrom shinyjs disable addClass
#' @param tab Character: tab to disable
disableTab <- function(tab) {
    # Style item as disabled
    addClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
             class = "disabled")
    # Disable link itself
    disable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Enable a tab from the navbar
#' @importFrom shinyjs removeClass enable
#' @param tab Character: tab to enable
enableTab <- function(tab) {
    # Style item as enabled
    removeClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
                class = "disabled")
    # Enable link itself
    enable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Create script for autocompletion of text input
#' 
#' Uses the JavaScript library jquery.textcomplete
#' 
#' @param id Character: input ID
#' @param words Character: words to suggest
#' @param novalue Character: string when there's no matching values
#' @param char Character to succeed accepted word
#'
#' @return HTML string with the JavaScript script prepared to run
#' 
#' @examples 
#' words <- c("tumor_stage", "age", "gender")
#' psichomics:::textSuggestions("textareaid", words)
textSuggestions <- function(id, words, novalue="No matching value", char=" ") {
    varId <- paste0(gsub("-", "_", id), "_words")
    var <- paste0(varId, ' = ["', paste(words, collapse = '", "'), '"];')
    
    js <- paste0('$("#', escape(id), '").textcomplete([{
        match: /([a-zA-Z0-9_\\.]{1,})$/,
        search: function(term, callback) {
            var words = ', varId, ', sorted = [];
            for (i = 0; i < words.length; i++) {
                sorted[i] = fuzzy(words[i], term);
            }
            sorted.sort(fuzzy.matchComparator);
            sorted = sorted.map(function(i) { return i.term; });
            callback(sorted);
        },
        index: 1,
        cache: true,
        replace: function(word) {
            return word + "', char ,'";
        }}], { noResultsMessage: "', novalue, '"});')
    js <- HTML("<script>", var, js, "</script>")
    return(js)
}

#' Create a textarea input control
#'
#' Create a resizable input control for entry of unstructured text values
#'
#' @param inputId The \code{input} slot that will be used to access the value.
#' @param label Display label for the control, or \code{NULL} for no label.
#' @param value Initial value.
#' @param width The width of the input, e.g. \code{'400px'}, or \code{'100\%'};
#'   see \code{\link{validateCssUnit}}.
#' @param placeholder A character string giving the user a hint as to what can
#'   be entered into the control. Internet Explorer 8 and 9 do not support this
#'   option.
#' @return A textarea input control that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link{updateTextAreaInput}}
#'
#' @importFrom shiny div validateCssUnit
#'
#' @examples
#' ## Only run examples in interactive R sessions
#' \dontrun{
#' if (interactive()) {
#'
#' ui <- fluidPage(
#'   psichomics:::textAreaInput("caption", "Caption", "Data Summary",
#'                              width = "1000px"),
#'   verbatimTextOutput("value")
#' )
#' server <- function(input, output) {
#'   output$value <- renderText({ input$caption })
#' }
#' shinyApp(ui, server)
#' }
#' }
textAreaInput <- function(inputId, label, value = "", width = NULL,
                          placeholder = NULL) {
    div(class = "form-group shiny-input-container",
        style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"),
        tags$label(label, `for` = inputId),
        tags$textarea(id = inputId, class = "form-control", placeholder = placeholder,
                      paste(value, collapse = "\n"))
    )
}

#' Change the value of a textarea input on the client
#'
#' @param value The value to set for the input object.
#'
#' @seealso \code{\link{textAreaInput}}
#' @inheritParams shiny::updateTextInput
#'
#' @examples
#' \dontrun{
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'
#' ui <- fluidPage(
#'   sliderInput("controller", "Controller", 0, 20, 10),
#'   textAreaInput("inText", "Input textarea"),
#'   textAreaInput("inText2", "Input textarea 2")
#' )
#'
#' server <- function(input, output, session) {
#'   observe({
#'     # We'll use the input$controller variable multiple times, so save it as x
#'     # for convenience.
#'     x <- input$controller
#'
#'     # This will change the value of input$inText, based on x
#'     updateTextAreaInput(session, "inText", value = paste("New text", x))
#'
#'     # Can also set the label, this time for input$inText2
#'     updateTextAreaInput(session, "inText2",
#'       label = paste("New label", x),
#'       value = paste("New text", x))
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#' }
updateTextAreaInput <- updateTextInput

#' Plot survival curves using Highcharts
#' 
#' @param object A survfit object as returned from the \code{survfit} function
#' @param ... Extra parameters to pass to \code{hc_add_series} function
#' @param fun Name of function or function used to transform the survival curve:
#' \code{log} will put y axis on log scale, \code{event} plots cumulative events
#' (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard function (f(y) =
#' -log(y)), and \code{cloglog} creates a complimentary log-log survival plot
#' (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param markTimes Label curves marked at each censoring time? TRUE by default
#' @param symbol Symbol to use as marker (plus sign by default)
#' @param markerColor Color of the marker ("black" by default); use NULL to use
#' the respective color of each series
#' @param ranges Plot interval ranges? FALSE by default
#' @param rangesOpacity Opacity of the interval ranges (0.3 by default)
#' 
#' @importFrom highcharter %>% hc_add_series highchart hc_tooltip hc_yAxis
#' hc_plotOptions fa_icon_mark JS
#' @importFrom rlist list.parse
#' @return Highcharts object to plot survival curves
#' 
#' @examples
#' 
#' # Plot Kaplan-Meier curves
#' require("survival")
#' require("highcharter")
#' leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
#' hchart(leukemia.surv)
#' 
#' # Plot the cumulative hazard function
#' lsurv2 <- survfit(Surv(time, status) ~ x, aml, type='fleming') 
#' hchart(lsurv2, fun="cumhaz")
#' 
#' # Plot the fit of a Cox proportional hazards regression model
#' fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian)
#' ovarian.surv <- survfit(fit, newdata=data.frame(age=60))
#' hchart(ovarian.surv, ranges = TRUE)
hchart.survfit <- function(object, ..., fun = NULL, markTimes = TRUE,
                           symbol = "plus", markerColor = "black",
                           ranges = FALSE, rangesOpacity = 0.3) {
    groups <- NULL
    # Check if there are groups
    if (is.null(object$strata))
        strata <- c("Series 1" = length(object$time))
    else
        strata <- object$strata
    
    # Modify data according to functions (adapted from survival:::plot.survfit)
    if (is.character(fun)) {
        tfun <- switch(fun,
                       log = function(x) x,
                       event = function(x) 1 - x,
                       cumhaz = function(x) -log(x),
                       cloglog = function(x) log(-log(x)),
                       pct = function(x) x * 100,
                       logpct = function(x) 100 * x,
                       identity = function(x) x,
                       function(x) x)
    } else if (is.function(fun)) {
        tfun <- fun
    } else {
        tfun <- function(x) x
    }
    
    firsty <- tfun(1)
    object$surv <- tfun(object$surv)
    if (ranges && !is.null(object$upper)) {
        object$upper <- tfun(object$upper)
        object$lower <- tfun(object$lower)
    }
    
    # Prepare data
    data <- data.frame(x=object$time, y=object$surv,
                       up=object$upper, low=object$lower,
                       group=rep(names(strata), strata), 
                       stringsAsFactors = FALSE)
    # Data markers
    marker <- list(list(fillColor=markerColor, symbol=symbol, enabled=TRUE))
    if(markTimes)
        mark <- object$n.censor == 1
    else
        mark <- FALSE
    
    # Adjust Y axis range
    yValues <- object$surv
    ymin <- ifelse(min(yValues) >= 0, 0, min(yValues))
    ymax <- ifelse(max(yValues) <= 1, 1, max(yValues))
    
    hc <- highchart() %>%
        hc_tooltip(pointFormat="{point.y}") %>%
        hc_yAxis(min=ymin, max=ymax) %>%
        hc_plotOptions(line = list(marker = list(enabled = FALSE)))
    
    count <- 0
    
    # Process groups by columns (CoxPH-like) or in a single column
    if(!is.null(ncol(object$surv))) {
        groups <- seq(ncol(object$surv))
    } else {
        groups <- names(strata)
    }
    
    summ <- summary(object)$table
    for (name in groups) {
        if (!is.null(ncol(object$surv))) {
            df <- df[c("x", paste(c("y", "low", "up"), col, sep="."))]
            names(df) <- c("x", "y", "low", "up")
            submark <- mark
        } else {
            df <- subset(data, group == name)
            submark <- mark[data$group == name]
        }
        
        # Add first value if there is no value for time at 0 in the data
        if (!0 %in% df$x)
            first <- list(list(x=0, y=firsty))
        else
            first <- NULL
        
        # Mark events
        ls <- list.parse(df)
        names(ls) <- NULL
        if (markTimes)
            ls[submark] <- lapply(ls[submark], c, marker=marker)
        
        if (is.matrix(summ))
            curveSumm <- summ[name, ]
        else
            curveSumm <- summ
        
        hc <- do.call(hc_add_series, c(list(
            hc, data=c(first, ls), step="left", name=name, zIndex=1,
            color=JS("Highcharts.getOptions().colors[", count, "]"), ...),
            curveSumm))
        
        if (ranges && !is.null(object$upper)) {
            # Add interval range
            range <- lapply(ls, function(i) 
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=rangesOpacity, 
                lineWidth=0,
                color=JS("Highcharts.getOptions().colors[", count, "]"),
                ...)
        }
        count <- count + 1
    }
    
    return(hc)
}

#' Render a data table with Sparkline HTML elements
#' 
#' @details This slighlty modified version of \code{\link{renderDataTable}}
#' calls a JavaScript function to convert the Sparkline HTML elements to
#' interactive Highcharts
#' 
#' @importFrom DT renderDataTable JS
#' 
#' @param ... Arguments to pass to \code{\link{renderDataTable}}
#' @param options List of options to pass to \code{\link{renderDataTable}}
renderDataTableSparklines <- function(..., options=NULL) {
    # Escape is set to FALSE to render the Sparkline HTML elements
    renderDataTable(..., escape=FALSE, env=parent.frame(n=1), options=c(
        list(drawCallback=JS("drawSparklines")), options))
}