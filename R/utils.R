badge <- function(inputId, label)
    span(class = "badge", id = inputId, label)

buttonGroups <- function(...) {
    div(class="btn-group", role="group", ...)
}

progressbar <- function(value, min = 0, max = 100, label = paste0(value, "%"),
                        striped = T) {
    stripedText <- ifelse(striped, "progress-bar-striped", "")
    div(class="progress",
        div(class=paste("progress-bar", stripedText),
            role="progressbar",
            "aria-valuenow"=value,
            "aria-valuemin"=min,
            "aria-valuemax"=max,
            style=paste0("width: ", value, "%;"),
            label
        )
    )
}

alertNew <- function(...)
    div(class="alert alert-warning alert-dismissible", role="alert", ...)

dropdown <- function(inputId) {
    div(class="dropdown",
        tags$button(class="btn btn-default dropdown-toggle",
                    type="button", id=inputId, "data-toggle"="dropdown",
                    "aria-haspopup"="true", "aria-expanded"="true",
                    list("Dropdown", span(class="caret"))),
        tags$ul(class="dropdown-menu", "aria-labelledby"="dropdownMenu1",
                tags$li(a(href="#", "Action")),
                tags$li(a(href="#", "Another action")),
                tags$li(role="separator", class="divider"),
                tags$li(a(href="#", "Something else")))
    )
}

#' @export
bsModal2 <- function (id, title, trigger, ..., style, size = NULL, footer = NULL)  {
    if (!missing(style)) {
        modalHeader <- paste("modal-header", style)
    } else {
        modalHeader <- "modal-header"
    }
    if (!is.null(size)) {
        if (size == "large") size = "modal-lg"
        else if (size == "small") size = "modal-sm"
        size <- paste("modal-dialog", size)
    }
    else size <- "modal-dialog"
    bsTag <- shiny::tags$div(
        class = "modal sbs-modal fade", 
        id = id, tabindex = "-1", `data-sbs-trigger` = trigger, 
        shiny::tags$div(
            class = size, shiny::tags$div(
                class = "modal-content", 
                shiny::tags$div(
                    class = modalHeader,
                    shiny::tags$button(
                        type = "button",
                        class = "close", `data-dismiss` = "modal",
                        shiny::tags$span(shiny::HTML("&times;"))), 
                    shiny::tags$h4(class = "modal-title", title)), 
                shiny::tags$div(class = "modal-body", list(...)), 
                shiny::tags$div(
                    class = "modal-footer",
                    shiny::tags$button(type = "button", 
                                       class = "btn btn-default",
                                       `data-dismiss` = "modal", 
                                       "Close"),
                    footer))))
    htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

#' Create a progress object
#' @export
startProgress <- function(message, divisions, global = sharedData) {
    print(message)
    global$progress.divisions <- divisions
    global$progress <- shiny::Progress$new()
    global$progress$set(message = message, value = 0)
}

#' Update a progress object
#' @export
updateProgress <- function(message = "Hang in there", value = NULL, max = NULL,
                           detail = NULL, divisions = NULL, 
                           global = sharedData) {
    if (!is.null(divisions)) {
        startProgress(message, divisions, global)
        return(NULL)
    }
    divisions <- global$progress.divisions
    if (is.null(value)) {
        value <- global$progress$getValue()
        max <- global$progress$getMax()
        value <- value + (max - value)
    }
    amount <- ifelse(is.null(max), value/divisions, 1/max/divisions)
    global$progress$inc(amount = amount, message = message, detail = detail)
    print(paste(message, detail))
}

#' Close the progress even if there's an error
#' @export
closeProgress <- function(message=NULL, global = sharedData) {
    # Close the progress even if there's an error
    if (!is.null(message)) print(message)
    global$progress$close()
}

#' Allows to add id to an image
#' @export
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
        htmltools::htmlDependencies(iconTag) <- htmltools::htmlDependency(
            "font-awesome", "4.5.0", c(href = "shared/font-awesome"), 
            stylesheet = "css/font-awesome.min.css")
    }
    iconTag
}

#' Disable a tab from the navbar
#' @export
disableTab <- function(tab) {
    # Style item as disabled
    addClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
             class = "disabled")
    # Disable link itself
    disable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Enable a tab from the navbar
#' @export
enableTab <- function(tab) {
    # Style item as enabled
    removeClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
                class = "disabled")
    # Enable link itself
    enable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Add scatter
#' @export
hc_scatter <- function (hc, x, y, z = NULL, color = NULL, label = NULL, 
                        showInLegend = FALSE, viridis.option = "D", ...) {
    assertthat::assert_that(highcharter:::.is_highchart(hc), 
                            length(x) == length(y), 
                            is.numeric(x), is.numeric(y))
    df <- data_frame(x, y)
    if (!is.null(z)) {
        assert_that(length(x) == length(z))
        df <- df %>% mutate(z = z)
    }
    if (!is.null(color)) {
        assert_that(length(x) == length(color))
        assert_that(viridis.option %in% c("A", "B", "C", "D"))
        cols <- colorize_vector(color, option = viridis.option)
        df <- df %>% mutate(valuecolor = color, color = cols)
    }
    if (!is.null(label)) {
        assert_that(length(x) == length(label))
        df <- df %>% mutate(label = label)
    }
    # Add arguments to data points if they match the length of the data
    args <- list(...)
    for (i in seq_along(args)) {
        if (length(x) == length(args[[i]])) {
            df <- cbind(df, setNames(list(args[i]), names(args)[i]))
            args[[i]] <- character(0)
        }
    }
    
    ds <- list.parse3(df)
    type <- ifelse(!is.null(z), "bubble", "scatter")
    if (!is.null(label)) {
        dlopts <- list(enabled = TRUE, format = "{point.label}")
    }
    else {
        dlopts <- list(enabled = FALSE)
    }
    args <- Filter(length, args)
    do.call("hc_add_series", c(list(hc, data = ds, type = type, 
                                    showInLegend = showInLegend, 
                                    dataLabels = dlopts), args))
}

#' Create script for autocompletion for text input
#' 
#' Uses the JavaScript library jquery.textcomplete
#' 
#' @param id Character: input ID
#' @param words Character: words to suggest
#' @param novalue Character: string when there's no matching values
#'
#' @return HTML string with the JavaScript script prepared to run
#' @export
#' 
#' @examples 
#' words <- c("tumor_stage", "age", "gender")
#' textComplete("textareaid", words)
textComplete <- function(id, words, novalue = "No matching value") {
    words <- paste0('["', paste(words, collapse = '", "'), '"]')
    js <- paste0('$("#', id, '").textcomplete([{
        match: /([a-zA-Z0-9_\\.]{1,})$/,
        search: function(term, callback) {
            term = term.toLowerCase();
            callback($.map(', words, ',
                function(word) {
                    word = word.toLowerCase();
                    var comp = term.length < 4 ? 1 : term.length * 2;
                    var match = fuzzy(word, term);
                    return match.score >= comp ? word : null;
                }
            ));
        },
        index: 1,
        replace: function(word) {
            return word + " ";
        }}], { noResultsMessage: "', novalue, '"});')
    js <- paste("<script>", js, "</script>")
    return(HTML(js))
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
#' @examples
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'
#' ui <- fluidPage(
#'   textAreaInput("caption", "Caption", "Data Summary", width = "1000px"),
#'   verbatimTextOutput("value")
#' )
#' server <- function(input, output) {
#'   output$value <- renderText({ input$caption })
#' }
#' shinyApp(ui, server)
#' }
#' @export
textAreaInput <- function(inputId, label, value = "", width = NULL,
                          placeholder = NULL) {
    
    div(class = "form-group shiny-input-container",
        style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"),
        tags$label(label, `for` = inputId),
        tags$textarea(id = inputId, class = "form-control", value = value,
                      placeholder = placeholder)
    )
}