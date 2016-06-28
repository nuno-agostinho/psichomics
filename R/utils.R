#' Sample variance by row
#' 
#' Calculate the sample variance of each row in the given matrix
#' 
#' @param x Matrix
#' @param na.rm Boolean: should the NAs be ignored? FALSE by default
#' 
#' @return Variance for each row
rowVar <- function (x, na.rm = FALSE) {
    means <- rowMeans(x, na.rm = na.rm)
    meansSqDev <- (x - means)^2
    squaresSum <- rowSums(meansSqDev, na.rm = na.rm)
    nas <- rowSums(is.na(x))
    return(squaresSum/(ncol(x) - nas - 1))
}

#' Return the type of a given sample
#' 
#' @param sample Character: ID of the sample
#' @param filename Character: path to RDS file containing corresponding type
getSampleTypes <- function(sample, 
                           filename = system.file("extdata",  
                                                  "TCGAsampleType.RDS",
                                                  package="psichomics")) {
    typeList <- readRDS(filename)
    type <- gsub(".*?-([0-9]{2}).-.*", "\\1", sample, perl = TRUE)
    return(typeList[type])
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
#' @export
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
                        class = "close", `data-dismiss` = "modal",
                        tags$span(HTML("&times;"))), 
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
                                 src=c(href="sbs"),
                                 script="shinyBS.js", 
                                 stylesheet="shinyBS.css")
    attachDependencies(bsTag, shinyBSDep)
}

#' Create a progress object
#' 
#' @param message Character: progress message
#' @param divisions Integer: number of divisions in the progress bar
#' @param global Shiny's global variable
#' 
#' @export
startProgress <- function(message, divisions, global = sharedData) {
    print(message)
    global$progress.divisions <- divisions
    global$progress <- Progress$new()
    global$progress$set(message = message, value = 0)
}

#' Update a progress object
#' 
#' @details If \code{divisions} isn't NULL, a progress bar is started with the 
#' given divisions. If \code{value} is NULL, the progress bar will be 
#' incremented by one; otherwise, the progress bar will be incremented by the
#' integer given in value.
#' 
#' @inheritParams startProgress
#' @param value Integer: current progress value
#' @param max Integer: maximum progress value
#' @param detail Character: detailed message
#' 
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
    
    if (!is.null(detail))
        print(paste(message, detail, sep=": "))
    else
        print(message)
    return(invisible(TRUE))
}

#' Close the progress even if there's an error
#' 
#' @param message Character: message to show in progress bar
#' @param global Global Shiny variable where all data is stored
#' 
#' @export
closeProgress <- function(message=NULL, global = sharedData) {
    # Close the progress even if there's an error
    if (!is.null(message)) print(message)
    global$progress$close()
}

#' Allows to add id to an image
#' 
#' @param name Character: name of the icon
#' @param class Additional classes to customise the icon style
#' @param lib Icon library to use (either "font-awesome" or "glyphicon")
#' @param ... Extra arguments to the icon tag
#' 
#' @importFrom htmltools htmlDependencies
#' 
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
        htmlDependencies(iconTag) <- htmlDependency(
            "font-awesome", "4.5.0", c(href = "shared/font-awesome"), 
            stylesheet = "css/font-awesome.min.css")
    }
    iconTag
}

#' Disable a tab from the navbar
#' @importFrom shinyjs disable addClass
#' @param tab Character: tab to disable
#' @export
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
#' @export
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
#' @export
#' 
#' @examples 
#' words <- c("tumor_stage", "age", "gender")
#' textComplete("textareaid", words)
textComplete <- function(id, words, novalue = "No matching value", char=" ") {
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
#' @export
updateTextAreaInput <- updateTextInput