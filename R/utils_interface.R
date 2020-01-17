#' Sidebar without a well
#' 
#' Modified version of \code{shiny::sidebarPanel} without a well
#' 
#' @importFrom shiny div tags
#' 
#' @inherit shiny::sidebarPanel
#' 
#' @return HTML elements
#' @keywords internal
sidebar <- function(..., width=4) {
    div(class = paste0("col-sm-", width), tags$form(...))
}

#' Link to run arbitrary JavaScript code
#' 
#' @param text Character: text label
#' @param code Character: JavaScript code
#' 
#' @return HTML elements
#' @keywords internal
linkToRunJS <- function(text, code) {
    HTML(sprintf('<a href="#" onclick="%s; return false;">%s</a>', code, text))
}

#' Create a row for a HTML table
#' 
#' @param ... Elements to include in the row
#' @param th Boolean: is this row the table head?
#' 
#' @return HTML elements
#' @keywords internal
tableRow <- function (..., th=FALSE) {
    args <- list(...)
    if (th) row <- tags$th
    else    row <- tags$td
    do.call(tags$tr, lapply(args, row))
}

#' Modified colour input with 100\% width
#' 
#' @inheritDotParams colourpicker::colourInput
#' @importFrom colourpicker colourInput
#' 
#' @return HTML elements
#' @keywords internal
colourInputMod <- function(...) {
    colourSelector <- colourInput(...)
    colourSelector[[2]][["style"]] <- "width: 100%;"
    return(colourSelector)
}


#' Create word break opportunities (for HTML) using given characters
#' 
#' @param str Character: text
#' @param pattern Character: pattern(s) of interest to be used as word break
#' opportunities
#' 
#' @importFrom shiny HTML
#' 
#' @return String containing HTML elements
#' @keywords internal
prepareWordBreak <- function(str, pattern=c(".", "-", "\\", "/", "_", ",", 
                                            " ")) {
    res <- str
    # wbr: word break opportunity
    for (p in pattern) res <- gsub(p, paste0(p, "<wbr>"), res, fixed=TRUE)
    return(HTML(res))
}

#' Convert vector of values to JavaScript array
#' 
#' @param values Character vector
#' 
#' @return Character with valid JavaScript array
#' @keywords internal
toJSarray <- function(values) {
    paste0("[", paste0(paste0("\'", values, "\'"), collapse=", "), "]")
}

#' Style button used to initiate a process
#' 
#' @param id Character: button identifier
#' @param label Character: label
#' @inheritDotParams shiny::actionButton -inputId -label
#' @param class Character: class
#' 
#' @importFrom shinyjs hidden
#' @importFrom shiny tags actionButton
#' 
#' @return HTML for a button
#' @keywords internal
processButton <- function(id, label, ..., class="btn-primary") {
    spinner <- tags$i(id=paste0(id, "Loading"), class="fa fa-spinner fa-spin")
    button  <- actionButton(id, class=class, type="button", 
                            label=div(hidden(spinner), label), ...)
    return(button)
}

#' Create an icon based on set operations
#' 
#' Based on the \code{\link[shiny]{icon}()} function
#' 
#' @param name Character: icon name
#' @param class Character: additional classes to customise the icon element
#' @param ... Extra arguments for the icon HTML element
#' 
#' @importFrom shiny icon
#' @importFrom htmltools htmlDependency htmlDependencies htmlDependencies<-
#' 
#' @return Icon element
#' @keywords internal
setOperationIcon <- function (name, class=NULL, ...) {
    if (length(list(...)) == 0) {
        style <- paste("font-size: 20px;", "line-height: 0;",
                       "vertical-align: bottom;", "display: inline-block;")
    } else {
        style <- NULL
    }
    
    prefix <- "set"
    iconClass <- ""
    if (!is.null(name)) 
        iconClass <- paste0(prefix, " ", prefix, "-", name)
    if (!is.null(class)) 
        iconClass <- paste(iconClass, class)
    iconTag <- tags$i(class=iconClass, style=style, ...)
    htmlDependencies(iconTag) <- htmlDependency(
        "set-operations", "1.0",
        c(href="set-operations"), stylesheet = "css/set-operations.css")
    return(iconTag)
}



#' Set the status of a process to style a given button
#' 
#' \itemize{
#'   \item{\code{startProcess}: Style button to show a process is in progress}
#'   \item{\code{endProcess}: Style button to show a process finished; also, 
#'   closes the progress bar (if \code{closeProgressbar = TRUE}) and 
#'   prints the difference between the current time and \code{time}}
#' }
#' 
#' @param id Character: button identifier
#' @importFrom shinyjs show
#' 
#' @return \code{startProcess} returns the start time of the process (may be 
#' used as the \code{time} argument to \code{endProcess}), whereas
#' \code{endProcess} returns the difference between current time and \code{time}
#' (or \code{NULL} if \code{time} is not specified)
#' @keywords internal
startProcess <- function(id) {
    disable(id)
    show(paste0(id, "Loading"))
    return(Sys.time())
}

#' @rdname startProcess
#' 
#' @param time \code{POSIXct} object: start time needed to show the interval
#' time (if \code{NULL}, the time interval is not displayed)
#' @param closeProgressBar Boolean: close progress bar?
#' 
#' @importFrom shinyjs enable hide
endProcess <- function(id, time=NULL, closeProgressBar=TRUE) {
    enable(id)
    hide(paste0(id, "Loading"))
    if (closeProgressBar) suppressWarnings(closeProgress())
    if (!is.null(time)) {
        diffTime <- Sys.time() - time
        display(diffTime, "Process finished in")
        return(diffTime)
    }
    return(NULL)
}


#' Create a modal window
#'
#' @param session Shiny session
#' @param title Character: title
#' @inheritDotParams shiny::modalDialog -title -size -footer
#' @param style Character: style (\code{NULL}, \code{warning}, \code{error} or
#'   \code{info})
#' @param iconName Character: icon name
#' @param footer HTML elements to use in footer
#' @param echo Boolean: print to console?
#' @param size Character: size of the modal (\code{small}, \code{medium} or
#'   \code{large})
#' @param dismissButton Boolean: show dismiss button in footer?
#' @param caller Character: caller module identifier
#'
#' @importFrom shiny renderUI div icon showModal modalButton modalDialog
#' @importFrom shinyBS toggleModal
#' @importFrom R.utils capitalize
#'
#' @seealso \code{\link{showAlert}()}
#' @inherit psichomics return
#' @keywords internal
styleModal <- function(session, title, ..., style=NULL,
                       iconName="exclamation-circle", footer=NULL, echo=FALSE, 
                       size="medium", dismissButton=TRUE, caller=NULL) {
    size <- switch(size, "small"="s", "large"="l", "medium"="m")
    if (dismissButton) footer <- tagList(modalButton("Dismiss"), footer)
    
    modal <- modalDialog(..., title=div(icon(iconName), title), size=size,
                         footer=footer, easyClose=FALSE)
    if (!is.null(style)) {
        style <- match.arg(style, c("info", "warning", "error"))
        modal[[3]][[1]][[3]][[1]][[3]][[1]] <-
            tagAppendAttributes(modal[[3]][[1]][[3]][[1]][[3]][[1]],
                                class=style)
    }
    showModal(modal, session)
    if (echo) {
        if (style == "info") style <- "Information"
        msg <- sprintf("%s: %s", capitalize(style), title)
        if (!is.null(caller)) msg <- sprintf('%s [in "%s"]', msg, caller)
        message(msg)
    }
    return(invisible(TRUE))
}

#' @rdname styleModal
errorModal <- function(session, title, ..., size="small", footer=NULL, 
                       caller=NULL) {
    styleModal(session, title, ..., footer=footer, style="error", size=size,
               echo=TRUE, iconName="times-circle", caller=caller)
}

#' @rdname styleModal
warningModal <- function(session, title, ..., size="small", footer=NULL,
                         caller=NULL) {
    styleModal(session, title, ..., footer=footer, style="warning", size=size,
               echo=TRUE, iconName="exclamation-circle", caller=caller)
}

#' @rdname styleModal
infoModal <- function(session, title, ..., size="small", footer=NULL,
                      caller=NULL) {
    styleModal(session, title, ..., footer=footer, style="info", size=size,
               echo=TRUE, iconName="info-circle", caller=caller)
}

#' Show or remove an alert
#' 
#' @inheritParams styleModal
#' @param ... Arguments to render as elements of alert
#' @param style Character: style (\code{error}, \code{warning} or \code{NULL})
#' @param dismissible Boolean: is the alert dismissible?
#' @param alertId Character: identifier
#' 
#' @seealso \code{\link{showModal}()}
#' @importFrom shiny span h3 renderUI div tagList
#' 
#' @inherit psichomics return
#' @keywords internal
showAlert <- function(session, ..., title, style=NULL, dismissible=TRUE, 
                      alertId="alert", iconName=NULL, caller=NULL) {
    if (dismissible) {
        dismissible <- "alert-dismissible"
        dismiss <- tags$button(type="button", class="close",
                               "data-dismiss"="alert", "aria-label"="Close",
                               span("aria-hidden"="true", "\u00D7"))
    } else {
        dismissible <- NULL
        dismiss <- NULL
    }
    
    # Log information
    args <- list(...)
    if (style == "info") style <- "Information"
    msg <- sprintf("%s: %s", capitalize(style), title)
    if (!is.null(caller)) msg <- sprintf('%s [in "%s"]', msg, caller)
    body <- paste(lapply(args, format), collapse=" ")
    
    newline <- "\n  "
    processHTML <- function(str, newline) {
        str <- gsub("\n[ ]*", " ", str) # Strip forced newlines
        str <- gsub("[ ]*<br[/]{0,1}>[ ]*", newline, str) # Convert newline
        str <- gsub("[ ]*<.*?>[ ]*", "", str) # Strip other HTML tags
        return(str)
    }
    body <- processHTML(body, newline)
    message(msg, newline, body)
    
    style <- switch(style, "error"="alert-danger", "warning"="alert-warning",
                    "success"="alert-success")
    
    output <- session$output
    output[[alertId]] <- renderUI({
        tagList(div(h4(icon(iconName), title), id="myAlert", class="alert",
                    class=style, role="alert", class="animated bounceInUp", 
                    class=dismissible, dismiss, ...))
    })
}

#' @rdname showAlert
successAlert <- function(session, ..., title=NULL, dismissible=TRUE,
                         alertId="success", caller=NULL) {
    showAlert(session, ..., style="success", title=title, 
              iconName="check-circle", dismissible=dismissible, 
              alertId=alertId, caller=caller)
}

#' @rdname showAlert
errorAlert <- function(session, ..., title=NULL, dismissible=TRUE,
                       alertId="alert", caller=NULL) {
    showAlert(session, ..., style="error", title=title, 
              iconName="times-circle", dismissible=dismissible, 
              alertId=alertId, caller=caller)
}

#' @rdname showAlert
warningAlert <- function(session, ..., title=NULL, dismissible=TRUE,
                         alertId="alert", caller=NULL) {
    showAlert(session, ..., style="warning", title=title, 
              iconName="exclamation-circle", dismissible=dismissible, 
              alertId=alertId, caller=caller)
}

#' @rdname showAlert
#' 
#' @param output Shiny output
removeAlert <- function(output, alertId="alert") {
    output[[alertId]] <- renderUI(NULL)
}

#' Alert in the style of a dialogue box with a button
#' 
#' @param id Character: identifier
#' @param description Character: description
#' @param buttonId Character: button identifier
#' @param buttonLabel Character: button label
#' @param buttonIcon Character: button icon
#' @param ... Extra parameters when creating the alert
#' @param type Character: type of alert (error or warning)
#' @param bigger Boolean: wrap the \code{description} in a \code{h4} tag?
#'
#' @importFrom shiny icon div actionButton
#'
#' @return HTML elements
#' @keywords internal
inlineDialog <- function(description, ..., buttonLabel=NULL, buttonIcon=NULL, 
                         buttonId=NULL, id=NULL, type=c("error", "warning"),
                         bigger=FALSE) {
    type <- match.arg(type)
    if (identical(type, "error")) type <- "danger"
    typeIcon <- switch(type, danger="exclamation-circle",
                       warning="exclamation-triangle")
    
    if (!is.null(buttonLabel)) {
        if (!is.null(buttonIcon))
            icon <- icon(buttonIcon)
        else
            icon <- NULL
        
        typeClass <- sprintf("btn-%s btn-block", type)
        button <- tagList(br(), br(), actionButton(buttonId, icon=icon, 
                                                   buttonLabel,
                                                   class=typeClass))
    } else {
        button <- NULL
    }
    
    if (bigger) {
        description <- h4(style="margin-top: 5px !important;",
                          icon(typeIcon), description)
    } else {
        description <- tagList(icon(typeIcon), description)
    }
    
    typeClass <- sprintf("alert alert-%s", type)
    div(id=id, class=typeClass, role="alert", style="margin-bottom: 0px;",
        description, button, ...)
}

#' @rdname inlineDialog
errorDialog <- function(description, ...)
    inlineDialog(description, ..., type="error")

#' @rdname inlineDialog
warningDialog <- function(description, ...)
    inlineDialog(description, ..., type="warning")


#' Display characters in the command-line
#' 
#' @param char Character: message
#' @param timeStr Character: message when a \code{difftime} object is passed to
#' the \code{char} argument
#' 
#' @importFrom shiny isRunning
#' 
#' @return \code{NULL} (display message in command-line)
#' @keywords internal
display <- function(char, timeStr="Time difference of") {
    if (!isRunning()) cat("", fill=TRUE)
    if (is(char, "difftime")) {
        message(timeStr, " ", format(unclass(char), digits=3), " ", 
                attr(char, "units"))
    } else {
        cat(char, fill=TRUE)
    }
}

#' Create, set and terminate a progress object
#' 
#' @param message Character: progress message
#' @param divisions Integer: number of divisions in the progress bar
#' @param global Shiny's global variable
#' 
#' @importFrom shiny isRunning Progress
#' @importFrom utils txtProgressBar
#' 
#' @inherit psichomics return
#' @keywords internal
startProgress <- function(message, divisions,
                          global=if (isRunning()) sharedData else getHidden()) {
    display(message)
    if (isRunning()) {
        global$progress <- Progress$new()
        global$progress$set(message = message, value = 0)
    } else {
        global$progress <- txtProgressBar(style=3)
    }
    global$progress.divisions <- divisions
    return(invisible(global))
}

#' @rdname startProgress
#' 
#' @details If \code{divisions} is not \code{NULL}, a progress bar starts with 
#' the given divisions. If \code{value = NULL}, the progress bar increments one
#' unit; otherwise, the progress bar increments \code{value}.
#' 
#' @param value Integer: current progress value
#' @param max Integer: maximum progress value
#' @param detail Character: detailed message
#' @param console Boolean: print message to console?
#' 
#' @importFrom shiny isRunning Progress
#' @importFrom utils setTxtProgressBar
updateProgress <- function(message="Loading...", value=NULL, max=NULL, 
                           detail=NULL, divisions=NULL, 
                           global=if (isRunning()) sharedData else getHidden(),
                           console=TRUE) {
    isGUIversion <- isRunning()
    if (!interactive()) return(NULL)
    if (!is.null(divisions)) {
        if (!isGUIversion)
            setHidden(startProgress(message, divisions, new.env()))
        else
            startProgress(message, divisions, global)
        return(NULL)
    }
    
    divisions <- global$progress.divisions
    if (is.null(value)) {
        if (!isRunning()) { # CLI version
            currentValue <- global$progress$getVal()
            max   <- 1
        } else {
            currentValue <- global$progress$getValue()
            max          <- global$progress$getMax()
        }
        value <- currentValue + (max - currentValue)
    }
    amount <- ifelse(is.null(max), value/divisions, 1/max/divisions)
    
    # Print message to console
    if (console) {
        msg <- message
        if (!is.null(detail) && !identical(detail, ""))
            msg <- paste(msg, detail, sep=": ")
        display(msg)
    }
    
    # Increment progress
    if (!isGUIversion) {
        if (!is.null(global)) {
            value <- min(global$progress$getVal() + amount, 1)
            setTxtProgressBar(global$progress, value)
            setHidden(global)
        }
    } else {
        if (is.null(detail)) detail <- ""
        global$progress$inc(amount=amount, message=message, detail=detail)
    }
    return(invisible(TRUE))
}

#' @rdname startProgress
#' @importFrom shiny isRunning Progress
closeProgress <- function(message=NULL, 
                          global=if (isRunning()) sharedData else getHidden()) {
    # Close the progress even if there's an error
    if (!is.null(message)) display(message)
    
    isGUIversion <- isRunning()
    if (isGUIversion)
        global$progress$close()
    else
        close(global$progress)
}

#' Modified version of \code{shinyBS::bsModal}
#' 
#' \code{bsModal} is used within the UI to create a modal window. This allows to
#' modify the modal footer.
#' 
#' @inheritParams shinyBS::bsModal
#' @param footer UI set: List of elements to include in the footer
#' @param style Character: message style can be \code{warning}, \code{error}, 
#' \code{info} or \code{NULL}
#' @param size Character: Modal size (\code{small}, \code{default} or 
#' \code{large})
#' 
#' @importFrom shiny tagAppendAttributes
#' @importFrom shinyBS bsModal
#' 
#' @return HTML elements
#' @keywords internal
bsModal2 <- function (id, title, trigger, ..., size=NULL, footer=NULL, 
                      style = NULL)  {
    if (is.null(size))
        modal <- bsModal(id, title, trigger, ...)
    else
        modal <- bsModal(id, title, trigger, ..., size=size)
    
    if (!is.null(style)) {
        style <- match.arg(style, c("info", "warning", "error"))
        modal[[3]][[1]][[3]][[1]][[3]][[1]] <-
            tagAppendAttributes(modal[[3]][[1]][[3]][[1]][[3]][[1]],
                                class=style)
    }
    
    modal[[3]][[1]][[3]][[1]][[3]][[3]] <-
        tagAppendChild(modal[[3]][[1]][[3]][[1]][[3]][[3]], footer)
    return(modal)
}

#' Enable or disable a tab from the \code{navbar}
#' 
#' @param tab Character: tab
#' 
#' @importFrom shinyjs disable addClass
#' 
#' @inherit psichomics return
#' @keywords internal
disableTab <- function(tab) {
    # Style item as disabled
    addClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
             class = "disabled")
    # Disable link itself
    disable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' @rdname disableTab
#' @importFrom shinyjs removeClass enable
enableTab <- function(tab) {
    # Style item as enabled
    removeClass(selector = paste0(".navbar li:has(a[data-value=", tab, "])"),
                class = "disabled")
    # Enable link itself
    enable(selector = paste0(".navbar li a[data-value=", tab, "]"))
}

#' Create script for auto-completion of text input
#' 
#' Uses the JavaScript library \code{jquery.textcomplete}
#' 
#' @param id Character: input ID
#' @param words Character: words to suggest
#' @param novalue Character: string when there's no matching values
#' @param char Character to succeed accepted word
#'
#' @return HTML string with the JavaScript script prepared to run
#' @keywords internal
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

#' Render a data table with sparkline HTML elements
#' 
#' @details This slightly modified version of \code{\link{renderDataTable}()}
#' calls a JavaScript function to convert the sparkline HTML elements to an
#' interactive \code{highchart} object
#' 
#' @inheritDotParams shiny::renderDataTable -options -escape -env
#' @param options List of options to pass to \code{\link{renderDataTable}()}
#' 
#' @importFrom DT renderDataTable JS
#' 
#' @inherit psichomics return
#' @keywords internal
renderDataTableSparklines <- function(..., options=NULL) {
    # Escape is set to FALSE to render the Sparkline HTML elements
    renderDataTable(..., escape=FALSE, env=parent.frame(n=1), options=c(
        list(drawCallback=JS("drawSparklines")), options))
}
