# Interactive tests to perform:
# - fileBrowser()
# - fileBrowser("~/Desktop")
# - fileBrowser(caption="Select a file")
# - fileBrowser(directory=TRUE)
# - fileBrowser(multiple=TRUE) # select one file
# - fileBrowser(multiple=TRUE) # select two or more files
# - fileBrowser(multiple=TRUE, directory=TRUE) # select one directory
# - fileBrowser(multiple=TRUE, directory=TRUE) # select two or more directories

chooseFilesLinux <- function(default, caption, multiple, directory) {
    directory <- ifelse(directory, "--directory", "")
    multiple  <- ifelse(multiple,  "--multiple",  "")
    sep       <- " &sep& "
    separator <- sprintf('--separator="%s"', sep)

    # Default location
    if (!is.null(default) && nzchar(default)) {
        default <- gsub("\"", "\\\\\"", path.expand(default))
        default <- sprintf('--filename="%s/"', default)
    } else {
        default <- ""
    }

    prompt <- ""
    if (!is.null(caption) && nzchar(caption)) {
        prompt <- sprintf("--title='%s'", caption)
    }
    args <- " --file-selection %s %s %s %s %s"
    args <- sprintf(args, directory, multiple, prompt, separator, default)
    path <- suppressWarnings(system2("zenity", args=args, stderr=TRUE))

    # Return NA if user cancels the action
    if (!is.null(attr(path, "status")) && attr(path, "status")) return(NA)

    # Error: Gtk-Message: GtkDialog mapped without a transient parent
    if(length(path) == 2) path <- path[2]

    # Vector of multiple files/directories
    path <- strsplit(path, sep)[[1]]
    return(path)
}

chooseFilesMac <- function(default, caption, multiple, directory) {
    directory <- ifelse(directory, "folder", "file")
    multiple  <- ifelse(multiple, "with multiple selections allowed", "")

    if (!is.null(caption) && nzchar(caption)) {
        prompt <- sprintf("with prompt \\\"%s\\\"", caption)
    } else {
        prompt <- ""
    }

    # Default location
    if (!is.null(default) && nzchar(default)) {
        default <- sprintf("default location \\\"%s\\\"",
                           path.expand(default))
    } else {
        default <- ""
    }

    app  <- "path to frontmost application as text"
    args <- '-e "tell app (%s) to set thePaths to (choose %s %s %s %s)"'
    # Get POSIX paths of selected files
    args <- paste(args,
                  '-e "if class of thePaths is not list"',
                  '-e "log POSIX path of thePaths"',
                  '-e "else"',
                  '-e "repeat with eachPath in thePaths"',
                  '-e "log POSIX path of eachPath"',
                  '-e "end repeat"',
                  '-e "end if"')
    args <- sprintf(args, app, directory, multiple, prompt, default)
    args <- trimWhitespace(args)

    path <- suppressWarnings(system2("osascript", args=args, stderr=TRUE))

    # Return NA if the user cancels the action
    if (!is.null(attr(path, "status")) && attr(path, "status")) return(NA)

    return(path)
}

#' Interactive folder selection using a native dialogue
#'
#' @param default Character: path to initial folder
#' @param caption Character: caption on the selection dialogue
#' @param multiple Boolean: allow to select multiple files?
#' @param directory Boolean: allow to select directories instead of files?
#'
#' @details
#' Platform-dependent implementation:
#' \itemize{
#'  \item{\strong{Windows}: calls the \code{utils::choose.files} R function.}
#'  \item{\strong{macOS}: uses AppleScript to display a folder selection
#'  dialogue. If \code{default = NA}, folder selection falls back to the
#'  default behaviour of the \code{choose folder} AppleScript command.
#'  Otherwise, paths are expanded with \code{\link{path.expand}()}.}
#'  \item{\strong{Linux}: calls the \code{zenity} system command.}
#' }
#'
#' @source \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @return A length one character vector, character NA if 'Cancel' was selected
#' @keywords internal
fileBrowser <- function(default=NULL, caption=NULL, multiple=FALSE,
                        directory=FALSE) {
    system <- Sys.info()['sysname']
    if (is.null(system)) {
        stop("File browser is unsupported in this system")
    } else if (isRStudioServer()) {
        stop("File browser is currently unsupported for RStudio Server")
    } else if (system == 'Darwin') {
        path <- chooseFilesMac(default, caption, multiple, directory)
    } else if (system == 'Linux') {
        path <- chooseFilesLinux(default, caption, multiple, directory)
    } else if (system == "Windows") {
        if (is.null(default)) default <- ""
        if (is.null(caption)) caption <- ""

        if (directory) {
            path <- utils::choose.dir(default, caption)
        } else {
            path <- utils::choose.files(default, caption, multiple)
        }
    }

    if (identical(path, "")) path <- NULL
    if (length(path) > 1) path <- paste(path, collapse=" && ")
    return(path)
}

#' @importFrom shinyBS bsPopover
bsPopoverMod <- function(...) {
    showInfo <- bsPopover(...)
    constructor <- paste(sprintf(
        "$.fn.popover.Constructor.DEFAULTS.whiteList.%s = [];",
        c("kbd", "table", "tr", "td", "th", "thead", "tbody")), collapse=" ")
    showInfo[[3]][[1]] <- gsub(
        "$(document).ready(function() {",
        paste("$(document).ready(function() {", constructor),
        showInfo[[3]][[1]], fixed=TRUE)
    return(showInfo)
}

createInfoAndClearButtons <- function(id, info=FALSE, clearable=FALSE) {
    if (info) {
        infoId   <- paste0(id, "-info")
        infoElem <- actionButton(inputId=infoId, label=NULL,
                                 style="background: #eee",
                                 icon=icon("question-circle"))
    } else {
        infoId   <- id
        infoElem <- NULL
    }

    if (clearable) {
        clearId   <- paste0(id, "-clear")
        clearElem <- actionButton(inputId=clearId, label=NULL,
                                  style="background: #eee",
                                  icon=icon("times-circle"))
    } else {
        clearId   <- id
        clearElem <- NULL
    }

    if (info || clearable) {
        buttonsElem <- div(class="input-group-btn", clearElem, infoElem)
    } else {
        buttonsElem <- NULL
    }
    return(buttonsElem)
}

#' @importFrom htmltools tagQuery
showExtraInfo <- function(id, FUN, title, content, placement) {
    id      <- paste0(id, "-info")
    title   <- gsub("\n", "", as.character(title),   fixed=TRUE)
    content <- gsub("\n", "", as.character(content), fixed=TRUE)
    if (identical(FUN, bsPopover)) {
        showInfo <- bsPopoverMod(id, placement=placement,
                                 options=list(container="body"),
                                 title=title, content=content)
    } else if (identical(FUN, bsTooltip)) {
        showInfo <- FUN(id, placement=placement,
                        options=list(container="body"), title=title)
    } else {
        showInfo <- NULL
    }
    return(showInfo)
}

#' File browser input
#'
#' Input to interactively select a file or directory on the server
#'
#' @param id Character: input identifier
#' @param label Character: input label (if \code{NULL}, no labels are displayed)
#' @param value Character: initial value (paths are expanded via
#' \code{\link{path.expand}()})
#' @param placeholder Character: placeholder when no file or folder is selected
#' @param info Boolean: add information icon for tooltips and pop-overs
#' @param infoFUN Function to use to provide information (e.g.
#' \code{shinyBS::bsTooltip} and \code{shinyBS::bsPopover})
#' @param infoPlacement Character: placement of the information (top, bottom,
#' right or left)
#' @param infoTitle Character: text to show as title of information
#' @param infoContent Character: text to show as content of information
#' @param clearable Boolean: allow to clear selected file or directory?
#'
#' @details
#' To show the dialog for file input, the \code{\link{prepareFileBrowser}()}
#' function needs to be included in the server logic.
#'
#' This widget relies on \code{\link{fileBrowser}()} to present an interactive
#' dialogue to users for selecting a directory on the local filesystem.
#' Therefore, this widget is intended for shiny apps that are run locally - i.e.
#' on the same system that files/directories are to be accessed - and not from
#' hosted applications (e.g. from \url{https://www.shinyapps.io}).
#'
#' @importFrom shinyBS bsPopover bsTooltip
#'
#' @source \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @return HTML elements for a file browser input
#' @keywords internal
#'
#' @seealso
#' \code{\link{updateFileBrowserInput}()} and \code{\link{prepareFileBrowser}()}
fileBrowserInput <- function(id, label, value=NULL, placeholder=NULL,
                             info=FALSE, infoFUN=NULL, infoPlacement="right",
                             infoTitle="", infoContent="", clearable=FALSE) {
    if (!is.null(value) && !is.na(value)) value <- path.expand(value)
    if (is.null(placeholder)) placeholder <- ""

    buttonsElem <- createInfoAndClearButtons(id, info, clearable)
    if (info) {
        showInfo <- showExtraInfo(id, infoFUN, infoTitle, infoContent,
                                  infoPlacement)
    } else {
        showInfo <- NULL
    }

    fileBrowserButton <- div(class="btn btn-info fileBrowser-input",
                             id=sprintf("%sButton", id), 'Browse...')
    if (isRStudioServer()) fileBrowserButton <- disabled(fileBrowserButton)
    fileBrowserButton <- div(class="input-group-btn", fileBrowserButton)
    filepathInput <- tags$input(
        id=id, value=value, type='text', placeholder=placeholder,
        # readonly = if (!isRStudioServer()) 'readonly' else NULL,
        class='form-control fileBrowser-input-chosen-dir')

    check <- function (x, y) {
        # Based on shiny:::`%AND%`
        if (!is.null(x) && !is.na(x)) {
            if (!is.null(y) && !is.na(y$children[[1]])) {
                return(y)
            }
        }
        return(NULL)
    }

    tagList(
        div(class='form-group fileBrowser-input-container',
            check(label, tags$label(label)),
            div(class='input-group shiny-input-container', style='width:100%;',
                fileBrowserButton, filepathInput, buttonsElem)),
        showInfo)
}

fileBrowserShinyproxyInput <- function(id, label, info=FALSE, infoFUN=NULL,
                                       infoPlacement="right", infoTitle="",
                                       infoContent="", clearable=FALSE) {
    buttonsElem <- createInfoAndClearButtons(id, info, clearable)
    if (info) {
        showInfo <- showExtraInfo(id, infoFUN, infoTitle, infoContent,
                                  infoPlacement)
    } else {
        showInfo <- NULL
    }

    input <- fileInput(id, label)
    input <- replaceStrInList(input, "btn-default", "btn-info")

    input <- tagQuery(input)$find(".input-group")$append(buttonsElem)$allTags()
    return(tagList(input, showInfo))
}

#' Change the value of a \code{\link{fileBrowserInput}()} on the client
#'
#' @param session Shiny session
#' @param id Character: identifier
#' @param ... Additional arguments passed to \code{\link{fileBrowser}()}. Only
#' used if \code{value = NULL}.
#' @param value Character: file or directory path
#' @param ask Boolean: ask user to pick a file using file browser?
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object. For \code{\link{fileBrowserInput}()} objects, this changes the value
#' displayed in the text-field and triggers a client-side change event. A
#' directory selection dialogue is not displayed.
#'
#' @source \url{https://github.com/wleepang/shiny-directory-input}
#'
#' @inherit psichomics return
#' @keywords internal
updateFileBrowserInput <- function(session, id, ..., value=NULL, ask=FALSE) {
    if (ask) value <- fileBrowser(...)
    button <- sprintf("%sButton", id)
    session$sendInputMessage(button, list(path=value))
}

#' Prepare file browser dialogue and update the input's value accordingly to
#' selected file or directory
#'
#' @inheritParams appServer
#' @param id Character: input identifier
#' @inheritDotParams fileBrowser
#' @param modalId Character: modal window identifier
#'
#' @inherit psichomics return
#' @keywords internal
prepareFileBrowser <- function(session, input, id, modalId="modal", ...) {
    buttonId <- sprintf("%sButton", id)
    observeEvent(input[[buttonId]], {
        if (input[[buttonId]] > 0) { # Prevent execution on initial launch
            errorTitle <- NULL
            if (is.null(Sys.info())) {
                errorTitle <- c(
                    "The file browser is not supported for this system")
            } else if (isRStudioServer()) {
                errorTitle <- c(
                    "The file browser is not supported in RStudio Server")
            } else {
                updateFileBrowserInput(session, id, ..., ask=TRUE)
            }
            if (!is.null(errorTitle)) {
                errorModal(session, errorTitle,
                           "Please use instead the text input field to type",
                           "the full path to the file or folder of interest.",
                           modalId=modalId, caller="File browser")
            }
        }
    })

    # Clear file selection
    clearId <- paste0(id, "-clear")
    observeEvent(input[[clearId]], {
        if (input[[clearId]] > 0) { # Prevent execution on initial launch
            updateTextInput(session, id, value="")
        }
    })
}
