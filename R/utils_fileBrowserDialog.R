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
#'  dialogue. If \code{default = NA}, folder selection fallbacks to the
#'  default behaviour of the \code{choose folder} AppleScript command.
#'  Otherwise, paths are expanded with \code{\link{path.expand}()}.}
#'  \item{\strong{Linux}: calls the \code{zenity} system command.}
#' }
#' 
#' If for some reason an error occurs (e.g. when using a remote server), the
#' dialog fallbacks to an alternative, non-native file browser.
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
        directory <- ifelse(directory, "folder", "file")
        multiple  <- ifelse(multiple, "with multiple selections allowed", "")
        
        if (!is.null(caption) && nzchar(caption))
            prompt <- sprintf("with prompt \\\"%s\\\"", caption)
        else
            prompt <- ""
        
        # Default location
        if (!is.null(default) && nzchar(default)) {
            default <- sprintf("default location \\\"%s\\\"",
                               path.expand(default))
        } else {
            default <- ""
        }
        
        app  <- "path to frontmost application as text"
        args <- '-e "tell app (%s) to POSIX path of (choose %s %s %s %s)"'
        args <- sprintf(args, app, directory, multiple, prompt, default)
        
        path <- suppressWarnings(system2("osascript", args=args, stderr=TRUE))
        
        # Return NA if the user cancels the action
        if (!is.null(attr(path, "status")) && attr(path, "status")) return(NA)
    } else if (system == 'Linux') {
        directory <- ifelse(directory, "--directory", "")
        multiple  <- ifelse(multiple,  "--multiple", "")
        
        prompt <- ""
        if (!is.null(caption) && nzchar(caption))
            prompt <- sprintf("--title='%s'", caption)
        
        args <- " --file-selection %s %s %s"
        args <- sprintf(args, directory, multiple, prompt)
        path <- suppressWarnings(system2("zenity", args=args, stderr=TRUE))
        
        # Return NA if user cancels the action
        if (!is.null(attr(path, "status")) && attr(path, "status")) return(NA) 
        
        # Error: Gtk-Message: GtkDialog mapped without a transient parent
        if(length(path) == 2) path <- path[2]
    } else if (system == "Windows") {
        if (is.null(default)) default <- ""
        if (is.null(caption)) caption <- ""
        path <- utils::choose.files(default, caption, !directory && multiple)
    }
    
    if (identical(path, "")) path <- NULL
    return(path)
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
    
    check <- function (x, y) {
        # Based on shiny:::`%AND%`
        if (!is.null(x) && !is.na(x)) 
            if (!is.null(y) && !is.na(y)) 
                return(y)
        return(NULL)
    }
    
    if (info) {
        infoId   <- paste0(id, "-info")
        infoElem <- div(class="input-group-addon", id=infoId,
                        icon("question-circle"))
    } else {
        infoId   <- id
        infoElem <- NULL
    }
    
    if (clearable) {
        clearId   <- paste0(id, "-clear")
        clearElem <- div(class="input-group-addon", id=clearId,
                         icon("times-circle"))
    } else {
        clearId   <- id
        clearElem <- NULL
    }
    
    infoTitle   <- gsub("\n", "", as.character(infoTitle),   fixed=TRUE)
    infoContent <- gsub("\n", "", as.character(infoContent), fixed=TRUE)
    if (identical(infoFUN, bsPopover)) {
        showInfo <- infoFUN(infoId, placement=infoPlacement, 
                            options=list(container="body"), 
                            title=infoTitle, content=infoContent)
    } else if (identical(infoFUN, bsTooltip)) {
        showInfo <- infoFUN(infoId, placement=infoPlacement, 
                            options=list(container="body"), title=infoTitle)
    } else {
        showInfo <- NULL
    }
    
    fileBrowserButton <- div(class="btn btn-default fileBrowser-input",
                             id=sprintf("%sButton", id), 'Browse...')
    if (isRStudioServer()) fileBrowserButton <- disabled(fileBrowserButton)
    fileBrowserButton <- div(class="input-group-btn", fileBrowserButton)
    filepathInput <- tags$input(
        id=id, value=value, type='text', placeholder=placeholder,
        # readonly = if (!isRStudioServer()) 'readonly' else NULL,
        class='form-control fileBrowser-input-chosen-dir')
    
    tagList(
        div(class='form-group fileBrowser-input-container',
            check(label, tags$label(label)),
            div(class='input-group shiny-input-container', style='width:100%;',
                fileBrowserButton, filepathInput, clearElem, infoElem)),
        showInfo)
}

#' Change the value of a \code{\link{fileBrowserInput}()} on the client
#'
#' @param session Shiny session
#' @param id Character: input identifier
#' @param value Character: file or directory path
#' @param ... Additional arguments passed to \code{\link{fileBrowser}()}. Only
#' used if \code{value = NULL}.
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
updateFileBrowserInput <- function(session, id, ..., value=NULL) {
    if (is.null(value)) value <- fileBrowser(...)
    
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
            if (is.null(Sys.info()))
                errorTitle <- c(
                    "The file browser is not supported for this system")
            else if (isRStudioServer())
                errorTitle <- c(
                    "The file browser is not supported in RStudio Server")
            else
                updateFileBrowserInput(session, id, ...)
            
            if (!is.null(errorTitle)) {
                errorModal(session, errorTitle, 
                           "Please use instead the text input field to type",
                           "the full path to the file or folder of interest.", 
                           modalId=modalId, caller="File browser")
            }
        }
    })
}
