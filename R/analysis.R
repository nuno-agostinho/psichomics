#' Missing information modal template
#'  
#' @param session Shiny session
#' @param dataType Character: type of data missing
#' @param buttonId Character: identifier of button to take user to load missing 
#' data
#' 
#' @examples
#' \dontrun{
#'  session <- session$ns
#'  buttonInput <- "takeMeThere"
#'  buttonId <- ns(buttonInput)
#'  dataType <- "Inclusion levels"
#'  missingDataModal(session, buttonId, dataType)
#'  observeEvent(input[[buttonInput]], missingDataGuide(dataType))
#' }
#' @return NULL (this function is used to modify the Shiny session's state)
missingDataModal <- function(session, dataType, buttonId) {
    template <- function(buttonLabel) {
        errorModal(
            session, paste("Load", tolower(dataType)),
            "This analysis requires", tolower(dataType), "to proceed.",
            footer=actionButton(buttonId, buttonLabel, "data-dismiss"="modal",
                                class="btn-danger"))
    }
    
    switch(dataType,
           "Clinical data"=template("Load"),
           "Junction quantification"=template("Load"),
           "Inclusion levels"=template("Load or calculate"))
}

#' @rdname missingDataModal
loadRequiredData <- function(dataType) {
    panel <- switch(dataType,
                    "Clinical data"="TCGA",
                    "Junction quantification"="TCGA",
                    "Inclusion levels"="Alternative splicing"
    )
    
    return(sprintf("showDataPanel('%s');", panel))
}

#' @rdname missingDataModal
missingDataGuide <- function(dataType) {
    js <- loadRequiredData(dataType)
    runjs(js)
}

#' User interface for the data analyses
#' 
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column selectizeInput conditionalPanel
#' navbarMenu
#' 
#' @return HTML element as character
analysesUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "analysis")
    
    # Load available analyses
    ui <- lapply(uiList, function(ui) tabPanel(attr(ui, "name"), ui) )
    do.call(tab, c(list(icon="flask", title="Analyses", menu=TRUE), ui))
}

#' Process survival data to calculate survival curves
#' 
#' @param timeStart Character: name of column containing starting time of the
#' interval or follow up time
#' @param timeStop Character: name of column containing ending time of the 
#' interval
#' @param event Character: name of column containing time of the event of
#' interest
#' @param followup Character: name of column containing follow up time
#' @param group Character: group of each individual
#' @param clinical Data frame: clinical data
#' 
#' @details The event time will only be used to determine whether the event has
#' happened (1) or not in case of NAs (0)
#' 
#' @return Data frame with terms needed to calculate survival curves
processSurvData <- function(event, timeStart, timeStop, followup, group, 
                            clinical) {
    cols <- c(followup = followup, start = timeStart,
              stop = timeStop, event = event)
    survTime <- lapply(cols, timePerPatient, clinical)
    survTime <- as.data.frame(survTime)
    
    # Create new time using the starting time replacing the NAs with
    # days to last follow up
    nas <- is.na(survTime$start)
    survTime$time <- survTime$start
    survTime$time[nas] <- survTime$followup[nas]
    
    # Indicate event of interest and groups
    survTime$event <- ifelse(!is.na(survTime$event), 1, 0)
    survTime$groups <- group
    
    if (!is.null(timeStop)) {
        # Create new time using the ending time replacing the NAs
        # with days to last follow up
        nas <- is.na(survTime$stop)
        survTime$time2 <- survTime$stop
        survTime$time2[nas] <- survTime$followup[nas]
    }
    return(survTime)
}

#' Get all columns matching a given string and return a single vector with the
#' max time for each patient if available
#'
#' @param col Character: column of interest
#' @param clinical Data.frame: clinical data
#'
#' @return Numeric vector with days recorded for columns of interest
timePerPatient <- function(col, clinical) {
    cols <- grep(col, names(clinical))
    row <- apply(clinical[cols], 1, function(i)
        if(!all(is.na(i))) max(as.numeric(i), na.rm = TRUE) else NA)
    return(row)
}

#' Update available clinical attributes when the clinical data changes
#' @param session Shiny session
#' @importFrom shiny observe updateSelectizeInput
#' @return NULL (this function is used to modify the Shiny session's state)
updateClinicalParams <- function(session) {
    observe({
        clinical <- getClinicalData()
        if (!is.null(clinical)) {
            # Allow the user to select any "days_to" attribute available
            daysTo <- grep("days_to_", names(clinical), value=TRUE, fixed=TRUE)
            subDaysTo <- gsub(".*(days_to_.*)", "\\1", daysTo)
            choices <- unique(subDaysTo)
            names(choices) <- gsub("_", " ", choices, fixed=TRUE)
            names(choices) <- capitalize(names(choices))
            
            # Update choices for starting time (follow up time)
            updateSelectizeInput(session, "timeStart", choices=choices,
                                 selected="days_to_death")
            
            # Update choices for ending time
            updateSelectizeInput(session, "timeStop", choices=choices)
            
            # Update choices for events of interest
            names(choices) <- gsub("Days to ", "", names(choices), fixed=TRUE)
            names(choices) <- capitalize(names(choices))
            updateSelectizeInput(session, "event",
                                 choices=list(
                                     "Suggested events"=choices,
                                     "All clinical data columns"=names(clinical)),
                                 selected="days_to_death")
        }
    })
}

#' Process survival curves terms to calculate survival curves
#' 
#' @details \code{timeStop} is only considered if \code{censoring} is either
#' \code{interval} or \code{interval2}
#'
#' @inheritParams processSurvData
#' @param censoring Character: censor using "left", "right", "interval" or
#' "interval2"
#' @param scale Character: rescale the survival time to "days", "weeks",
#' "months" or "years"
#' @param formulaStr Character: formula to use
#' @param coxph Boolean: fit a Cox proportional hazards regression model? FALSE 
#' by default
#' 
#' @importFrom stats formula
#' @importFrom survival coxph Surv
#'
#' @return A list with a \code{formula} object and a data frame with terms
#' needed to calculate survival curves
#' @export
#' 
#' @examples 
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup", 
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
processSurvTerms <- function(clinical, censoring, event, timeStart, 
                             timeStop=NULL, group=NULL, formulaStr=NULL, 
                             coxph=FALSE, scale="days",
                             followup="days_to_last_followup") {
    # Ignore timeStop if interval-censoring is not selected
    if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "") 
        timeStop <- NULL
    
    # Check if using or not interval-censored data
    formulaSurv <- ifelse(is.null(timeStop),
                          "Surv(time/%s, event, type=censoring) ~", 
                          "Surv(time/%s, time2, event, type=censoring) ~")
    scale <- switch(scale, days=1, weeks=7, months=30.42, years=365.25)
    formulaSurv <- sprintf(formulaSurv, scale)
    
    survTime <- processSurvData(event, timeStart, timeStop, followup, group, 
                                clinical)
    
    # Estimate survival curves by groups or using formula
    if (formulaStr == "" || is.null(formulaStr)) {
        formulaTerms <- "groups"
    } else {
        formulaTerms <- formulaStr
        survTime <- cbind(survTime, clinical)
    }
    
    form <- formula(paste(formulaSurv, formulaTerms))
    
    if (coxph)
        res <- coxph(form, data=survTime)
    else
        res <- list(form=form, survTime=survTime)
    class(res) <- c("survTerms", class(res))
    return(res)
}

#' Compute estime of a survival curve using processed survival terms
#' 
#' @param survTerms survTerms object: processed survival terms
#' @param ... Extra arguments passed to \code{survfit}
#' 
#' @importFrom survival survfit
#' @method survfit survTerms
#' 
#' @return an object of class "survfit". See survfit.object for details. Methods
#' defined for survfit objects are print, plot, lines, and points.
#' @export survfit.survTerms
#'
#' @examples 
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup", 
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
#' require("survival")
#' survfit(survTerms)
survfit.survTerms <- function(survTerms, ...) {
    survfit(survTerms$form, data=survTerms$survTime, ...)
}

#' Test difference between two or more survival curves using processed survival 
#' terms
#' 
#' @param survTerms survTerms object: processed survival terms
#' @param ... Extra arguments passed to \code{survdiff}
#' 
#' @importFrom survival survdiff
#' 
#' @return an object of class "survfit". See survfit.object for details. Methods
#' defined for survfit objects are print, plot, lines, and points.
#' @export
#'
#' @examples
#' clinical <- read.table(text = "2549   NA ii  female
#'                                 840   NA i   female
#'                                  NA 1204 iv    male
#'                                  NA  383 iv  female
#'                                1293   NA iii   male
#'                                  NA 1355 ii    male")
#' names(clinical) <- c("patient.days_to_last_followup", 
#'                      "patient.days_to_death",
#'                      "patient.stage_event.pathologic_stage",
#'                      "patient.gender")
#' timeStart  <- "days_to_death"
#' event      <- "days_to_death"
#' formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
#' survTerms  <- processSurvTerms(clinical, censoring="right", event, timeStart,
#'                                formulaStr=formulaStr)
#' survdiff.survTerms(survTerms)
survdiff.survTerms <- function(survTerms, ...) {
    survdiff(survTerms$form, data=survTerms$survTime, ...)
}

#' Plot survival curves
#' 
#' @param surv Survival object
#' @param interval Boolean: show interval ranges? FALSE by default
#' @param mark Boolean: mark times? TRUE by default
#' @param title Character: plot title
#' @param pvalue Numeric: p-value of the survival curves
#' @param scale Character: time scale; default is "days"
#' 
#' @importFrom shiny tags br
#' 
#' @return Plot of survival curves
#' @export
#' @examples 
#' require("survival")
#' fit <- survfit(Surv(time, status) ~ x, data = aml)
#' plotSurvivalCurves(fit)
plotSurvivalCurves <- function(surv, mark=TRUE, interval=FALSE, pvalue=NULL, 
                               title="Survival analysis", scale="days") {
    hc <- hchart(surv, ranges=interval, markTimes=mark) %>%
        hc_chart(zoomType="xy") %>%
        hc_title(text=title) %>%
        hc_yAxis(title=list(text="Proportion of individuals"),
                 crosshair=TRUE) %>%
        hc_xAxis(title=list(text=paste("Time in", scale)), crosshair=TRUE) %>%
        hc_tooltip(
            headerFormat = paste(
                tags$small("{point.x}", scale), br(),
                span(style="color:{point.color}", "\u25CF "),
                tags$b("{series.name}"), br()),
            pointFormat = paste(
                "Records: {series.options.records}", br(),
                "Events: {series.options.events}", br(),
                "Median: {series.options.median}")) %>%
        hc_plotOptions(series=list(stickyTracking=FALSE))
    
    if (!is.null(pvalue))
        hc <- hc_subtitle(hc, text=paste("log-rank p-value:", pvalue))
    return(hc)
}

#' Check if survival analyses successfully completed or returned errors
#' 
#' @param session Shiny session
#' @param ... Arguments to pass to function \code{processSurvTerms}
#' 
#' @importFrom shiny tags
#' @return List with survival analysis results
processSurvival <- function(session, ...) {
    # Calculate survival curves
    survTerms <- tryCatch(processSurvTerms(...), error=return)
    if ("simpleError" %in% class(survTerms)) {
        if (survTerms[[1]] == "The formula field can't be empty") {
            errorModal(session, "Error in formula",
                       "The formula field can't be empty.")
        } else {
            errorModal(session, "Formula error",
                       "Maybe you misplaced a ", tags$kbd("+"), ", ",
                       tags$kbd(":"), " or ", tags$kbd("*"), "?", br(),
                       br(), "The following error was raised:", br(), 
                       tags$code(survTerms$message))
        }
        return(NULL)
    }
    return(survTerms)
}

#' Test the survival difference between survival groups
#' 
#' @inheritParams survdiff.survTerms
#' 
#' @note Instead of raising errors, an NA is returned
#' 
#' @return p-value of the survival difference or NA
#' @export
#' 
#' @examples 
#' require("survival")
#' data <- aml
#' timeStart  <- "event"
#' event      <- "event"
#' followup   <- "time"
#' data$event  <- NA
#' data$event[aml$status == 1] <- aml$time[aml$status == 1]
#' censoring  <- "right"
#' formulaStr <- "x"
#' survTerms <- processSurvTerms(data, censoring=censoring, event=event, 
#'                               timeStart=timeStart, followup=followup,
#'                               formulaStr=formulaStr)
#' testSurvival(survTerms)
testSurvival <- function (survTerms, ...) {
    # If there's an error with survdiff, return NA
    pvalue <- tryCatch({
        # Test the difference between survival curves
        diff <- survdiff.survTerms(survTerms, ...)
        
        # Calculate p-value with 5 significant numbers
        pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
        return(as.numeric(signifDigits(pvalue)))
    }, error = function(e) NA)
    return(pvalue)
}

#' Label groups based on a given cut-off
#' 
#' @param data Numeric: test data
#' @param cutoff Numeric: test cutoff
#' @param label Character: label to prefix group names (NULL by default)
#' @param gte Boolean: test with greater than or equal to cutoff (TRUE) or use
#' less than or equal to cutoff (FALSE)? TRUE by default
#' 
#' @return Labeled groups
#' @export
#' 
#' @examples
#' labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5)
#' 
#' labelBasedOnCutoff(data=c(1, 0, 0, 1, 0, 1), cutoff=0.5, "Ratio")
#'                    
#' # Use "greater than" instead of "greater than or equal to"
#' labelBasedOnCutoff(data=c(1, 0, 0, 0.5, 0, 1), cutoff=0.5, gte=FALSE)
labelBasedOnCutoff <- function (data, cutoff, label=NULL, gte=TRUE) {
    len <- length(data)
    group <- rep(NA, len)
    
    if (gte) {
        comp <- `>=`
        str1 <- ">="
        str2 <- "<"
    } else {
        comp <- `>`
        str1 <- ">"
        str2 <- "<="
    }
    group <- comp(data, cutoff)

    # Assign a value based on the inclusion levels cut-off
    if (is.null(label)) {
        group[group == "TRUE"]  <- paste(str1, cutoff)
        group[group == "FALSE"] <- paste(str2, cutoff)
    } else {
        group[group == "TRUE"]  <- paste(label, str1, cutoff)
        group[group == "FALSE"] <- paste(label, str2, cutoff)
    }
    
    length(group) <- len
    return(group)
}

#' Test the survival difference between two survival groups given a cutoff
#' 
#' @inheritParams processSurvTerms
#' @param cutoff Numeric: Cut-off of interest
#' @param data Numeric: elements of interest to test against the cut-off
#' @param group Pre-filled vector of missing values with the length of data
#' @param filter Boolean or numeric: elements to use (all by default)
#' @param ... Arguments to pass to \code{processSurvTerms}
#' @param session Shiny session
#' 
#' @importFrom survival survdiff
#' @return p-value of the survival difference
testSurvivalCutoff <- function(cutoff, data, filter=TRUE, clinical, ...,
                               group=NULL, session=NULL) {
    if (is.null(group)) group <- rep(NA, nrow(clinical))
    group[filter] <- data[filter] >= cutoff
    
    # Assign a value based on the inclusion levels cut-off
    group[group == "TRUE"]  <- paste("Inclusion levels >=", cutoff)
    group[group == "FALSE"] <- paste("Inclusion levels <", cutoff)
    
    # Calculate survival curves
    if (!is.null(session)) {
        survTerms <- processSurvival(session, group=group, clinical=clinical,
                                     ...)
        if (is.null(survTerms)) return(NULL)
    } else {
        survTerms <- tryCatch(processSurvTerms(group=group, clinical=clinical, 
                                               ...), 
                              error=return)
        if ("simpleError" %in% class(survTerms)) return(NA)
    }
    
    pvalue <- testSurvival(survTerms)
    return(pvalue)
}

#' Server logic for the analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny observe observeEvent updateSelectizeInput
#' @importFrom shinyjs hide show
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
analysesServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("analysis")
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"