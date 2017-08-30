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
#' @param modal Character: modal identifier
loadRequiredData <- function( modal=NULL ) {
    modal <- ifelse(is.null(modal), "null", modal)
    return(sprintf("showDataPanel('#%s');", modal))
}

#' @rdname missingDataModal
missingDataGuide <- function(dataType) {
    js <- loadRequiredData(dataType)
    runjs(js)
}

#' @rdname appUI
#' 
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column tags
analysesUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(
        ns, "analysis", 
        priority=c("pcaUI", "diffSplicingUI", "survivalUI", "infoUI"))
    
    # Load available analyses
    analyses <- tagList()
    for (k in seq(uiList)) {
        ui <- uiList[[k]]
        if ( is(ui, "shiny.tag.list") ) {
            tabName  <- attr(ui, "name")
            analyses <- c(analyses, list(tabPanel(tabName, ui)))
        } else {
            headerName <- attr(ui, "name")
            analyses   <- c(analyses, list(headerName))
            for (each in seq(ui)) {
                name     <- attr(ui[[each]], "name")
                analyses <- c(analyses, list(tabPanel(name, ui[each])))
            }
        }
        analyses <- c(analyses, list("----")) # Add a divider after each section
    }
    analyses <- analyses[-length(analyses)]
    do.call(tab, c(list(icon="flask", title="Analyses", menu=TRUE), analyses))
}

#' Retrieve clinical data based on attributes required for survival analysis
#' 
#' @param ... Character: names of columns to retrieve
#' @param formulaStr Character: right-side of the formula for survival analysis
#' 
#' @return Filtered clinical data
getClinicalDataForSurvival <- function(..., formulaStr=NULL) {
    cols <- unlist(list(...))
    if (!is.null(formulaStr) && formulaStr != "") {
        form <- formula(paste("1 ~", formulaStr))
        cols <- c(cols, all.vars(form))
    }
    clinical <- getClinicalData(cols)
    return(clinical)
}

#' Assign the value from one of the patient's samples to that patient
#' 
#' Assign a value to patients based on the frequency of the respective type of 
#' their samples
#' 
#' @details
#' Match filtered samples with patients to retrieve values per patient. One
#' single sample is matched to a patient based on the sample type frequency. For
#' instance, imagine that:
#' 
#' \itemize{
#'     \item{10 patients have a tumour and control sample;}
#'     \item{5 patients have a tumour sample;}
#'     \item{2 patients have only a control sample;}
#'     \item{2 patients have only a metastasis sample.}
#' }
#'  
#' In total, there are 15 tumour, 12 control and 2 metastasis samples. As tumour
#' samples are the majority, tumour samples will be matched to patients. 
#' Patients without tumour samples will then be matched to control samples (2nd
#' most frequent sample type), if available. Finally, the remaining patients 
#' will be matched to metastasis samples.
#' 
#' @param data Data frame or matrix: values per sample
#' @param clinical Data frame or matrix: clinical dataset (only required if the
#' \code{patients} argument is not handed)
#' @param patients Character: patient identifiers (only required if the
#' \code{clinical} argument is not handed)
#' @param pattern Character: pattern to use when filtering sample types (NULL by
#' default, i.e. no filtering occurs)
#' @param filterOut Boolean: filter out (TRUE) or filter in (FALSE) sample types
#' based on a given pattern; by default, sample types are filtered out
#' 
#' @inheritParams matchPatientToSingleSample
#' 
#' @return Alternative splicing quantification per clinical patients
#' @export
getValuePerPatient <- function(data, match, clinical=NULL, patients=NULL,
                               pattern=NULL, filterOut=TRUE) {
    if (is.null(clinical) && is.null(patients)) {
        stop("You cannot leave both 'clinical' and 'patients' arguments ",
             "as NULL.")
    } else if (is.null(patients)) {
        patients <- rownames(clinical)
    }
    
    # Get sample identifiers of interest
    types <- parseSampleGroups(names(match))
    
    if (!is.null(pattern)) {
        # Filter sample types based on a user-defined pattern
        pattern <- paste(pattern, collapse="|")
        filter <- grepl(pattern, types)
        if (filterOut) filter <- !filter
    } else {
        filter <- TRUE
    }
    
    matchFiltered <- match[filter]
    matchFiltered <- matchFiltered[!is.na(matchFiltered)]
    
    # Assign only one sample per patient based on sample type frequency
    matchSingle <- matchPatientToSingleSample(matchFiltered)
    
    # Match samples with clinical patients (remove non-matching samples)
    clinicalValues <- data.frame(matrix(NA, nrow=nrow(data), 
                                        ncol=length(patients)))
    colnames(clinicalValues) <- patients
    rownames(clinicalValues) <- rownames(data)
    clinicalValues[ , matchSingle] <- data[ , names(matchSingle)]
    return(clinicalValues)
}

#' @rdname getValuePerPatient
getPSIperPatient <- function(psi, match, clinical=NULL, patients=NULL,
                             pattern=NULL, filterOut=TRUE) {
    .Deprecated("getValuePerPatient")
    getValuePerPatient(psi, match, clinical, patients, pattern, filterOut)
}

#' Match patients to a single sample according to sample type frequency
#'
#' Only one sample per patient is returned. For patients with more than one
#' sample, the attributed sample is chosen according to the frequency of its
#' type.
#'
#' @param match Matrix: match between samples and patients
#'
#' @return Integer containing the patient and the respective sample as its name
matchPatientToSingleSample <- function(match) {
    # Get frequency of sample types
    types <- parseSampleGroups(names(match))
    freq  <- names(sort(table(types), decreasing=TRUE))
    # Create a list of patient-sample matches based on sample types
    matchByType <- split(match, types)
    # Order the list based on the frequency of the sample types
    matchByType <- matchByType[freq]
    
    # Filter out duplicated items based on the items found on previous list
    # indexes
    filterDuplicatedItems <- function(i, aList) {
        if (i == 1) {
            diff <- TRUE
        } else {
            previous <- Reduce(union, aList[seq(i - 1)])
            diff <- !aList[[i]] %in% previous
        }
        return(aList[[i]][diff])
    }
    
    # Match patients to a single sample according to sample type frequency
    res <- unlist(lapply(seq(matchByType), filterDuplicatedItems, matchByType))
    # Remove potentially duplicated samples of the same sample type
    res <- res[!duplicated(res)]
    return(res)
}

#' Process survival data to calculate survival curves
#' 
#' @inheritParams getAttributesTime
#' @param group Character: group relative to each patient
#' @param clinical Data frame: clinical data
#' @param survTime \code{survTime} object: Times to follow up, time start, time 
#' stop and event (optional)
#' 
#' @details The event time will only be used to determine whether the event has
#' occurred (1) or not (0) in case of missing values.
#' 
#' If \code{survTime} is NULL, the survival times will be fetch from the
#' clinical dataset according to the names given in \code{timeStart},
#' \code{timeStop}, \code{event} and \code{followup}. This can became quite slow
#' when using the function in a for loop. If these variables are constant, 
#' consider running the function \code{\link{getAttributesTime}} to retrieve the
#' time of such columns once and hand the result to the \code{survTime} argument
#' of this function.
#' 
#' @return Data frame with terms needed to calculate survival curves
processSurvData <- function(event, timeStart, timeStop, followup, group, 
                            clinical, survTime=NULL) {
    if ( is.null(survTime) ) {
        survTime <- getAttributesTime(clinical, event, timeStart, timeStop,
                                      followup)
    }
    
    # Create new time using the starting time replacing the NAs with
    # days to last follow up
    nas <- is.na(survTime$start)
    survTime$time <- survTime$start
    survTime$time[nas] <- survTime$followup[nas]
    
    # Indicate event of interest and groups
    survTime$event <- ifelse(!is.na(survTime$event), 1, 0)
    if (!is.null(names(group)))
        survTime[names(group), "groups"] <- group
    else
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

#' Retrieve the time for given columns in a clinical dataset
#' 
#' @param timeStart Character: name of column containing starting time of the
#' interval or follow up time
#' @param timeStop Character: name of column containing ending time of the 
#' interval
#' @param event Character: name of column containing time of the event of
#' interest
#' @param followup Character: name of column containing follow up time
#' @param clinical Data frame: clinical data
#' 
#' @return Data frame containing the time for the given columns
#' 
#' @export
getAttributesTime <- function(clinical, event, timeStart, timeStop=NULL,
                              followup="days_to_last_followup") {
    cols <- c(followup=followup, start=timeStart, stop=timeStop, event=event)
    
    # Retrive time for given attributes
    timePerPatient <- function(col, clinical) {
        cols <- grep(col, colnames(clinical), value=TRUE)
        row  <- apply(clinical[cols], 1, function(i)
            if(!all(is.na(i))) max(as.numeric(i), na.rm = TRUE) else NA)
        return(row)
    }
    survTime <- lapply(cols, timePerPatient, clinical)
    
    survTime <- as.data.frame(survTime)
    class(survTime) <- c("data.frame", "survTime")
    return(survTime)
}

#' @inherit getAttributesTime
#' @export
getColumnsTime <- function(clinical, event, timeStart, timeStop=NULL,
                           followup="days_to_last_followup") {
    .Deprecated("getAttributesTime")
    getAttributesTime(clinical=clinical, event=event, timeStart=timeStart, 
                      timeStop=timeStop, followup=followup)
}

#' Update available clinical attributes when the clinical data changes
#' 
#' @param session Shiny session
#' @param attrs Character: patient attributes
#' 
#' @importFrom shiny observe updateSelectizeInput
#' @return NULL (this function is used to modify the Shiny session's state)
updateClinicalParams <- function(session, attrs) {
    if (!is.null(attrs)) {
        # Allow the user to select any "days_to" attribute available
        daysTo <- grep("days_to_", attrs, value=TRUE, fixed=TRUE)
        subDaysTo <- gsub(".*(days_to_.*)", "\\1", daysTo)
        choices <- unique(subDaysTo)
        names(choices) <- gsub("_", " ", choices, fixed=TRUE)
        names(choices) <- capitalize(names(choices))
        
        # Update choices for starting or follow up time
        updateSelectizeInput(
            session, "timeStart", choices=list(
                "Select a clinical attribute"="",
                "Suggested attributes"=choices,
                "All clinical attributes"=attrs),
            selected="days_to_death")
        
        # Update choices for ending time
        updateSelectizeInput(
            session, "timeStop", choices=list(
                "Select a clinical attribute"="",
                "Suggested attributes"=choices,
                "All clinical attributes"=attrs))
        
        # Update choices for events of interest
        names(choices) <- gsub("Days to ", "", names(choices), fixed=TRUE)
        names(choices) <- capitalize(names(choices))
        updateSelectizeInput(
            session, "event", choices=list(
                "Select a clinical attribute"="",
                "Suggested events"=choices,
                "All clinical attributes"=attrs),
            selected="days_to_death")
    } else {
        choices <- c("No clinical data loaded"="")
        updateSelectizeInput(session, "timeStart", choices=choices)
        updateSelectizeInput(session, "timeStop",  choices=choices)
        updateSelectizeInput(session, "event",     choices=choices)
    }
}

#' Process survival curves terms to calculate survival curves
#'
#' @inheritParams processSurvData
#' @param censoring Character: censor using "left", "right", "interval" or
#' "interval2"
#' @param scale Character: rescale the survival time to "days", "weeks",
#' "months" or "years"
#' @param formulaStr Character: formula to use
#' @param coxph Boolean: fit a Cox proportional hazards regression model? FALSE 
#' by default
#' @param survTime survTime object: times to follow up, time start, time stop
#' and event (optional)
#' 
#' @importFrom stats formula
#' @importFrom survival coxph Surv
#'
#' @details \code{timeStop} is only considered if \code{censoring} is either
#' \code{interval} or \code{interval2}
#'
#' If \code{survTime} is NULL, the survival times will be fetch from the
#' clinical dataset according to the names given in \code{timeStart},
#' \code{timeStop}, \code{event} and \code{followup}. This can became quite slow
#' when using the function in a for loop. If these variables are constant, 
#' consider running the function \code{\link{getAttributesTime}} to retrieve the
#' time of such columns once and hand the result to the \code{survTime} argument
#' of this function.
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
                             followup="days_to_last_followup", survTime=NULL) {
    # Ignore timeStop if interval-censoring is not selected
    if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "") 
        timeStop <- NULL
    
    # Check if using or not interval-censored data
    formulaSurv <- ifelse(is.null(timeStop),
                          "Surv(time/%s, event, type=censoring) ~", 
                          "Surv(time/%s, time2, event, type=censoring) ~")
    scaleStr <- scale
    scale <- switch(scaleStr, days=1, weeks=7, months=30.42, years=365.25)
    formulaSurv <- sprintf(formulaSurv, scale)
    
    survData <- processSurvData(event, timeStart, timeStop, followup, group, 
                                clinical, survTime)
    
    # Estimate survival curves by groups or using formula
    if (formulaStr == "" || is.null(formulaStr)) {
        formulaTerms <- "groups"
    } else {
        formulaTerms <- formulaStr
        factor <- sapply(clinical, is.factor)
        clinical[factor] <- lapply(clinical[factor], as.character)
        survData <- cbind(survData, clinical)
    }
    
    form <- formula(paste(formulaSurv, formulaTerms))
    
    if (coxph) {
        res <- coxph(form, data=survData)
        
        if (!is.null(res$xlevels$groups)) {
            # Correct group names
            name <- res$xlevels$groups[-1]
            names(res$means) <- name
            names(res$coefficients)  <- name
            # names(res$wald.test)  <- name
        }
    } else {
        res <- list(form=form, survTime=survData)
    }
    
    res$scale <- scaleStr
    class(res) <- c("survTerms", class(res))
    return(res)
}

#' Compute estimate of a survival curve using processed survival terms
#' 
#' @param survTerms \code{survTerms} object: processed survival terms
#' @inheritDotParams survival::survfit.formula -formula -data
#' 
#' @importFrom survival survfit
#' @method survfit survTerms
#' 
#' @return \code{survfit} object. See \code{survfit.object} for details. Methods
#' defined for survfit objects are \code{print}, \code{plot}, \code{lines}, and 
#' \code{points}.
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
    res <- survfit(survTerms$form, data=survTerms$survTime, ...)
    res$scale <- survTerms$scale
    
    # Correct group names
    groups <- deparse(survTerms$form[[3]])
    if (!is.null(res$strata) && groups == "groups") {
        name <- paste0("^", groups, "=")
        names(res$strata) <- gsub(name, "", names(res$strata))
    }
    return(res)
}

#' Test difference between two or more survival curves using processed survival 
#' terms
#' 
#' @param survTerms survTerms object: processed survival terms
#' @inheritDotParams survival::survdiff -formula -data
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
#' @param auto Boolean: return the plot automatically prepared (TRUE) or only
#' the bare minimum (FALSE)? TRUE by default
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
                               title="Survival analysis", scale=NULL,
                               auto=TRUE) {
    hc <- hchart(surv, ranges=interval, markTimes=mark)
    if (auto) {
        if (is.null(scale)) {
            if (is.null(surv$scale))
                scale <- "days"
            else
                scale <- surv$scale
        }
        
        hc <- hc %>%
            hc_chart(zoomType="xy") %>%
            hc_title(text=title) %>%
            hc_yAxis(title=list(text="Survival proportion"),
                     crosshair=TRUE) %>%
            hc_xAxis(title=list(text=paste("Time in", scale)),
                     crosshair=TRUE) %>%
            hc_tooltip(
                headerFormat = paste(
                    tags$small("{point.x:.3f}", scale), br(),
                    span(style="color:{point.color}", "\u25CF "),
                    tags$b("{series.name}"), br()),
                pointFormat = paste(
                    "Survival proportion: {point.y:.3f}", br(),
                    "Records: {series.options.records}", br(),
                    "Events: {series.options.events}", br(),
                    "Median: {series.options.median}")) %>%
            hc_plotOptions(series=list(stickyTracking=FALSE))
    }
    
    if (!is.null(pvalue))
        hc <- hc_subtitle(hc, text=paste("log-rank p-value:", pvalue))
    return(hc)
}

#' Check if survival analyses successfully completed or returned errors
#' 
#' @param session Shiny session
#' @inheritDotParams processSurvTerms
#' 
#' @importFrom shiny tags
#' @return List with survival analysis results
processSurvival <- function(session, ...) {
    # Calculate survival curves
    survTerms <- tryCatch(processSurvTerms(...), error=return)
    if ("simpleError" %in% class(survTerms)) {
        if (survTerms[[1]] == paste("contrasts can be applied only to",
                                    "factors with 2 or more levels")) {
            errorModal(session, "Formula error",
                       "Cox models can only be applied to 2 or more groups.")
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
#' @inheritDotParams survival::survdiff -formula -data
#' 
#' @note Instead of raising errors, an \code{NA} is returned
#' 
#' @return p-value of the survival difference or \code{NA}
#' @export
#' 
#' @examples 
#' require("survival")
#' data <- aml
#' timeStart  <- "event"
#' event      <- "event"
#' followup   <- "time"
#' data$event <- NA
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
        
        # Calculate p-value with given significant digits
        pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
        return(as.numeric(signifDigits(pvalue)))
    }, error = function(e) NA)
    return(pvalue)
}

#' Label groups based on a given cutoff
#' 
#' @param data Numeric: test data
#' @param cutoff Numeric: test cutoff
#' @param label Character: label to prefix group names (NULL by default)
#' @param gte Boolean: test with greater than or equal to cutoff (TRUE) or use
#' less than or equal to cutoff (FALSE)? TRUE by default
#' 
#' @return Labelled groups
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
        str1 <- "&gt;="
        str2 <- "&lt;"
    } else {
        comp <- `>`
        str1 <- "&gt;"
        str2 <- "&lt;="
    }
    group <- comp(data, cutoff)
    
    # Assign a value based on the inclusion levels cutoff
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
#' @param cutoff Numeric: Cutoff of interest
#' @param data Numeric: elements of interest to test against the cutoff
#' @param filter Boolean or numeric: elements to use (all by default)
#' @inheritDotParams processSurvTerms -group -clinical
#' @param session Shiny session
#' @param survivalInfo Boolean: return extra survival information
#' 
#' @importFrom survival survdiff
#' @return p-value of the survival difference
testSurvivalCutoff <- function(cutoff, data, filter=TRUE, clinical, ...,
                               session=NULL, survivalInfo=FALSE) {
    group <- labelBasedOnCutoff(data, cutoff, label="Inclusion levels")
    
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
    if (is.na(pvalue)) pvalue <- 1
    if (survivalInfo) attr(pvalue, "info") <- survfit(survTerms)
    return(pvalue)
}

#' Calculate optimal alternative splicing quantification cutoff to separate
#' survival curves
#'
#' @details \code{timeStop} is only considered if \code{censoring} is either
#' \code{interval} or \code{interval2}
#'
#' @inheritParams processSurvTerms
#' @inheritParams testSurvivalCutoff
#' @param psi Numeric: PSI values to test against the cutoff
#' @param session Shiny session (only used for the visual interface)
#' 
#' @return Optimal alternative splicing quantification cutoff
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
#' 
#' psi <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
#' opt <- optimalPSIcutoff(clinical, psi, "right", event, timeStart)
optimalPSIcutoff <- function(clinical, psi, censoring, event, timeStart, 
                             timeStop=NULL, followup="days_to_last_followup",
                             session=NULL, filter=TRUE, survTime=NULL) {
    if ( is.null(survTime) ) {
        survTime <- getAttributesTime(clinical, event, timeStart, timeStop,
                                      followup)
    }
    
    # Supress warnings from failed calculations while optimising
    opt <- suppressWarnings(
        optim(0, testSurvivalCutoff, data=psi, filter=filter, clinical=clinical,
              censoring=censoring, timeStart=timeStart, timeStop=timeStop, 
              event=event, followup=followup, survTime=survTime,
              session=session,
              # Method and parameters interval
              method="Brent", lower=0, upper=1))
    return(opt)
}

#' @rdname appServer
#' 
#' @importFrom shiny observe observeEvent
#' @importFrom shinyjs hide show
analysesServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions(
        "analysis", priority=c("pcaServer", "diffSplicingServer", 
                               "survivalServer", "infoServer"))
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"