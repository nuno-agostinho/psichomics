## TODO(NunoA): Should groups be merged if there's an intersection?
## TODO(NunoA): How to correctly do interval censoring?

#' User interface of survival analysis
#' 
#' @param id Character: namespace identifier
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS tagList uiOutput sidebarPanel radioButtons helpText hr 
#' selectizeInput checkboxInput sliderInput actionButton mainPanel
#' conditionalPanel
#' 
#' @return Character with HTML
survivalUI <- function(id) {
    ns <- NS(id)
    
    tagList(
        uiOutput(ns("modal")),
        sidebarPanel(
            radioButtons(ns("censoring"), "Data censoring", selected="right",
                         inline=TRUE, choices=c(Left="left",
                                                Right="right",
                                                Interval="interval",
                                                "Interval 2" = "interval2")),
            selectizeInput(ns("timeStart"), choices = NULL, "Follow up time"),
            # If the chosen censoring contains the word 'interval', show this input
            conditionalPanel(
                paste0("input[id='", ns("censoring"),
                       "'].indexOf('interval') > -1"),
                selectizeInput(ns("timeStop"), choices=NULL, "Ending time")),
            helpText("In case there's no record for a patient, the days to last",
                     "follow up will be used instead."),
            selectizeInput(ns("event"), choices = NULL, "Event of interest"),
            hr(),
            radioButtons(ns("modelTerms"), selected="groups", inline=TRUE,
                         div("Select groups for survival analysis",
                             icon("question-circle")),
                         choices=c("Clinical groups"="groups",
                                   "Clinical groups (interaction)"="formula",
                                   "Inclusion levels cut-off"="psiCutoff")),
            bsTooltip(ns("modelTerms"), placement="right", 
                      options = list(container = "body"),
                      paste(
                          "Perform survival analysis using:<br/>\u2022",
                          "User-created <b>clinical groups</b><br/>\u2022",
                          "A formula that can test clinical attributes with",
                          "<b>interactions</b><br/>\u2022 <b>Inclusion levels",
                          "cut-off</b> for the selected alternative splicing",
                          "event")),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "groups"),
                selectGroupsUI(ns("dataGroups"), "Clinical groups to plot"),
                checkboxInput(ns("showOutGroup"), 
                              "Show data outside chosen groups",
                              value = FALSE)),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "formula"),
                textAreaInput(
                    ns("formula"), "Formula with clinical attributes", 
                    placeholder="Start typing to suggest clinical attributes"),
                uiOutput(ns("formulaSuggestions")),
                helpText(
                    "To analyse a series of attributes, separate each attribute",
                    "with a", tags$kbd("+"), ". To analyse interactions, use", 
                    tags$kbd(":"), " (interactions are only usable with Cox",
                    "models). For example, ",
                    tags$kbd("pathologic_stage : gender + race"), br(), br(),
                    "Interesting attributes include", tags$b("pathologic_stage"), 
                    "to get tumour stages.")),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "psiCutoff"),
                uiOutput(ns("optimalPsi"))),
            hr(),
            radioButtons(ns("scale"), "Display time in", inline=TRUE,
                         c(Days="days", Weeks="weeks", Months="months",
                           Years="years")),
            checkboxInput(ns("markTimes"), "Show time marks", value = FALSE),
            checkboxInput(ns("ranges"), "Show interval ranges", value = FALSE),
            actionButton(ns("coxModel"), "Fit Cox PH model"),
            actionButton(ns("survivalCurves"), class="btn-primary",
                         "Plot survival curves")
        ),
        mainPanel(
            highchartOutput(ns("survival")),
            uiOutput(ns("coxphUI"))
        )
    )
}

#' Process survival data to calculate survival curves
#' 
#' @param timeStart Numeric: starting time of the interval or follow up time
#' @param timeStop Numeric: ending time of the interval
#' @param event Numeric: time of the event of interest
#' @param group Character: group of each individual
#' @param clinical Data.frame: clinical data
#' 
#' @details The event time will only be used to determine whether the event has
#' happened (1) or not in case of NAs (0)
#' 
#' @return Data frame with terms needed to calculate survival curves
processSurvData <- function(timeStart, timeStop, event, group, clinical) {
    cols <- c(followup = "days_to_last_followup", start = timeStart,
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

#' Process survival curves terms to calculate survival curves
#'
#' @param group Character: group of each individual 
#' @param clinical Data frame: clinical data
#' @param censoring Character: censor using "left", "right", "interval" or
#' "interval2"
#' @param timeStart Numeric: starting time
#' @param timeStop Numeric: ending time (needed only for interval-censored data)
#' @param scale Character: rescale the survival time to "days", "weeks",
#' "months" or "years"
#' @param dataEvent Character: event of interest
#' @param modelTerms Character: use "groups", "formula" or "psiCutoff" for the 
#' survival curves?
#' @param formulaStr Character: formula to use
#' @param coxph Boolean: fit a Cox proportional hazards regression model? FALSE 
#' by default
#' 
#' @importFrom stats formula
#' @importFrom survival coxph Surv
#'
#' @return A list with a \code{formula} object and a data frame with terms
#' needed to calculate survival curves
processSurvTerms <- function(group, clinical, censoring, timeStart, timeStop, 
                             dataEvent, modelTerms, formulaStr=NULL, 
                             coxph=FALSE, scale="days") {
    # Ignore timeStop if interval-censoring is not selected
    if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "") 
        timeStop <- NULL
    
    # Check if using or not interval-censored data
    formulaSurv <- ifelse(is.null(timeStop),
                          "Surv(time/%s, event, type=censoring) ~", 
                          "Surv(time/%s, time2, event, type=censoring) ~")
    scale <- switch(scale, days=1, weeks=7, months=30.42, years=365.25)
    formulaSurv <- sprintf(formulaSurv, scale)
    
    survTime <- processSurvData(timeStart, timeStop, dataEvent, group, clinical)
    
    # Estimate survival curves by groups or using formula
    if (modelTerms == "groups" || modelTerms == "psiCutoff") {
        formulaTerms <- "groups"
    } else if (modelTerms == "formula") {
        formulaTerms <- formulaStr
        
        if (formulaTerms == "" || is.null(formulaTerms))
            stop("The formula field can't be empty")
        
        survTime <- cbind(survTime, clinical)
    }
    
    form <- formula(paste(formulaSurv, formulaTerms))
    
    if (coxph)
        res <- coxph(form, data=survTime)
    else
        res <- list(form=form, survTime=survTime)
    return(res)
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
#' @param ... Arguments to pass to \code{survdiff}
#' 
#' @importFrom survival survdiff
#' @return p-value of the survival difference or NA any error occurs
testSurvival <- function(...) {
    # If there's an error with survdiff, return NA
    pvalue <- tryCatch({
        # Test the difference between survival curves
        diff <- survdiff(...)
        
        # Calculate p-value with 5 significant numbers
        pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
        return(as.numeric(signifDigits(pvalue)))
    }, error = function(e) NA)
    return(pvalue)
}

#' Test the survival difference between two survival groups given a cutoff
#' 
#' @inheritParams processSurvTerms
#' @param cutoff Numeric: Cut-off of interest
#' @param data Numeric: attribute of interest of the clinical data
#' @param group Pre-filled vector of NAs with the length of data
#' @param filter Boolean or numeric: interest of the data elements
#' @param ... Arguments to pass to \code{processSurvTerms}
#' @param session Shiny session
#' @param modals Boolean: show error dialogs? TRUE by default
#' 
#' @importFrom survival survdiff
#' @return p-value of the survival difference
testSurvivalCutoff <- function(cutoff, data, group, filter, ..., session=NULL,
                               modals=TRUE) {
    group[filter] <- data >= cutoff
    
    # Assign a value based on the inclusion levels cut-off
    group[group == "TRUE"]  <- paste("Inclusion levels >=", cutoff)
    group[group == "FALSE"] <- paste("Inclusion levels <", cutoff)
    
    # Calculate survival curves
    if (modals) {
        survTerms <- processSurvival(session, group, ...,
                                     modelTerms="psiCutoff")
        if (is.null(survTerms)) return(NULL)
    } else {
        survTerms <- tryCatch(processSurvTerms(group, ...,
                                               modelTerms="psiCutoff"),
                              error=return)
        if ("simpleError" %in% class(survTerms)) return(NA)
    }
    
    pvalue <- testSurvival(survTerms$form, data=survTerms$survTime)
    return(pvalue)
}

#' Server logic of survival analysis
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom R.utils capitalize
#' @importFrom shiny renderUI observe updateSelectizeInput observeEvent isolate 
#' br tagList hr div icon
#' @importFrom stats pchisq optim
#' @importFrom survival survfit survdiff
#' @importFrom highcharter hchart hc_chart hc_yAxis hc_xAxis hc_tooltip
#' hc_subtitle hc_tooltip renderHighchart
#' @importFrom DT dataTableOutput renderDataTable
survivalServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups", getClinicalData(), 
                       "Clinical data")
    
    # Update available clinical data attributes to use in a formula
    output$formulaSuggestions <- renderUI({
        attributes <- names(getClinicalData())
        textSuggestions(ns("formula"), attributes)
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        label <- "Follow up time"
        if (grepl("interval", input$censoring, fixed=TRUE))
            label <- "Starting time"
        updateSelectizeInput(session, "timeStart", label=label)
    })
    
    # Update available clinical attributes when the clinical data changes
    updateClinicalParams(session)
    
    observeEvent(input$missingClinical, missingDataGuide("Clinical data"))
    observeEvent(input$missingInclusionLevels,
                 missingDataGuide("Inclusion levels"))
    
    # Plot survival curves
    observeEvent(input$survivalCurves, {
        isolate({
            # Get user input
            clinical   <- getClinicalData()
            timeStart  <- input$timeStart
            timeStop   <- input$timeStop
            dataEvent  <- input$event
            censoring  <- input$censoring
            outGroup   <- input$showOutGroup
            modelTerms <- input$modelTerms
            formulaStr <- input$formula
            intRanges  <- input$ranges
            markTimes  <- input$markTimes
            psi        <- getInclusionLevels()
            event      <- getEvent()
            psiCutoff  <- input$psiCutoff
            scale      <- input$scale
            # Get chosen groups
            chosen <- input$dataGroups
            dataGroups <- getGroupsFrom("Clinical data")[chosen, , drop=FALSE]
        })
        
        if (is.null(clinical)) {
            missingDataModal(session, "Clinical data", ns("missingClinical"))
        } else if (modelTerms == "groups" && nrow(dataGroups) > 0 &&
                   anyDuplicated(unlist(dataGroups[, "Rows"])) > 0) {
            rows <- dataGroups[, "Rows"]
            len <- sapply(rows, length)
            ns <- rep(names(len), len)
            dup <- duplicated(unlist(rows)) | duplicated(unlist(rows), 
                                                         fromLast=TRUE)
            intersectingGroups <- unique(ns[dup])
            
            # If the chosen groups have any intersections
            errorModal(session, "Clinical groups intersect",
                       "There is an intersection between the selected groups.",
                       "If you'd like to perform complex queries, use the",
                       "'Formula' field.", br(), br(),
                       "The following groups have intersecting values:",
                       tags$kbd(paste(intersectingGroups, collapse = ", ")))
        } else {
            if (modelTerms == "groups") {
                # Assign one group for each clinical patient
                fillGroups <- groupPerPatient(dataGroups, nrow(clinical),
                                              outGroup)
            } else if (modelTerms == "psiCutoff") {
                if (is.null(psi)) {
                    missingDataModal(session, "Inclusion levels",
                                     ns("missingInclusionLevels"))
                    return(NULL)
                } else if (is.null(event)) {
                    errorModal(session, "No event selected",
                               "Select an alternative splicing event.")
                    return(NULL)
                }
                
                # Get tumour sample IDs (normal and control samples are not
                # interesting for survival analysis)
                match <- getClinicalMatchFrom("Inclusion levels")
                types <- getSampleTypes(names(match))
                tumour <- match[!grepl("Normal|Control", types)]
                
                # Group samples by the inclusion levels cut-off
                clinicalIDs <- nrow(clinical)
                groups <- rep(NA, clinicalIDs)
                psi <- as.numeric(psi[event, toupper(names(tumour))])
                groups[tumour] <- psi >= psiCutoff
                
                # Assign a value based on the inclusion levels cut-off
                # groups[is.na(groups)] <- "NA"
                groups[groups == "TRUE"]  <- paste("Inclusion levels >=",
                                                   psiCutoff)
                groups[groups == "FALSE"] <- paste("Inclusion levels <",
                                                   psiCutoff)
                fillGroups <- groups
            } else {
                fillGroups <- "All data"
            }
            
            survTerms <- processSurvival(session, fillGroups, clinical,
                                         censoring, timeStart, timeStop, 
                                         dataEvent, modelTerms, formulaStr, 
                                         scale=scale)
            if (is.null(survTerms)) return(NULL)
            surv <- tryCatch(survfit(survTerms$form, data=survTerms$survTime), 
                             error = return)
            if ("simpleError" %in% class(surv)) {
                errorModal(session, "Formula error",
                           "The following error was raised:", br(),
                           tags$code(surv$message))
                return(NULL)
            }
            
            pvalue <- testSurvival(survTerms$form, data=survTerms$survTime)
            
            # Plot survival curves
            output$survival <- renderHighchart({
                hc <- hchart(surv, ranges=intRanges, markTimes=markTimes) %>%
                    hc_chart(zoomType="xy") %>%
                    hc_yAxis(title=list(text="Proportion of individuals")) %>%
                    hc_xAxis(title=list(text=paste("Time in", scale))) %>%
                    hc_tooltip(headerFormat=paste(
                        "{point.x}", scale, "<br>")) %>%
                    hc_subtitle(text=paste("log-rank p-value:", pvalue)) %>%
                    hc_tooltip(crosshairs=TRUE)
            })
        }
    })
    
    # Plot cox model
    observeEvent(input$coxModel, {
        isolate({
            # Get user input
            clinical   <- getClinicalData()
            timeStart  <- input$timeStart
            timeStop   <- input$timeStop
            dataEvent  <- input$event
            censoring  <- input$censoring
            outGroup   <- input$showOutGroup
            modelTerms <- input$modelTerms
            formulaStr <- input$formula
            intRanges  <- input$ranges
            scale      <- input$scale
            # Get chosen groups
            chosen <- input$dataGroups
            dataGroups <- getGroupsFrom("Clinical data")[chosen, , drop=FALSE]
        })
        
        if (is.null(clinical)) {
            missingDataModal(session, "Clinical data", ns("missingClinical"))
        } else if (nrow(dataGroups) > 0 &&
                   anyDuplicated(unlist(dataGroups[, "Rows"])) > 0) {
            # If the chosen groups have any intersections
            errorModal(session, "Clinical groups intercept",
                       "There is an interception between clinical groups.")
        } else {
            if (modelTerms == "groups") {
                # Assign one group for each clinical patient
                dataGroups <- groupPerPatient(dataGroups, nrow(clinical),
                                              outGroup)
            } else {
                dataGroups <- NULL
            }
            
            # Calculate survival curves
            survTerms <- processSurvival(session, dataGroups, clinical, 
                                         censoring, timeStart, timeStop, 
                                         dataEvent, modelTerms, formulaStr,
                                         coxph=TRUE, scale=scale)
            if (is.null(survTerms)) return(NULL)
            summary <- summary(survTerms)
            print(summary)
            
            output$coxphUI <- renderUI({
                tagList(
                    dataTableOutput(ns("coxTests")), hr(),
                    dataTableOutput(ns("coxGroups"))
                )
            })
            
            output$coxGroups <- renderDataTable({
                groups <- cbind(signifDigits(summary$coefficients),
                                signifDigits(summary$conf.int[ , 2:4]))
                return(groups)
            }, style="bootstrap", selection='none',
            options=list(scrollX=TRUE))
            
            output$coxTests <- renderDataTable({
                tests <- rbind("Wald test"=summary$waldtest,
                               "Log test"=summary$logtest,
                               "Score (logrank) test"=summary$sctest)
                colnames(tests) <- c("Value", "Degrees of freedom", "p-value")
                return(tests)
            }, style="bootstrap", selection='none',
            options=list(info=FALSE, paging=FALSE, searching=FALSE, 
                         scrollX=TRUE))
        }
    })
    
    # Calculate optimal inclusion levels
    output$optimalPsi <- renderUI({
        # Get user input
        clinical   <- getClinicalData()
        timeStart  <- input$timeStart
        timeStop   <- input$timeStop
        dataEvent  <- input$event
        censoring  <- input$censoring
        psi        <- getInclusionLevels()
        event      <- getEvent()
        
        if (is.null(clinical)) {
            return(helpText(icon("exclamation-circle"), 
                            "Please, load clinical data."))
        } else if (is.null(getInclusionLevels())) {
            return(helpText(icon("exclamation-circle"),
                            "Please, load or calculate the quantification of",
                            "alternative splicing events."))
        } else if (is.null(getEvent()) || getEvent() == "") {
            return(helpText(icon("exclamation-circle"), 
                            "Please, select an alternative splicing event."))
        } else {
            # Get tumour sample IDs (normal and control samples are not
            # interesting for survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- getSampleTypes(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Group samples by the inclusion levels cut-off
            clinicalIDs <- nrow(clinical)
            groups <- rep(NA, clinicalIDs)
            psi <- as.numeric(psi[event, toupper(names(tumour))])
            
            # Supress warnings from failed calculations while optimising
            opt <- suppressWarnings(
                optim(0, testSurvivalCutoff, group=groups, filter=tumour,
                      data=psi, clinical=clinical, censoring=censoring, 
                      timeStart=timeStart, timeStop=timeStop, 
                      dataEvent=dataEvent, session=session,
                      # Method and parameters interval
                      method="Brent", lower=0, upper=1))
            
            slider <- tagList(
                sliderInput(ns("psiCutoff"), value = 0.5, min=0, max=1,
                            step=0.01, paste(
                                "Quantification cut-off for the selected",
                                "splicing event")),
                bsTooltip(ns("psiCutoff"), placement="right", 
                          options = list(container = "body"),
                          paste("You can click on the white circle and then",
                                "use the left and right arrows for finer",
                                "control.")))
            
            observe({
                if (!is.na(opt$value)) 
                    value <- opt$par
                else
                    value <- 0.5
                updateSliderInput(session, "psiCutoff", value=value)
            })
            
            if (!is.na(opt$value))
                return(tagList(
                    slider,
                    div(tags$b("Optimal quantification cut-off:"), opt$par, 
                        br(), tags$b("Minimal log-rank p-value:"), opt$value)))
            else
                return(tagList(
                    slider, div(icon("bell-o"), "No optimal cut-off was found",
                                "for this alternative splicing event.")))
        }
    })
}

attr(survivalUI, "loader") <- "analysis"
attr(survivalUI, "name") <- "Survival curves"
attr(survivalUI, "selectEvent") <- TRUE
attr(survivalServer, "loader") <- "analysis"