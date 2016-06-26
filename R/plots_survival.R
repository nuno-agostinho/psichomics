## TODO(NunoA): Should groups be merged if there's an intersection?
## TODO(NunoA): How to correctly do interval censoring?

#' User interface of the survival curves
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
                paste0("input[id='", ns("censoring"), "'].indexOf('interval') > -1"),
                selectizeInput(ns("timeStop"), choices = NULL, "Ending time")),
            helpText("In case there's no record for a patient, the days to last",
                     "follow up will be used instead."),
            selectizeInput(ns("event"), choices = NULL, "Event of interest"),
            radioButtons(ns("modelTerms"), selected="groups", inline=TRUE,
                         "Select model terms of the right-hand using",
                         choices=c("Clinical groups"="groups",
                                   "Formula"="formula",
                                   "Inclusion leves cutoff"="psiCutoff")),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "groups"),
                fluidRow(
                    column(10, selectizeInput(ns("dataGroups"),
                                              "Clinical groups to use",
                                              choices = NULL, multiple = TRUE)),
                    column(2, actionButton(ns("dataGroupsEdit"), "Edit",
                                           class="inline_selectize"))),
                checkboxInput(ns("showOutGroup"), "Show data outside chosen groups",
                              value = FALSE)),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "formula"),
                textAreaInput(ns("formula"), "Formula for right-hand side"),
                uiOutput(ns("formulaAutocomplete")),
                helpText("Interesting attributes include", 
                         tags$b("pathologic_stage"))),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "psiCutoff"),
                numericInput(ns("psiCutoff"),  value = 0.5, step=0.01,
                             "Cutoff value for the selected event")),
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
#' @param groups Character: group of each individual
#' @param clinical Data.frame: clinical data
#' 
#' @details The event time will only be used to determine whether the event has
#' happened (1) or not in case of NAs (0)
#' 
#' @return Data frame with terms needed to calculate survival curves
processSurvData <- function(timeStart, timeStop, event, groups, clinical) {
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
    survTime$groups <- groups
    
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
#' @param session Session object from Shiny function
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
processSurvTerms <- function(session, group, clinical, censoring, timeStart, 
                             timeStop, dataEvent, modelTerms, formulaStr, 
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
        if (formulaTerms == "" || is.null(formulaTerms)) {
            errorModal(session, "Error in formula",
                       "The formula field can't be empty.")
            return(NULL)
        }
        survTime <- cbind(survTime, clinical)
    }
    
    form <- tryCatch(formula(paste(formulaSurv, formulaTerms)), error = return)
    if ("simpleError" %in% class(form)) {
        errorModal(session, "Formula error",
                   "Maybe you misplaced a ", tags$kbd("+"), ", ", tags$kbd(":"), 
                   " or ", tags$kbd("*"), "?", br(), br(),  
                   "The following error was raised:", br(), 
                   tags$code(form$message))
        return(NULL)
    }
    
    if (coxph)
        res <- tryCatch(coxph(form, data=survTime), error=return)
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

updateClinicalFields <- function(session) {
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
            updateSelectizeInput(
                session, "timeStop", choices = choices, options=list(
                    onInitialize = I('function() { this.setValue(""); }')))
            
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

#' @importFrom R.utils capitalize
#' @importFrom stats pchisq
#' @importFrom survival survfit survdiff
#' @importFrom highcharter hchart hc_chart hc_yAxis hc_xAxis hc_tooltip
#' hc_subtitle hc_tooltip renderHighchart
survivalServer <- function(input, output, session) {
    ns <- session$ns
    
    # Update available clinical data attributes to use in a formula
    output$formulaAutocomplete <- renderUI({
        attributes <- names(getClinicalData())
        textComplete(ns("formula"), attributes)
    })
    
    # Update available group choices to select
    observe({
        groups <- getGroupsFrom("Clinical data")
        updateSelectizeInput(
            session, "dataGroups", choices=groups[, "Names"])
    #         options = list(placeholder =
    #                            ifelse(length(groups) > 0,
    #                                   "Click 'Select all' to select all groups",
    #                                   "No groups created")))
    })

    # # Select all data groups when pressing the respective "Select all" button
    # observeEvent(input$dataGroups_selectAll, {
    #     updateSelectizeInput(
    #         session, "dataGroups",
    #         selected = getGroupsFrom("Clinical data")[, "Names"])
    # })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        label <- "Follow up time"
        if (grepl("interval", input$censoring, fixed=TRUE))
            label <- "Starting time"
        updateSelectizeInput(session, "timeStart", label=label)
    })
    
    # Update available clinical attributes when the clinical data changes
    updateClinicalFields(session)
    
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
            errorModal(session, "Clinical data missing",
                       "Insert clinical data first.")
        } else if (modelTerms == "groups" && nrow(dataGroups) > 0 &&
                   anyDuplicated(unlist(dataGroups[, "Rows"])) > 0) {
            # If the chosen groups have any intersections
            errorModal(session, "Clinical groups intercept",
                       "There is an interception between clinical groups.")
        } else {
            if (modelTerms == "groups") {
                # Assign one group for each clinical patient
                fillGroups <- groupPerPatient(dataGroups, nrow(clinical),
                                              outGroup)
            } else if (modelTerms == "psiCutoff") {
                if (is.null(psi)) {
                    errorModal(session, "Inclusion levels missing",
                               "You need to calculate or load inclusion levels first.")
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
                
                # Group samples by the inclusion levels cutoff
                clinicalIDs <- nrow(clinical)
                groups <- rep(NA, clinicalIDs)
                psi <- as.numeric(psi[event, toupper(names(tumour))])
                groups[tumour] <- psi >= psiCutoff
                
                # Assign a value based on the inclusion levels cutoff
                # groups[is.na(groups)] <- "NA"
                groups[groups == "TRUE"]  <- paste("Inclusion levels >=",
                                                   psiCutoff)
                groups[groups == "FALSE"] <- paste("Inclusion levels <",
                                                   psiCutoff)
                fillGroups <- groups
            } else {
                fillGroups <- "All data"
            }
            
            # Calculate survival curves
            survTerms <- processSurvTerms(session, fillGroups, clinical,
                                          censoring, timeStart, timeStop,
                                          dataEvent, modelTerms, formulaStr,
                                          scale = scale)
            form <- survTerms$form
            data <- survTerms$survTime
            surv <- tryCatch(survfit(form, data = data),
                             error = return)
            
            if ("simpleError" %in% class(surv)) {
                errorModal(session, "Formula error",
                           "The following error was raised:", br(),
                           tags$code(surv$message))
                return(NULL)
            }
            
            # If there's an error with survdiff, show NA
            pvalue <- tryCatch({
                # Test the difference between survival curves
                diff <- survdiff(form, data = data)
                
                # Calculate p-value with 5 significant numbers
                pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
                signifDigits(pvalue)
            }, error = function(e) NA)
            
            # Plot survival curves
            output$survival <- renderHighchart({
                hc <- hchart(surv, ranges=intRanges, markTimes=markTimes) %>%
                    hc_chart(zoomType="xy") %>%
                    hc_yAxis(title=list(text="Proportion of individuals")) %>%
                    hc_xAxis(title=list(text=paste("Time in", scale))) %>%
                    hc_tooltip(headerFormat='Time: {point.x}<br>') %>%
                    hc_subtitle(text=paste("p-value:", pvalue)) %>%
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
            errorModal(session, "Clinical data missing",
                       "Insert clinical data first.")
        } else if (nrow(dataGroups) > 0 &&
                   anyDuplicated(unlist(dataGroups[, "Rows"])) > 0) {
            # If the chosen groups have any intersections
            errorModal(session, "Clinical groups intercept",
                       "There is an interception between clinical groups.")
        } else {
            # Calculate survival curves
            survTerms <- processSurvTerms(session, dataGroups, clinical,
                                          censoring, timeStart, timeStop,
                                          dataEvent, modelTerms, formulaStr,
                                          coxph=TRUE, scale=scale)
            if ("simpleError" %in% class(survTerms)) {
                errorModal(session, "Formula error",
                           "The following error was raised:", br(),
                           tags$code(survTerms$message))
                return(NULL)
            }
            
            surv <- survfit(survTerms)
            summary <- summary(survTerms)
            print(summary)
            
            output$coxphUI <- renderUI({
                # highchartOutput(ns("coxPlot"))
                tagList(
                    dataTableOutput(ns("coxGroups")),
                    dataTableOutput(ns("coxTests"))
                )
            })
            
            output$coxGroups <- renderDataTable({
                groups <- cbind(rownames(summary$coefficients),
                                signifDigits(summary$coefficients),
                                signifDigits(summary$conf.int[ , 2:4]))
                return(groups)
            }, options = list(scrollX = TRUE))
            
            output$coxTests <- renderDataTable({
                tests <- rbind("Wald test"=summary$waldtest,
                               "Log test"=summary$logtest,
                               "Score (logrank) test"=summary$sctest)
                tests <- cbind(rownames(tests), tests)
                colnames(tests) <- c("Statistical test", "Value",
                                     "Degrees of freedom", "p-value")
                return(tests)
            }, options = list(scrollX = TRUE))
            
            # output$coxPlot <- renderHighchart({
            #     # Plot survival curves
            #     hchart(surv, ranges = intRanges) %>%
            #         hc_chart(zoomType="xy") %>%
            #         hc_yAxis(title=list(text="Proportion of individuals")) %>%
            #         hc_xAxis(title=list(text="Time in days"))
            # })
        }
    })
}

attr(survivalUI, "loader") <- "plots"
attr(survivalUI, "name") <- "Survival curves"
attr(survivalServer, "loader") <- "plots"