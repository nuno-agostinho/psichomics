#' @rdname appUI
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS tagList uiOutput sidebarPanel radioButtons helpText hr 
#' selectizeInput checkboxInput sliderInput actionButton mainPanel textAreaInput
#' conditionalPanel
#' @importFrom shinyjs hidden
survivalUI <- function(id) {
    ns <- NS(id)
    
    kaplanMeierOptions <- tagList(
        checkboxInput(ns("markTimes"), "Show censored observations",
                      value = TRUE),
        checkboxInput(ns("ranges"), "Show interval ranges", value = FALSE))
    
    modelChoices <- c(
        "No groups"="none",
        "Clinical groups"="groups",
        "Clinical groups (interaction)"="formula",
        "Inclusion levels cutoff from the selected splicing event"="psiCutoff")
    
    survival <- div(
        id=ns("survivalOptions"),
        radioButtons(ns("censoring"), "Data censoring", selected="right",
                     inline=TRUE, choices=c(Left="left", Right="right",
                                            Interval="interval",
                                            "Interval 2" = "interval2")),
        selectizeInput(ns("timeStart"), "Follow up time",
                       choices=c("No clinical data loaded"="")),
        # If chosen censoring contains the word 'interval', ask end time
        conditionalPanel(
            paste0("input[id='", ns("censoring"),
                   "'].indexOf('interval') > -1"),
            selectizeInput(ns("timeStop"), "Ending time",
                           choices=c("No clinical data loaded"=""))),
        helpText("For patients for which there is no event reported, time",
                 "to last follow up is used instead."),
        selectizeInput(ns("event"), "Event of interest",
                       choices=c("No clinical data loaded"="")),
        radioButtons(ns("scale"), "Display time in", inline=TRUE,
                     c(Days="days", Weeks="weeks", Months="months", 
                       Years="years")),
        hr(),
        radioButtons(ns("modelTerms"), selected="groups",
                     div("Select groups for survival analysis",
                         icon("question-circle")),
                     choices=modelChoices),
        bsTooltip(ns("modelTerms"), placement="right", 
                  options = list(container = "body"),
                  paste(
                      "Perform survival analysis using:<br/>\u2022",
                      "User-created <b>clinical groups</b><br/>\u2022",
                      "A formula that can test clinical attributes with",
                      "<b>interactions</b><br/>\u2022 <b>Inclusion levels",
                      "cutoff</b> from the selected alternative splicing",
                      "event")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "groups"),
            selectGroupsUI(ns("dataGroups"), label=NULL),
            checkboxInput(ns("showOutGroup"), 
                          "Show data outside selected groups",
                          value = FALSE)),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "formula"),
            textAreaInput(
                ns("formula"), "Formula with clinical attributes", 
                placeholder="Start typing to suggest clinical attributes"),
            uiOutput(ns("formulaSuggestions")),
            helpText(
                "To analyse a series of attributes, separate each",
                "attribute with a", tags$kbd("+"), ". To analyse", 
                "interactions, use", tags$kbd(":"), " (interactions are",
                "only usable with Cox models). For example, ",
                tags$kbd("tumor_stage : gender + race"), br(), br(),
                "Interesting attributes include", tags$b("tumor_stage"), 
                "to get tumour stages.")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"),
                    "psiCutoff"),
            hidden(sliderInput(ns("psiCutoff"), value=0.5, min=0, max=1,
                               step=0.01, "Splicing quantification cutoff")),
            uiOutput(ns("pvaluePlot"))),
        hr(),
        bsCollapse(open="KM options",
                   bsCollapsePanel(tagList(icon("sliders"),
                                           "Kaplan-Meier plot options"),
                                   value="KM options",
                                   kaplanMeierOptions, style="info")),
        actionButton(ns("coxModel"), "Fit Cox PH model"),
        actionButton(ns("survivalCurves"), class="btn-primary",
                     "Plot survival curves"))
    
    tagList(
        uiOutput(ns("modal")),
        sidebarPanel(
            errorDialog("No clinical data is available.",
                        id=ns("survivalDialog"),
                        buttonId=ns("loadClinical"),
                        buttonLabel="Load clinical data"),
            hidden(survival)
        ),
        mainPanel(
            highchartOutput(ns("survival")),
            uiOutput(ns("coxphUI"))
        )
    )
}

#' Prepare survival terms in case of valid input
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param coxph Boolean: prepare data for Cox models? FALSE by default
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
checkSurvivalInput <- function (session, input, coxph=FALSE) {
    ns <- session$ns
    
    isolate({
        patients      <- getPatientId()
        psi           <- getInclusionLevels()
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart  <- input$timeStart
        timeStop   <- input$timeStop
        event      <- input$event
        censoring  <- input$censoring
        outGroup   <- input$showOutGroup
        modelTerms <- input$modelTerms
        formulaStr <- input$formula
        intRanges  <- input$ranges
        markTimes  <- input$markTimes
        psiCutoff  <- input$psiCutoff
        scale      <- input$scale
        # Get chosen groups
        chosen <- getSelectedGroups(input, "dataGroups")
        # Get clinical data for the required attributes
        followup <- "days_to_last_followup"
        clinical <- getClinicalDataForSurvival(timeStart, timeStop, event,
                                               followup, formulaStr=formulaStr)
    })
    
    if (outGroup)
        outGroupName <- "(Outer data)"
    else
        outGroupName <- NA
    
    if ( is.null(patients) ) {
        missingDataModal(session, "Clinical data", ns("missingClinical"))
        return(NULL)
    } else if (modelTerms == "none") {
        groups <- groupPerElem(NULL, patients, outGroupName)
        formulaStr <- NULL
    } else if (modelTerms == "groups") {
        # Assign one group for each clinical patient
        groups <- groupPerElem(chosen, patients, outGroupName)
        formulaStr <- NULL
    } else if (modelTerms == "psiCutoff") {
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        } else if (is.null(splicingEvent) || splicingEvent == "") {
            errorModal(session, "No event selected",
                       "Select an alternative splicing event.")
            return(NULL)
        }
        
        # Assign alternative splicing quantification to patients based on
        # their samples
        clinicalPSI <- getPSIperPatient(psi, match, patients=patients)
        eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
        
        # Assign a value based on the inclusion levels cutoff
        groups <- labelBasedOnCutoff(eventPSI, psiCutoff, "Inclusion levels")
        formulaStr <- NULL
    } else if (modelTerms == "formula") {
        if (input$formula == "" || is.null(input$formula)) {
            errorModal(session, "Empty formula",
                       "Please, fill the formula field.")
            return(NULL)
        } else {
            groups <- NULL
        }
    }
    
    interval <- grepl("interval", censoring)
    if (event == "") {
        errorModal(session, "Empty field for event",
                   "Please, select the event of interest.")
    } else if (timeStart == "") {
        if (!interval) {
            errorModal(session, "Empty field for follow up time",
                       "Please, select follow up time.")
        } else {
            errorModal(session, "Empty field for starting time",
                       "Please, select starting time.")
        }
    } else if (timeStop == "" && interval) {
        errorModal(session, "Empty field for ending time",
                   "Please, select ending time to use interval censoring.")
    } else {
        survTerms <- processSurvival(session, clinical, censoring, event, 
                                     timeStart, timeStop, groups, formulaStr, 
                                     scale=scale, coxph=coxph)
        return(survTerms)
    }
}

#' @rdname appServer
#' 
#' @importFrom R.utils capitalize
#' @importFrom shiny renderUI observe updateSelectizeInput observeEvent isolate 
#' br tagList hr div icon updateSliderInput uiOutput renderUI
#' @importFrom stats pchisq optim
#' @importFrom survival survdiff
#' @importFrom highcharter hchart hc_chart hc_yAxis hc_xAxis hc_tooltip
#' hc_subtitle hc_tooltip renderHighchart hc_plotOptions
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom utils write.table
#' @importFrom shinyjs show hide
survivalServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups")
    
    observe({
        if ( is.null(getPatientAttributes()) ) {
            show("survivalDialog")
            hide("survivalOptions")
        } else {
            hide("survivalDialog")
            show("survivalOptions")
        }
    })
    observeEvent(input$loadClinical, missingDataGuide("Clinical data"))
    
    # Update available clinical data attributes to use in a formula
    output$formulaSuggestions <- renderUI({
        textSuggestions(ns("formula"), getPatientAttributes())
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        label <- "Follow up time"
        if (grepl("interval", input$censoring, fixed=TRUE))
            label <- "Starting time"
        updateSelectizeInput(session, "timeStart", label=label)
    })
    
    # Update available clinical attributes when the clinical data changes
    observe( updateClinicalParams(session, getPatientAttributes()) )
    
    observeEvent(input$missingClinical, missingDataGuide("Clinical data"))
    observeEvent(input$missingInclusionLevels,
                 missingDataGuide("Inclusion levels"))
    
    # Plot survival curves
    observeEvent(input$survivalCurves, {
        isolate({
            splicingEvent <- getEvent()
            # Get user input
            modelTerms <- input$modelTerms
            intRanges  <- input$ranges
            markTimes  <- input$markTimes
            psiCutoff  <- input$psiCutoff
            scale      <- input$scale
        })
        
        survTerms <- checkSurvivalInput(session, input)
        if (is.null(survTerms)) return(NULL)
        
        surv <- tryCatch(survfit(survTerms), error=return)
        if ("simpleError" %in% class(surv)) {
            errorModal(session, "Formula error",
                       "The following error was raised:", br(),
                       tags$code(surv$message))
            return(NULL)
        }
        
        pvalue <- testSurvival(survTerms)
        
        if (modelTerms == "psiCutoff") {
            plotTitle <- parseSplicingEvent(splicingEvent, char=TRUE)
            sub <- paste0("Splicing quantification cutoff: ", psiCutoff,
                          "; Log-rank p-value: ", pvalue)
        } else {
            plotTitle <- "Survival analysis"
            sub <- NULL
        }
        
        # Plot survival curves
        hc <- plotSurvivalCurves(surv, markTimes, intRanges, pvalue, plotTitle, 
                                 scale) %>% export_highcharts()
        if (!is.null(sub)) hc <- hc_subtitle(hc, text=sub)
        output$survival <- renderHighchart(hc)
    })
    
    # Fit Cox Proportional Hazards model
    observeEvent(input$coxModel, {
        survTerms <- checkSurvivalInput(session, input, coxph=TRUE)
        
        if (!is.null(survTerms)) {
            # Properly set group names
            names(survTerms$coefficients) <- 
                gsub("&gt;", ">", names(survTerms$coefficients), fixed=TRUE)
            names(survTerms$coefficients) <- 
                gsub("&lt;", "<", names(survTerms$coefficients), fixed=TRUE)
            
            summary <- summary(survTerms)
            print(summary)
            
            if (is.null(survTerms$coef)) {
                warningModal(session, "Null Cox model",
                             "Obtained a null Cox model.")
                survTerms <- NULL
            } else {
                # General statistical tests
                tests <- rbind("Wald test"=summary$waldtest,
                               "Log test"=summary$logtest,
                               "Score (logrank) test"=summary$sctest)
                colnames(tests) <- c("Value", "Degrees of freedom", "p-value")
                
                # Groups statistics
                cox <- cbind(summary$coefficients,
                             summary$conf.int[ , 2:4, drop=FALSE])
                
                # Properly set colnames and order them
                colnames(cox) <- c("Coefficient", "Hazard ratio", 
                                   "Standard error", "Wald test", "p-value", 
                                   "1 / Hazard ratio", "Lower 95%", "Upper 95%")
                cox <- cox[ , c(1, 3:2, 7:8, 6, 4:5), drop=FALSE]
            }
        }
        
        output$coxphUI <- renderUI({
            if (is.null(survTerms)) return(NULL)
            
            len <- length(summary$na.action)
            tagList(
                hr(), 
                downloadButton(ns("download"), "Save Cox model information",
                               class="pull-right btn-info"),
                h3("Cox", tags$abbr("PH", title="Proportional Hazards"), 
                   "model", tags$small(
                       summary$n, " patients with ", summary$nevent, " events",
                       if (len > 0)
                           paste0(" (", len, " missing values removed)"))),
                tags$b("Concordance: "), roundDigits(summary$concordance[[1]]),
                tags$b("(standard error: "), 
                roundDigits(summary$concordance[[2]]), tags$b(")"),
                br(), tags$b("R\u00B2: "), roundDigits(summary$rsq[[1]]), 
                tags$b("(max possible: "), roundDigits(summary$rsq[[2]]),
                tags$b(")"), 
                dataTableOutput(ns("coxTests")), hr(),
                dataTableOutput(ns("coxGroups"))
            )
        })
        
        output$download <- downloadHandler(
            filename=paste(getCategory(), "Cox PH model", Sys.Date(), ".txt"),
            content=function(file) { 
                len <- length(summary$na.action)
                info <- paste0(
                    "Cox proportional hazards model\n",
                    summary$n, "patients,", summary$nevent, "events",
                    if (len > 0) paste0(" (", len, " missing values removed)"), 
                    "\nConcordance: ", summary$concordance[[1]],
                    "\nSE: ", summary$concordance[[2]],
                    "\nR\u00B2: ", summary$rsq[[1]],
                    "\nMax possible: " , summary$rsq[[2]], "\n")
                
                write(info, file=file)
                suppressWarnings(
                    write.table(tests, file=file, quote=FALSE, sep="\t", 
                                append=TRUE))
                write("", file=file, append=TRUE)
                suppressWarnings(
                    write.table(cox, file=file, quote=FALSE, sep="\t", 
                                append=TRUE, col.names = NA))
            }
        )
        
        output$coxTests <- renderDataTable({
            if (is.null(survTerms)) return(NULL)
            
            pvalue <- signifDigits(tests[, 3])
            tests[, 1:2] <- roundDigits(tests[ , 1:2, drop=FALSE])
            tests[, 3] <- pvalue
            return(tests)
        }, style="bootstrap", selection='none',
        options=list(info=FALSE, paging=FALSE, searching=FALSE, scrollX=TRUE))
        
        output$coxGroups <- renderDataTable({
            if (is.null(survTerms)) return(NULL)
            
            cox <- roundDigits(cox)
            cox[ , "Coefficient"] <- paste(cox[ , "Coefficient"],
                                           cox[ , "Standard error"], 
                                           sep = " \u00B1 ")
            cox[ , "Lower 95%"] <- sprintf("%s to %s",
                                           cox[ , "Lower 95%"],
                                           cox[ , "Upper 95%"])
            colnames(cox)[4] <- "95% CI"
            cox <- cox[ , -c(2, 5:6), drop=FALSE]
            return(cox)
        }, style="bootstrap", selection='none', options=list(scrollX=TRUE))
    })
    
    # Calculate optimal inclusion levels
    output$pvaluePlot <- renderUI({
        patients      <- getPatientId()
        psi           <- getInclusionLevels()
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart     <- input$timeStart
        timeStop      <- input$timeStop
        event         <- input$event
        censoring     <- input$censoring
        # Get clinical data for the required attributes
        followup <- "days_to_last_followup"
        clinical <- getClinicalDataForSurvival(timeStart, timeStop, event,
                                               followup)
        
        if (is.null(patients)) {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"), 
                            "Please, load clinical data."))
        } else if (is.null(getInclusionLevels())) {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"),
                            "Please, load or calculate the quantification of",
                            "alternative splicing events."))
        } else if (is.null(getEvent()) || getEvent() == "") {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"), 
                            "Please, select an alternative splicing event."))
        } else {
            output$survival  <- renderHighchart(NULL)
            output$coxphUI   <- renderUI(NULL)
            # output$coxTests  <- renderDataTable(NULL)
            # output$coxGroups <- renderDataTable(NULL)
            
            # Assign alternative splicing quantification to patients based on
            # their samples
            clinicalPSI <- getPSIperPatient(psi, match, patients=patients)
            eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
            
            # Calculate optimal alternative splicing quantification cutoff
            opt <- optimalPSIcutoff(clinical, eventPSI, censoring=censoring, 
                                    event=event, timeStart=timeStart, 
                                    timeStop=timeStop, session=session)
            
            observe({
                value <- 0.5
                if (!is.na(opt$value) && opt$value < 1) value <- opt$par
                updateSliderInput(session, "psiCutoff", value=value)
            })
            
            show("psiCutoff")
            slider <- uiOutput(ns("cutoffPvalue"))
            categories <- seq(0, 0.99, 0.01)
            
            survTime <- getAttributesTime(clinical, event, timeStart, timeStop)
            pvalues <- lapply(
                categories, testSurvivalCutoff, data=eventPSI,
                clinical=clinical, censoring=censoring, timeStart=timeStart, 
                timeStop=timeStop, event=event, survTime=survTime, 
                session=session, survivalInfo=TRUE)
            
            patients     <- lapply(pvalues, function(n) attr(n, "info")$n)
            noSeparation <- vapply(patients, length, numeric(1)) == 1
            patients[noSeparation] <- NA
            patients1 <- vapply(patients, "[[", 1, FUN.VALUE = numeric(1))
            patients2 <- NA
            patients2[!noSeparation] <- vapply(patients[!noSeparation], 
                                               "[[", 2, FUN.VALUE = numeric(1))
            
            pvalues      <- -log10(unlist(pvalues))
            significance <- -log10(0.05)
            
            data <- data.frame(x=categories, y=pvalues, 
                               patients1=patients1, patients2=patients2)
            data <- list_parse(data)
            
            label <- tags$label(class="control-label",
                                "-log\u2081\u2080(p-value) plot by cutoff")
            pvaluePlot <- highchart(height="100px") %>%
                hc_add_series(data=data,
                              zones=list(list(value=significance,
                                              color="lightgray"))) %>%
                hc_chart(zoomType="x") %>%
                hc_xAxis(tickInterval=0.1, showLastLabel=TRUE, endOnTick=TRUE,
                         min=0, max=1) %>%
                hc_yAxis(crosshair=list(color="gray", width=1, 
                                        dashStyle="shortdash"),
                         labels=list(enabled=FALSE)) %>%
                hc_legend(NULL) %>% 
                hc_tooltip(formatter=JS(
                    "function() { return getPvaluePlotTooltip(this); }")) %>%
                hc_plotOptions(series=list(
                    cursor="pointer",
                    point=list(events=list(click=JS(
                        "function () { setPSIcutoffSlider(this.x) }"))),
                    marker=list(radius=2)))
            
            if (!is.na(opt$value) && opt$value < 1) {
                return(tagList(slider, label, pvaluePlot))
            } else {
                return(tagList(
                    slider, label, pvaluePlot,
                    div(class="alert alert-warning", "No adequate",
                        "cutoff was found for this splicing event.")))
            }
        }
    })
    
    # Update contextual information for selected PSI cutoff
    observeEvent(input$psiCutoff, {
        patients      <- getPatientId()
        psi           <- getInclusionLevels()
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart     <- input$timeStart
        timeStop      <- input$timeStop
        event         <- input$event
        censoring     <- input$censoring
        psiCutoff     <- input$psiCutoff
        # Get clinical data for the required attributes
        followup <- "days_to_last_followup"
        clinical <- getClinicalDataForSurvival(timeStart, timeStop, event,
                                               followup)
        
        if (is.null(splicingEvent) || splicingEvent == "" || 
            is.null(psi) || is.null(patients)) return(NULL)
        
        # Assign alternative splicing quantification to patients based on their
        # samples
        clinicalPSI <- getPSIperPatient(psi, match, patients=patients)
        eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
        
        pvalue <- testSurvivalCutoff(
            psiCutoff, data=eventPSI, clinical=clinical, censoring=censoring, 
            timeStart=timeStart, timeStop=timeStop, event=event, 
            session=session, survivalInfo = TRUE)
        surv <- attr(pvalue, "info")
        
        patients <- NULL
        if (!is.na(pvalue) && pvalue < 1)
            patients <- paste0("(", surv$n[1], " vs ", surv$n[2], " patients)")
        
        output$cutoffPvalue <- renderUI(
            tagList(div(style="text-align:right; font-size:small",
                        tags$b("p-value of selected cutoff:"), round(pvalue, 3),
                        patients),
                    tags$br()))
    })
}

attr(survivalUI, "loader") <- "analysis"
attr(survivalUI, "name") <- "Survival analysis"
attr(survivalServer, "loader") <- "analysis"