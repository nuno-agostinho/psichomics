#' User interface of survival analysis
#' 
#' @param id Character: namespace identifier
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS tagList uiOutput sidebarPanel radioButtons helpText hr 
#' selectizeInput checkboxInput sliderInput actionButton mainPanel textAreaInput
#' conditionalPanel
#' 
#' @return Character with HTML
survivalUI <- function(id) {
    ns <- NS(id)
    
    kaplanMeierOptions <- tagList(
        checkboxInput(ns("markTimes"), "Show censored observations",
                      value = TRUE),
        checkboxInput(ns("ranges"), "Show interval ranges", value = FALSE)
    )
    
    modelChoices <- c(
        "Clinical groups"="groups",
        "Clinical groups (interaction)"="formula",
        "Inclusion levels cut-off from the selected splicing event"="psiCutoff")
    
    tagList(
        uiOutput(ns("modal")),
        sidebarPanel(
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
                          "cut-off</b> from the selected alternative splicing",
                          "event")),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "groups"),
                selectGroupsUI(ns("dataGroups"), label=NULL),
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
                    tags$kbd("tumor_stage : gender + race"), br(), br(),
                    "Interesting attributes include", tags$b("tumor_stage"), 
                    "to get tumour stages.")),
            conditionalPanel(
                sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "psiCutoff"),
                uiOutput(ns("optimalPsi"))),
            hr(),
            bsCollapse(open="KM options",
                       bsCollapsePanel(tagList(icon("sliders"),
                                               "Kaplan-Meier plot options"),
                                       value="KM options",
                                       kaplanMeierOptions, style="info")),
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

#' Server logic of survival analysis
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
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
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
survivalServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups", "Clinical data")
    
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
            clinical      <- getClinicalData()
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
            dataGroups <- input$dataGroups
            chosen <- getGroupsFrom("Clinical data")[dataGroups]
        })
        
        if (is.null(clinical)) {
            missingDataModal(session, "Clinical data", ns("missingClinical"))
        } else if (modelTerms == "groups") {
            # Assign one group for each clinical patient
            groups <- groupPerPatient(chosen, nrow(clinical), outGroup)
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
            clinicalPSI <- getPSIperPatient(psi, match, clinical)
            eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
            
            # Assign a value based on the inclusion levels cut-off
            groups <- labelBasedOnCutoff(eventPSI, psiCutoff,
                                         "Inclusion levels")
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
        
        survTerms <- processSurvival(session, clinical, censoring, event, 
                                     timeStart, timeStop, groups, formulaStr, 
                                     scale=scale)
        if (is.null(survTerms)) return(NULL)
        surv <- tryCatch(survfit(survTerms), error = return)
        if ("simpleError" %in% class(surv)) {
            errorModal(session, "Formula error",
                       "The following error was raised:", br(),
                       tags$code(surv$message))
            return(NULL)
        }
        
        pvalue <- testSurvival(survTerms)
        
        if (modelTerms == "psiCutoff") {
            plotTitle <- splicingEvent
            sub <- paste0("Splicing quantification cut-off: ", psiCutoff,
                          "; Log-rank p-value: ", pvalue)
        } else {
            plotTitle <- "Survival analysis"
            sub <- NULL
        }
        
        # Plot survival curves
        hc <- plotSurvivalCurves(surv, markTimes, intRanges, pvalue, plotTitle, 
                                 scale) %>%
            export_highcharts()
        if (!is.null(sub)) hc <- hc_subtitle(hc, text=sub)
        output$survival <- renderHighchart(hc)
    })
    
    # Fit Cox Proportional Hazards model
    observeEvent(input$coxModel, {
        isolate({
            clinical      <- getClinicalData()
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
            psiCutoff  <- input$psiCutoff
            scale      <- input$scale
            # Get chosen groups
            dataGroups <- input$dataGroups
            chosen <- getGroupsFrom("Clinical data")[dataGroups]
        })
        
        if (is.null(clinical)) {
            missingDataModal(session, "Clinical data", ns("missingClinical"))
            return(NULL)
        } else if (modelTerms == "groups") {
            # Assign one group for each clinical patient
            groups <- groupPerPatient(chosen, nrow(clinical), outGroup)
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
            clinicalPSI <- getPSIperPatient(psi, match, clinical)
            eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
            
            # Assign a value based on the inclusion levels cut-off
            groups <- labelBasedOnCutoff(eventPSI, psiCutoff,
                                         "Inclusion levels")
            
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
        
        # Calculate survival curves
        survTerms <- processSurvival(session, clinical, censoring, event,
                                     timeStart, timeStop, groups, formulaStr, 
                                     coxph=TRUE, scale=scale)
        if (!is.null(survTerms)) {
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
            }
        }
        
        output$coxphUI <- renderUI({
            if (is.null(survTerms)) return(NULL)
            
            len <- length(summary$na.action)
            tagList(
                hr(), 
                downloadButton(ns("download"), "Download Cox model information",
                               class="pull-right"),
                h3("Cox", tags$abbr("PH", title="Proportional Hazards"), 
                   "model", tags$small(
                       summary$n, " patients, ", summary$nevent, " events",
                       if (len > 0) 
                           paste0(" (", len, " missing values removed)"))),
                tags$b("Concordance: "), roundDigits(summary$concordance[[1]]),
                tags$b("(SE: "), roundDigits(summary$concordance[[2]]),
                tags$b(")"),
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
            if (is.null(survTerms))
                return(NULL)
            else
                return(roundDigits(cox))
        }, style="bootstrap", selection='none', options=list(scrollX=TRUE))
    })
    
    # Calculate optimal inclusion levels
    output$optimalPsi <- renderUI({
        clinical      <- getClinicalData()
        psi           <- getInclusionLevels()
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart     <- input$timeStart
        timeStop      <- input$timeStop
        event         <- input$event
        censoring     <- input$censoring
        
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
            output$survival  <- renderHighchart(NULL)
            output$coxphUI   <- renderUI(NULL)
            # output$coxTests  <- renderDataTable(NULL)
            # output$coxGroups <- renderDataTable(NULL)
            
            # Assign alternative splicing quantification to patients based on
            # their samples
            clinicalPSI <- getPSIperPatient(psi, match, clinical)
            eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
            
            # Calculate optimal alternative splicing quantification cut-off
            opt <- optimalPSIcutoff(clinical, eventPSI, censoring=censoring, 
                                    event=event, timeStart=timeStart, 
                                    timeStop=timeStop, session=session)
            
            observe({
                value <- 0.5
                if (!is.na(opt$value) && opt$value < 1) value <- opt$par
                updateSliderInput(session, "psiCutoff", value=value)
            })
            
            slider <- tagList(
                sliderInput(ns("psiCutoff"), value = 0.5, min=0, max=1,
                            step=0.01, "Splicing quantification cut-off"),
                uiOutput(ns("thisPvalue")))
            
            if (!is.na(opt$value) && opt$value < 1) {
                return(tagList(
                    slider, div(class="alert alert-success",
                                tags$b("Optimal cut-off:"), round(opt$par, 5), 
                                br(), tags$b("Minimal log-rank p-value:"),
                                round(opt$value, 3))))
            } else {
                return(tagList(
                    slider, div(class="alert alert-warning", "No optimal",
                                "cut-off was found for this splicing event.")))
            }
        }
    })
    
    # Update contextual information for selected PSI cut-off
    observeEvent(input$psiCutoff, {
        clinical      <- getClinicalData()
        psi           <- getInclusionLevels()
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart     <- input$timeStart
        timeStop      <- input$timeStop
        event         <- input$event
        censoring     <- input$censoring
        psiCutoff     <- input$psiCutoff
        
        if (is.null(getEvent()) || getEvent() == "" || 
            is.null(getInclusionLevels()) || is.null(clinical)) return(NULL)
        
        # Assign alternative splicing quantification to patients based on their
        # samples
        clinicalPSI <- getPSIperPatient(psi, match, clinical)
        eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
        
        # Assign a value based on the inclusion levels cut-off
        groups <- labelBasedOnCutoff(eventPSI, psiCutoff, "Inclusion levels")
        
        survTerms <- processSurvTerms(clinical, censoring, event, timeStart,
                                      timeStop, groups)
        surv <- survfit(survTerms)
        pvalue <- testSurvival(survTerms)
        
        patients <- NULL
        if (!is.na(pvalue) && pvalue < 1)
            patients <- paste0("(", surv$n[1], " vs ", surv$n[2], " patients)")
        
        output$thisPvalue <- renderUI(
            tagList(
                div(style="text-align:right; font-size:small",
                    tags$b("p-value of selected cut-off:"), round(pvalue, 3),
                    patients),
                tags$br()))
    })
}

attr(survivalUI, "loader") <- "analysis"
attr(survivalUI, "name") <- "Survival analysis"
attr(survivalServer, "loader") <- "analysis"