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
    
    kaplanMeierOptions <- tagList(
        checkboxInput(ns("markTimes"), "Show time marks", value = TRUE),
        checkboxInput(ns("ranges"), "Show interval ranges", value = FALSE)
    )
    
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
            helpText("In case there's no record for a patient, the days to last",
                     "follow up will be used instead."),
            selectizeInput(ns("event"), "Event of interest",
                           choices=c("No clinical data loaded"="")),
            radioButtons(ns("scale"), "Display time in", inline=TRUE,
                         c(Days="days", Weeks="weeks", Months="months", 
                           Years="years")),
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

#' Calculate optimal alternative splicing quantification cut-off to separate
#' survival curves
#'
#' @details \code{timeStop} is only considered if \code{censoring} is either
#' \code{interval} or \code{interval2}
#'
#' @inheritParams processSurvTerms
#' @inheritParams testSurvivalCutoff
#' @param session Shiny session (only used for the visual interface)
#' 
#' @return Optimal alternative splicing quantification cut-off
#' @export
optimalPSIcutoff <- function(clinical, data, filter, censoring, event,
                             timeStart, timeStop=NULL, session=NULL) {
    groups <- rep(NA, nrow(clinical))
    
    # Supress warnings from failed calculations while optimising
    opt <- suppressWarnings(
        optim(0, testSurvivalCutoff, group=groups, data=data, filter=filter,
              clinical=clinical, censoring=censoring, timeStart=timeStart, 
              timeStop=timeStop, event=event, session=session,
              # Method and parameters interval
              method="Brent", lower=0, upper=1))
    return(opt)
}

#' Server logic of survival analysis
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom R.utils capitalize
#' @importFrom shiny renderUI observe updateSelectizeInput observeEvent isolate 
#' br tagList hr div icon updateSliderInput
#' @importFrom stats pchisq optim
#' @importFrom survival survfit survdiff
#' @importFrom highcharter hchart hc_chart hc_yAxis hc_xAxis hc_tooltip
#' hc_subtitle hc_tooltip renderHighchart hc_title hc_plotOptions
#' @importFrom DT dataTableOutput renderDataTable
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
            # Get user input
            clinical   <- getClinicalData()
            timeStart  <- input$timeStart
            timeStop   <- input$timeStop
            event      <- input$event
            censoring  <- input$censoring
            outGroup   <- input$showOutGroup
            modelTerms <- input$modelTerms
            formulaStr <- input$formula
            intRanges  <- input$ranges
            markTimes  <- input$markTimes
            psi        <- getInclusionLevels()
            splicingEvent <- getEvent()
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
            
            # Get tumour sample IDs (matched normal and control samples are not
            # interesting for this survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Retrieve numeric PSIs from tumour samples
            eventPSI <- as.numeric(psi[splicingEvent, toupper(names(tumour))])
            groups   <- rep(NA, nrow(clinical))
            groups[tumour] <- eventPSI >= psiCutoff
            
            # Assign a value based on the inclusion levels cut-off
            # groups[is.na(groups)] <- "NA"
            groups[groups == "TRUE"]  <- paste("Inclusion levels >=", psiCutoff)
            groups[groups == "FALSE"] <- paste("Inclusion levels <", psiCutoff)
            formulaStr <- NULL
        } else {
            groups <- NULL
        }
        
        survTerms <- processSurvival(session, clinical, censoring, event, 
                                     timeStart, timeStop, groups, formulaStr, 
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
                                 scale)
        if (!is.null(sub)) hc <- hc_subtitle(hc, text=sub)
        output$survival <- renderHighchart(hc)
    })
    
    # Fit Cox Proportional Hazards model
    observeEvent(input$coxModel, {
        isolate({
            # Get user input
            clinical   <- getClinicalData()
            timeStart  <- input$timeStart
            timeStop   <- input$timeStop
            event      <- input$event
            censoring  <- input$censoring
            outGroup   <- input$showOutGroup
            modelTerms <- input$modelTerms
            formulaStr <- input$formula
            intRanges  <- input$ranges
            psi        <- getInclusionLevels()
            splicingEvent <- getEvent()
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
            
            # Get tumour sample IDs (matched normal and control samples are not
            # interesting for this survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Retrieve numeric PSIs from tumour samples
            eventPSI <- as.numeric(psi[splicingEvent, toupper(names(tumour))])
            groups   <- rep(NA, nrow(clinical))
            groups[tumour] <- eventPSI >= psiCutoff
            
            # Assign a value based on the inclusion levels cut-off
            # groups[is.na(groups)] <- "NA"
            groups[groups == "TRUE"]  <- paste("Inclusion levels >=", psiCutoff)
            groups[groups == "FALSE"] <- paste("Inclusion levels <", psiCutoff)
            formulaStr <- NULL
        } else {
            groups <- NULL
        }
        
        # Calculate survival curves
        survTerms <- processSurvival(session, clinical, censoring, event,
                                     timeStart, timeStop, groups, formulaStr, 
                                     coxph=TRUE, scale=scale)
        if (!is.null(survTerms)) {
            summary <- summary(survTerms)
            print(summary)
            
            # General statistical tests
            tests <- rbind("Wald test"=summary$waldtest,
                           "Log test"=summary$logtest,
                           "Score (logrank) test"=summary$sctest)
            colnames(tests) <- c("Value", "Degrees of freedom", "p-value")
            
            # Groups statistics
            cox <- cbind(summary$coefficients,
                         summary$conf.int[ , 2:4, drop=FALSE])
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
            if (is.null(survTerms)) return(NULL)
            return(roundDigits(cox))
        }, style="bootstrap", selection='none', options=list(scrollX=TRUE))
    })
    
    # Calculate optimal inclusion levels
    output$optimalPsi <- renderUI({
        # Get user input
        clinical      <- getClinicalData()
        timeStart     <- input$timeStart
        timeStop      <- input$timeStop
        event         <- input$event
        censoring     <- input$censoring
        psi           <- getInclusionLevels()
        splicingEvent <- getEvent()
        
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
            
            # Get tumour sample IDs (matched normal and control samples are not
            # interesting for this survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Retrieve numeric PSIs from tumour samples
            psi <- as.numeric(psi[splicingEvent, toupper(names(tumour))])
            
            # Calculate optimal alternative splicing quantification cut-off
            opt <- optimalPSIcutoff(clinical, data=psi, filter=tumour, 
                                    censoring=censoring, event=event,
                                    timeStart=timeStart, timeStop=timeStop,
                                    session=session)
            slider <- tagList(
                sliderInput(ns("psiCutoff"), value = 0.5, min=0, max=1,
                            step=0.01, paste(
                                "Splicing quantification cut-off for the",
                                "selected splicing event")))
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
                    div(tags$b("Optimal splicing quantification cut-off:"),
                        opt$par, br(), 
                        tags$b("Minimal log-rank p-value:"), opt$value)))
            else
                return(tagList(
                    slider, div(icon("bell-o"), "No optimal cut-off was found",
                                "for this alternative splicing event.")))
        }
    })
}

attr(survivalUI, "loader") <- "analysis"
attr(survivalUI, "name") <- "Survival analysis"
attr(survivalServer, "loader") <- "analysis"