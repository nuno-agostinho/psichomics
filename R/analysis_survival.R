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
        "Clinical groups (simple)"="groups",
        "Clinical groups (including their interactions)"="formula",
        "Gene expression cutoff from the selected gene"="geCutoff",
        "Inclusion levels cutoff from the selected splicing event"="psiCutoff")
    
    survivalTimeOptions <- tagList(
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
        radioButtons(
            ns("scale"), "Display time in", inline=TRUE,
            c(Days="days", Weeks="weeks", Months="months", Years="years")))
    
    survivalGroups <- tagList(
        radioButtons(ns("modelTerms"), selected="groups",
                     div("Select groups for survival analysis", 
                         icon("question-circle")), choices=modelChoices),
        bsTooltip(
            ns("modelTerms"), placement="right", 
            options = list(container = "body"),
            paste(
                "Perform survival analysis using:", tags$br(), "\u2022",
                "User-created", tags$b("clinical groups"), tags$br(), 
                "\u2022 A formula that can test clinical attributes with",
                tags$b("interactions"), tags$br(), 
                "\u2022", tags$b("Gene expression cutoff"), "based on the",
                "selected gene", tags$br(),
                "\u2022", tags$b("Inclusion levels cutoff"), "from the",
                "selected alternative splicing event")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "groups"),
            selectGroupsUI(
                ns("dataGroups"), label=NULL, returnAllDataValue=FALSE,
                returnAllDataLabel="Display data outside selected groups")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "formula"),
            textAreaInput(
                ns("formula"), "Formula with clinical attributes", 
                placeholder="Start typing for suggested clinical attributes"),
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
            sprintf("input[id='%s'].includes('%s')", ns("modelTerms"), 
                    "Cutoff"),
            selectGroupsUI(
                ns("sampleFiltering"),
                label=div(id=ns("helpFiltering"), "Sample filtering", 
                          icon("question-circle"))),
            bsTooltip(ns("helpFiltering"), options=list(container="body"),
                      placement="right", patientMultiMatchWarning())),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "geCutoff"),
            hidden(selectizeInput(
                ns("geneExpr"), "Gene expression", width="100%",
                choices=c("No gene expression available"=""))),
            hidden(div(id=ns("loadingGenes"), class="progress",
                       div(class="progress-bar progress-bar-striped active",
                           role="progressbar", style="width:100%",
                           "Loading available genes"))),
            hidden(selectizeGeneInput(ns("gene"))),
            hidden(sliderInput(ns("geCutoff"), value=0.5, min=0, max=1,
                               step=0.01, round=-2, "Gene expression cutoff")),
            hidden(uiOutput(ns("geInfo"))),
            uiOutput(ns("gePvaluePlot"))),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("modelTerms"), "psiCutoff"),
            hidden(sliderInput(ns("psiCutoff"), value=0.5, min=0, max=1,
                               step=0.01, "Splicing quantification cutoff")),
            uiOutput(ns("pvaluePlot"))))
    
    survival <- div(
        id=ns("survivalOptions"),
        bsCollapse(open=c("survivalTimeOptions", "survivalGroups", "KMoptions"),
                   multiple=TRUE,
                   bsCollapsePanel(
                       tagList(icon("calendar-times-o"),
                               "Selection of time features"),
                       value="survivalTimeOptions", style="info",
                       survivalTimeOptions),
                   bsCollapsePanel(
                       tagList(icon("users"), "Groups for survival analysis"),
                       value="survivalGroups", style="info", survivalGroups),
                   bsCollapsePanel(
                       tagList(icon("sliders"), "Kaplan-Meier plot options"),
                       value="KMoptions", style="info", kaplanMeierOptions)),
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
        match         <- getClinicalMatchFrom("Inclusion levels")
        splicingEvent <- getEvent()
        # Get user input
        timeStart  <- input$timeStart
        timeStop   <- input$timeStop
        event      <- input$event
        censoring  <- input$censoring
        outGroup   <- input$dataGroupsShowAllData
        modelTerms <- input$modelTerms
        formulaStr <- input$formula
        intRanges  <- input$ranges
        markTimes  <- input$markTimes
        scale      <- input$scale
        # Get chosen groups
        chosen  <- getSelectedGroups(input, "dataGroups", "Patients")
        samples <- getSelectedGroups(input, "sampleFiltering", "Samples")
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
    } else if (modelTerms == "geCutoff") {
        isolate({
            geneExpr <- getGeneExpression()[[input$geneExpr]]
            gene     <- input$gene
            geCutoff <- input$geCutoff
        })
        
        if (is.null(geneExpr)) {
            missingDataModal(session, "Gene Expression",
                             ns("missingGeneExpression"))
            return(NULL)
        } else if (is.null(gene) || gene == "") {
            errorModal(session, "No gene selected", "Please select a gene.",
                       caller="Survival analysis")
            return(NULL)
        }
        
        # Assign values to patients based on their samples
        eventGE <- assignValuePerPatient(geneExpr[gene, ], match, 
                                         patients=patients, 
                                         samples=unlist(samples))
        
        # Assign a value based on the inclusion levels cutoff
        groups <- labelBasedOnCutoff(eventGE, geCutoff, "Gene expression")
        formulaStr <- NULL
    } else if (modelTerms == "psiCutoff") {
        isolate({
            psi       <- getInclusionLevels()
            psiCutoff <- input$psiCutoff
        })
        
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        } else if (is.null(splicingEvent) || splicingEvent == "") {
            errorModal(session, "No event selected",
                       "Select an alternative splicing event.",
                       caller="Survival analysis")
            return(NULL)
        }
        
        # Assign values to patients based on their samples
        eventPSI <- assignValuePerPatient(psi[splicingEvent, ], match, 
                                          patients=patients, 
                                          samples=unlist(samples))
        
        # Assign a value based on the inclusion levels cutoff
        groups <- labelBasedOnCutoff(eventPSI, psiCutoff, "Inclusion levels")
        formulaStr <- NULL
    } else if (modelTerms == "formula") {
        if (input$formula == "" || is.null(input$formula)) {
            errorModal(session, "Empty formula", 
                       "Please fill the formula field.", 
                       caller="Survival analysis")
            return(NULL)
        } else {
            groups <- NULL
        }
    }
    
    interval <- grepl("interval", censoring)
    if (event == "") {
        errorModal(session, "No event selected", 
                   "Please select the event of interest.",
                   caller="Survival analysis")
    } else if (timeStart == "") {
        if (!interval) {
            errorModal(session, "No follow up time selected",
                       "Please select follow up time.",
                       caller="Survival analysis")
        } else {
            errorModal(session, "No starting time selected",
                       "Please select starting time.", 
                       caller="Survival analysis")
        }
    } else if (timeStop == "" && interval) {
        errorModal(session, "No ending time selected",
                   "Please select ending time to use interval censoring.",
                   caller="Survival analysis")
    } else {
        survTerms <- processSurvival(session, clinical, censoring, event, 
                                     timeStart, timeStop, groups, formulaStr, 
                                     scale=scale, coxph=coxph)
        attr(survTerms, "Colour") <- attr(groups, "Colour")
        return(survTerms)
    }
}

#' Logic set to perform survival analysis based on gene expression cut-offs
#' 
#' @inheritParams survivalServer
#' 
#' @importFrom shinyjs show hide
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
geneExprSurvSet <- function(session, input, output) {
    # Update available gene expression data choices
    observe({
        geneExpr <- getGeneExpression()
        if (!is.null(geneExpr)) {
            updateSelectizeInput(session, "geneExpr", 
                                 choices=rev(names(geneExpr)))
            show("geneExpr")
        } else {
            hide("geneExpr")
        }
    })
    
    # Update available gene choices depending on gene expression data loaded
    # Reactive avoids updating if the input remains the same
    updateGeneChoices <- reactive({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        genes <- rownames(geneExpr)
        updateSelectizeInput(session, "gene", choices=genes, server=TRUE)
    })
    
    # Update available gene choices depending on gene expression data loaded
    observe({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        if (!is.null(geneExpr) && input$modelTerms == "geCutoff") {
            show("loadingGenes")
            hide("gene")
            
            updateGeneChoices()
            
            hide("loadingGenes")
            show("gene")
            show("geCutoff")
            show("geInfo")
        } else {
            hide("loadingGenes")
            hide("gene")
            hide("geCutoff")
            hide("geInfo")
        }
    })
    
    # # Update selected gene based on currently selected splicing event
    # observe({
    #     geneExpr <- getGeneExpression()[[input$geneExpr]]
    #     event    <- getEvent()
    #     if (isolate(input$modelTerms) == "geCutoff" && !is.null(geneExpr) &&
    #         !is.null(event)) {
    #         
    #         gene <- parseSplicingEvent(event)$gene[[1]][[1]]
    #         gene <- grep(gene, rownames(geneExpr), value=TRUE)[[1]]
    #         updateSelectizeInput(session, "gene", selected=gene)
    #     }
    # })

    # Update gene expression cutoff values based on selected gene
    # Reactive avoids updating if the input remains the same
    updateGEcutoffSlider <- reactive({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        ge <- as.numeric(geneExpr[input$gene, ])
        updateSliderInput(session, "geCutoff", min=roundMinDown(ge, 2), 
                          max=roundMaxUp(ge, 2), value=round(mean(ge), 2))
    })

    # Update gene expression cutoff values based on selected gene
    observeEvent(input$gene, {
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        if (!is.null(geneExpr) && input$gene != "" &&
            input$modelTerms == "geCutoff") {
            updateGEcutoffSlider()
            enable("geCutoff")
        } else {
            disable("geCutoff")
        }
    })
    
    # Update gene information based on selected gene
    observe({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        gene     <- input$gene
        terms    <- input$modelTerms
        
        patients <- getPatientId()
        match    <- getClinicalMatchFrom("Inclusion levels")
        # Get user input
        timeStart <- input$timeStart
        timeStop  <- input$timeStop
        event     <- input$event
        censoring <- input$censoring
        samples   <- getSelectedGroups(input, "sampleFiltering", "Samples")
        # Get clinical data for the required attributes
        followup <- "days_to_last_followup"
        clinical <- getClinicalDataForSurvival(timeStart, timeStop, event,
                                               followup)
        
        ui <- NULL
        if (!is.null(geneExpr) && !is.null(gene) && !identical(gene, "") &&
            terms == "geCutoff" && !is.null(clinical)) {
            
            # Assign gene expression values to patients based on their samples
            eventGE <- assignValuePerPatient(geneExpr[gene, ], match, 
                                             patients=patients, 
                                             samples=unlist(samples))
            
            # Mean gene expression cutoff
            meanGEcutoff <- round(mean(eventGE, na.rm=TRUE), 2)
            label        <- labelBasedOnCutoff(eventGE, meanGEcutoff, 
                                               label="Gene expression")
            survTerms    <- processSurvTerms(clinical, censoring=censoring,
                                             event=event, timeStart=timeStart,
                                             timeStop=timeStop, 
                                             followup=followup, group=label)
            meanGEpvalue <- testSurvival(survTerms)
            
            updateSliderInput(session, "geCutoff", value=meanGEcutoff)
            
            # Optimal gene expression cutoff
            opt       <- optimalSurvivalCutoff(clinical, eventGE,
                                               censoring=censoring, 
                                               event=event, timeStart=timeStart, 
                                               timeStop=timeStop, 
                                               session=session)
            
            optimal   <- round(opt$par, 2)
            label     <- labelBasedOnCutoff(eventGE, optimal, 
                                            label="Gene expression")
            survTerms <- processSurvTerms(clinical, censoring=censoring,
                                          event=event, timeStart=timeStart,
                                          timeStop=timeStop, followup=followup,
                                          group=label)
            optPvalue <- testSurvival(survTerms)
            
            df <- data.frame(c("Mean expression", "Optimal cutoff"),
                             c(meanGEcutoff, optimal),
                             paste("p-value:", c(meanGEpvalue, optPvalue)))
            ui <- table2html(df, rownames=FALSE, colnames=FALSE, class="table")
            
            addLinkToUpdateSliderValue <- function(val) {
                val  <- format(val, nsmall=2)
                link <- linkToRunJS(val, sprintf("setGEcutoffSlider(%s)", val))
                gsub(val, link, ui, fixed=TRUE)
            }
            
            ui <- addLinkToUpdateSliderValue(meanGEcutoff)
            ui <- addLinkToUpdateSliderValue(optimal)
        }
        
        output$geInfo <- renderUI(tags$html(ui))
        # output$geInfo <- renderUI(tagList(helpText(meanGEtext),
        #                                   helpText(optGEtext)))
    })

    observe({
        patients <- getPatientId()
        geneExpr <- getGeneExpression()[input$geneExpr]
        gene     <- input$gene

        if (is.null(patients)) {
            hide("geOptions")
            info <- helpText(icon("exclamation-circle"),
                             "Please load clinical data.")
        } else if (is.null(geneExpr)) {
            hide("geOptions")
            info <- helpText(icon("exclamation-circle"),
                             "Please load gene expression data.")
        } else {
            info <- NULL
        }
        
        output$gePvaluePlot <- renderUI(info)
    })
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
    
    selectGroupsServer(session, "dataGroups", "Samples")
    selectGroupsServer(session, "sampleFiltering", "Samples",
                       # Prefer TCGA tumour samples
                       preference="Primary solid Tumor")
    
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
    observeEvent(input$missingGeneExpression, 
                 missingDataGuide("Gene expression"))
    
    # Plot survival curves
    observeEvent(input$survivalCurves, {
        isolate({
            splicingEvent <- getEvent()
            assembly      <- getAssemblyVersion()
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
                       tags$code(surv$message),
                       caller="Survival analysis")
            return(NULL)
        }
        
        pvalue <- testSurvival(survTerms)
        
        if (modelTerms == "psiCutoff") {
            plotTitle <- parseSplicingEvent(splicingEvent, char=TRUE,
                                            pretty=TRUE, extra=assembly)
            sub <- paste0("PSI cutoff: ", psiCutoff,
                          "; Log-rank p-value: ", pvalue)
        } else {
            plotTitle <- "Survival analysis"
            sub <- NULL
        }
        
        # Plot survival curves
        attr(surv, "Colour") <- attr(survTerms, "Colour")
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
                             "Obtained a null Cox model.",
                             caller="Survival analysis")
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
        samples <- getSelectedGroups(input, "sampleFiltering", "Samples")
        
        if (is.null(patients)) {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"), 
                            "Please load clinical data."))
        } else if (is.null(getInclusionLevels())) {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"),
                            "Please load or calculate the quantification of",
                            "alternative splicing events."))
        } else if (is.null(getEvent()) || getEvent() == "") {
            hide("psiCutoff")
            return(helpText(icon("exclamation-circle"), 
                            "Please select an alternative splicing event."))
        } else {
            output$survival  <- renderHighchart(NULL)
            output$coxphUI   <- renderUI(NULL)
            # output$coxTests  <- renderDataTable(NULL)
            # output$coxGroups <- renderDataTable(NULL)
            
            # Assign values to patients based on their samples
            eventPSI <- assignValuePerPatient(psi[splicingEvent, ], match, 
                                              patients=patients,
                                              samples=unlist(samples))
            
            show("psiCutoff")
            slider <- uiOutput(ns("cutoffPvalue"))
            categories <- seq(0, 0.99, 0.01)
            
            survTime <- getAttributesTime(clinical, event, timeStart, timeStop)
            pvalues <- lapply(
                categories, testSurvivalCutoff, data=eventPSI,
                clinical=clinical, censoring=censoring, timeStart=timeStart, 
                timeStop=timeStop, event=event, survTime=survTime, 
                session=session, survivalInfo=TRUE)
            
            # Automatically set minimal p-value
            value <- categories[which.min(unlist(pvalues))]
            observe({
                if (is.na(value)) value <- 0.5
                updateSliderInput(session, "psiCutoff", value=value)
            })
            
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
            
            firstSeriesColour <- JS("Highcharts.getOptions().colors[0]")
            label <- tags$label(class="control-label",
                                "-log\u2081\u2080(p-value) plot by cutoff")
            
            # Put the label of p-value plot to the right when there are many
            # significant points to the left
            signif <- pvalues >= -log10(0.05)
            labelAlign <- "left"
            if (sum(signif[1:50]) > sum(signif[51:100])) labelAlign <- "right"
            
            pvaluePlot <- highchart(height="100px") %>%
                hc_add_series(data=data,
                              zones=list(list(value=significance,
                                              color="lightgray"))) %>%
                hc_chart(zoomType="x") %>%
                hc_xAxis(tickInterval=0.1, showLastLabel=TRUE, endOnTick=TRUE,
                         min=0, max=1, minorGridLineWidth=0, 
                         gridLineWidth=0) %>%
                hc_yAxis(crosshair=list(color="gray", width=1,
                                        dashStyle="shortdash"),
                         labels=list(enabled=FALSE), gridLineWidth=0,
                         plotLines=list(list(
                             value=-log10(0.05), color=firstSeriesColour,
                             dashStyle="shortdash", width=1,
                             label=list(
                                 align=labelAlign, text="p < 0.05",
                                 style=list(color=firstSeriesColour))))) %>%
                hc_legend(NULL) %>% 
                hc_tooltip(formatter=JS(
                    "function() { return getPvaluePlotTooltip(this); }")) %>%
                hc_plotOptions(series=list(
                    cursor="pointer",
                    point=list(events=list(click=JS(
                        "function () { setPSIcutoffSlider(this.x) }"))),
                    marker=list(radius=2)))
            
            if (!is.na(value) && value < 1) {
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
        samples <- getSelectedGroups(input, "sampleFiltering", "Samples")
        
        if (is.null(splicingEvent) || splicingEvent == "" || 
            is.null(psi) || is.null(patients)) return(NULL)
        
        # Assign values to patients based on their samples
        eventPSI <- assignValuePerPatient(psi[splicingEvent, ], match, 
                                          patients=patients, 
                                          samples=unlist(samples))
        
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
    
    geneExprSurvSet(session, input, output)
}

attr(survivalUI, "loader") <- "analysis"
attr(survivalUI, "name") <- "Survival analysis"
attr(survivalServer, "loader") <- "analysis"