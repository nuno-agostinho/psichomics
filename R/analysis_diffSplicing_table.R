#' @rdname appUI
#' 
#' @importFrom shinyjs hidden disabled
#' @importFrom shiny actionLink downloadLink selectizeInput uiOutput tags
#' actionButton checkboxGroupInput helpText tagList sidebarLayout mainPanel
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom DT dataTableOutput
#' @importFrom highcharter highchartOutput
diffSplicingTableUI <- function(id) {
    ns <- NS(id)
    pvalueAdjust <- list("No p-value adjustment"="none",
                         "False Discovery Rate"=c(
                             "Benjamini-Hochberg's method"="BH",
                             "Benjamini-Yekutieli's method"="BY"),
                         "Family-wise Error Rate"=c(
                             "Bonferroni correction"="bonferroni",
                             "Holm's method"="holm",
                             "Hochberg's method"="hochberg",
                             "Hommel's method"="hommel"))
    
    statAnalysesOptions <- div(
        id=ns("statAnalysesOptions"),
        selectGroupsUI(ns("diffGroups"), label="Groups of samples to analyse",
                       noGroupsLabel="All samples as one group",
                       groupsLabel="Samples by selected groups"),
        selectGroupsUI(ns("diffASevents"), label="Splicing events to analyse",
                       noGroupsLabel="All splicing events",
                       groupsLabel="Splicing events from selected groups"),
        checkboxGroupInput(
            ns("statsChoices"), width="100%",
            "Choose statistical analyses to perform:",
            # Basic stats is on and disabled by JavaScript
            c("Variance and median"="basicStats",
              "Unpaired t-test (2 groups)"="ttest",
              "Wilcoxon rank sum test (2 groups)"="wilcoxRankSum",
              "Kruskal-Wallis rank sum test (2 or more groups)"="kruskal", 
              "Levene's test (2 or more groups)"="levene",
              "Fligner-Killeen test (2 or more groups)"="fligner",
              "Distribution of alternative splicing quantification per group"=
                  "density"),
            selected=c("basicStats", "kruskal", "levene", "density", "ttest",
                       "fligner", "wilcoxRankSum")),
        # Disable checkbox of basic statistics and of PSI distribution
        tags$script('$("[value=basicStats]").attr("disabled", true);'),
        tags$script('$("[value=density]").attr("disabled", true);'),
        helpText("For each alternative splicing event, groups with one or less",
                 "non-missing values are discarded."), hr(),
        selectizeInput(ns("pvalueAdjust"), selected="BH", width="100%",
                       "P-value adjustment", pvalueAdjust),
        processButton(ns("startAnalyses"), "Perform analyses"))
    
    labelsPanel <- tabPanel(
        "Labels",
        bsCollapse(
            bsCollapsePanel(
                "Label top differentially spliced events", value="top", 
                checkboxInput(
                    ns("labelTopEnable"), width="100%",
                    "Enable labelling of top differentially spliced events"),
                div(id=ns("labelTopOptions"),
                    selectizeInput(
                        ns("labelSortBy"), choices=NULL, width="100%",
                        "Sort top differentially spliced events by"),
                    radioButtons(ns("labelOrder"), "Sorting order",
                                 choices=c("Decreasing order"="decreasing",
                                           "Increasing order"="increasing")),
                    sliderInput(
                        ns("labelTop"), value=10, min=1, max=1000, 
                        width="100%", "Number of top events to label"))),
            bsCollapsePanel(
                "Label selected alternative splicing events", value="events",
                checkboxInput(
                    ns("labelEventEnable"), width="100%",
                    "Enable labelling of selected alternative splicing events"),
                selectizeInput(
                    width="100%", multiple=TRUE, choices=c(
                        "Type to search for alternative splicing events..."=""),
                    ns("labelEvents"), "Alternative splicing events to label")),
            bsCollapsePanel(
                "Label alternative splicing events from selected genes",
                value="genes", checkboxInput(
                    ns("labelGeneEnable"), width="100%",
                    "Enable labelling of events from selected genes"),
                selectizeInput(ns("labelGenes"), "Genes to label", width="100%",
                               choices=c("Type to search for a gene..."=""),
                               multiple=TRUE))),
        checkboxInput(ns("labelOnlyGene"), value=TRUE, width="100%",
                      "Label points using the gene symbol only"),
        actionButton(ns("unlabelPoints"), "Remove labels"),
        processButton(ns("labelPoints"), "Label points"))
    eventOptions <- prepareEventPlotOptions(ns("eventOptions"), ns, labelsPanel)
    
    survivalOptions <- tagList(
        helpText("For each splicing event, find the PSI cutoff that maximizes",
                 "differences in survival."),
        radioButtons(ns("censoring"), "Data censoring", width="100%",
                     selected="right", inline=TRUE, choices=c(
                         Left="left", Right="right",
                         Interval="interval", "Interval 2"="interval2")),
        selectizeInput(ns("timeStart"), choices=character(), "Follow up time",
                       width="100%"),
        # If the chosen censoring contains the word 'interval', show this input
        conditionalPanel(
            sprintf("input[id='%s'].indexOf('interval') > -1", ns("censoring")),
            selectizeInput(ns("timeStop"), choices=character(), "Ending time",
                           width="100%")),
        helpText("For patients for which there is no event reported, time",
                 "to last follow up is used instead."),
        selectizeInput(ns("event"), choices=NULL, width="100%",
                       "Event of interest"),
        selectGroupsUI(
            ns("sampleFiltering"),
            label=div(id=ns("helpFiltering"), "Sample filtering", 
                      icon("question-circle"))),
        bsTooltip(ns("helpFiltering"), options=list(container="body"),
                  placement="right", patientMultiMatchWarning()),
        radioButtons(
            ns("selected"), "Perform survival analysis based on", width="100%",
            choices=c(
                "Splicing events in current page of the table"="shown",
                "All splicing events in the table (slower)"="filtered",
                "All splicing events (slowest)"="all")),
        processButton(ns("survival"), "Plot survival curves")
    )
    
    sidebar <- sidebar(
        bsCollapse(
            id=ns("diffSplicingCollapse"), open="statAnalyses",
            bsCollapsePanel(
                list(icon("cogs"), "Perform statistical analyses"),
                value="statAnalyses", style="info",
                errorDialog(
                    paste("Alternative splicing quantification is required for",
                          "differential splicing analysis."),
                    id=ns("missingIncLevels"),
                    buttonLabel="Alternative splicing quantification",
                    buttonIcon="calculator",
                    buttonId=ns("loadIncLevels")),
                hidden(statAnalysesOptions)),
            bsCollapsePanel(
                list(icon("sliders"), "Plot options and table filtering"),
                style="info", value="plotEvents",
                errorDialog("Differential splicing analysis not yet performed.",
                            id=ns("missingDiffAnalyses")),
                hidden(eventOptions)),
            bsCollapsePanel(
                list(icon("heartbeat"),
                     "Survival analyses by splicing quantification cutoff"),
                style="info", value="survivalOptionsPanel",
                hidden(div(id=ns("survivalOptions"), survivalOptions)),
                errorDialog("Differential splicing analysis not yet performed.",
                            id=ns("survivalOptions-missingDiffSplicing")),
                errorDialog("Clinical data is not loaded.",
                            id=ns("survivalOptions-missingClinicalData"),
                            buttonLabel="Load survival data",
                            buttonId=ns("loadClinical")))))
    
    downloadTable <- div(
        class="btn-group dropup",
        tags$button(class="btn btn-default dropdown-toggle", type="button",
                    "data-toggle"="dropdown", "aria-haspopup"="true", 
                    "aria-expanded"="false", icon("download"), 
                    "Save table", tags$span(class="caret")),
        tags$ul(class="dropdown-menu", 
                tags$li(downloadLink(ns("downloadAll"), "All data")),
                tags$li(downloadLink(ns("downloadSubset"), "Filtered data"))))
    
    groupCreation <- div(
        class="btn-group dropup",
        tags$button(class="btn btn-default dropdown-toggle", type="button",
                    "data-toggle"="dropdown", "aria-haspopup"="true",
                    "aria-expanded"="false", icon("object-group"),
                    "Create group based on...", tags$span(class="caret")),
        tags$ul(class="dropdown-menu",
                disabled(tags$li(id=ns("groupBySelectedContainer"),
                                 actionLink(ns("groupBySelected"),
                                            "Selected splicing events"))),
                tags$li(actionLink(ns("groupByDisplayed"),
                                   "Splicing events displayed in the table"))))
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebar, mainPanel(
                ggplotUI(ns("psi-volcano")),
                dataTableOutput(ns("statsTable")),
                hidden(div(id=ns("tableToolbar"), class="btn-toolbar",
                           role="toolbar", downloadTable, groupCreation)),
                highchartOutput(ns("highchartsSparklines"), 0, 0))))
}

#' Create survival data based on a PSI cutoff
#' 
#' Data is presented in the table for statistical analyses
#' 
#' @inheritParams optimalSurvivalCutoff
#' @param eventPSI Numeric: alternative splicing quantification for multiple
#' samples relative to a single splicing event
#' @inheritParams assignValuePerPatient
#' 
#' @importFrom shiny tags
#' @importFrom jsonlite toJSON
#' @importFrom highcharter hc_title hc_legend hc_xAxis hc_yAxis hc_tooltip
#' hc_chart hc_plotOptions
#' 
#' @return Survival data including optimal PSI cutoff, minimal survival p-value
#' and HTML element required to plot survival curves
createOptimalSurvData <- function(eventPSI, clinical, censoring, event, 
                                  timeStart, timeStop, match, patients, 
                                  samples) {
    # Assign a value to patients based on their respective samples
    eventPSI <- assignValuePerPatient(eventPSI, match, patients=patients,
                                      samples=samples)
    opt <- optimalSurvivalCutoff(clinical, eventPSI, censoring, event, 
                                 timeStart, timeStop)
    
    # Assign a value based on the inclusion levels cutoff
    cutoff <- opt$par
    group  <- labelBasedOnCutoff(eventPSI, cutoff, label="")
    
    survTerms <- processSurvTerms(clinical, censoring, event, timeStart, 
                                  timeStop, group)
    surv <- survfit(survTerms)
    hc <- plotSurvivalCurves(surv, mark=FALSE, auto=FALSE)
    
    # Remove JavaScript used for colouring each series
    for (i in seq(hc$x$hc_opts$series))
        hc$x$hc_opts$series[[i]]$color <- NULL
    
    hc <- as.character(toJSON(hc$x$hc_opts$series, auto_unbox=TRUE))
    
    updateProgress("Survival analysis", console=FALSE)
    return(c("Optimal survival PSI cutoff"=cutoff,
             "Minimal survival p-value"=opt$value,
             "Survival curves"=hc))
}

#' Optimal survival difference given an inclusion level cutoff for a specific
#' alternative splicing event
#' 
#' @importFrom shinyjs runjs show hide
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
optimSurvDiffSet <- function(session, input, output) {
    ns <- session$ns
    
    # Interface of survival analyses
    observe({
        attrs  <- getPatientAttributes()
        diffAn <- getDifferentialSplicing()
        
        if (!is.null(attrs) && !is.null(diffAn)) {
            hide("survivalOptions-missingClinicalData")
            hide("survivalOptions-missingDiffSplicing")
            show("survivalOptions")
            updateClinicalParams(session, attrs)
        } else {
            hide("survivalOptions")
            if (is.null(attrs)) {
                show("survivalOptions-missingClinicalData")
                hide("survivalOptions-missingDiffSplicing")
            } else if (is.null(diffAn)) {
                hide("survivalOptions-missingClinicalData")
                show("survivalOptions-missingDiffSplicing")
            }
        }
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        anyDiffSplicing <- !is.null( getDifferentialSplicing() )
        anyPatients     <- !is.null( getPatientId() )
        anyCensoring    <- !is.null( input$censoring )
        
        if (anyDiffSplicing && anyPatients && anyCensoring) {
            if (grepl("interval", input$censoring, fixed=TRUE)) {
                label <- "Starting time"
            } else {
                label <- "Follow up time"
            }
            updateSelectizeInput(session, "timeStart", label=label)
        }
    })
    
    #' Calculate optimal survival cutoff for the inclusion levels of a given
    #' alternative splicing event
    observeEvent(input$survival, {
        time <- startProcess("survival")
        isolate({
            patients      <- getPatientId()
            psi           <- getInclusionLevels()
            match         <- getClinicalMatchFrom("Inclusion levels")
            statsTable    <- getDifferentialSplicing()
            statsFiltered <- getDifferentialSplicingFiltered()
            optimSurv     <- getDifferentialSplicingSurvival()
            # User input
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            display   <- input$statsTable_rows_current
            filtered  <- input$statsTable_rows_all
            selected  <- input$selected
            samples   <- getSelectedGroups(input, "sampleFiltering", "Samples")
            # Get clinical data for the required attributes
            followup <- "days_to_last_followup"
            clinical <- getClinicalDataForSurvival(timeStart, timeStop, event,
                                                   followup)
        })
        
        if ("shown" %in% selected) {
            if (!is.null(display)) {
                events <- rownames(statsFiltered)[display]
                subset <- psi[events, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it was not possible to get the", 
                           "events shown in the table.",
                           caller="Differential splicing analysis")
                endProcess("survival")
                return(NULL)
            }
        } else if ("filtered" %in% selected) {
            if (!is.null(filtered)) {
                events <- rownames(statsFiltered)[filtered]
                subset <- psi[events, ]
            } else {
                errorModal(
                    session, "Error with selected events",
                    "Unfortunately, it was not possible to get the events from",
                    "the table.", caller="Differential splicing analysis")
                endProcess("survival")
                return(NULL)
            }
        } else if ("all" %in% selected) {
            subset <- psi
        }
        startProgress("Performing survival analysis", nrow(subset))
        opt <- apply(subset, 1, createOptimalSurvData, clinical, censoring, 
                     event, timeStart, timeStop, match, patients, 
                     unlist(samples))
        
        if (length(opt) == 0) {
            errorModal(session, "No survival analyses",
                       "Optimal PSI cutoff for the selected alternative",
                       "splicing events returned no survival analyses.",
                       caller="Differential splicing analysis")
        } else {
            df <- data.frame(t(opt), stringsAsFactors=FALSE)
            if (is.null(optimSurv)) {
                # Prepare survival table
                nas <- rep(NA, nrow(statsTable))
                optimSurv <- data.frame(as.numeric(nas), as.numeric(nas),
                                        as.character(nas),
                                        stringsAsFactors=FALSE)
                rownames(optimSurv) <- rownames(statsTable)
                colnames(optimSurv) <- colnames(df)
            }
            
            optimSurv[rownames(df), 1] <- as.numeric(df[ , 1])
            optimSurv[rownames(df), 2] <- as.numeric(df[ , 2])
            
            # Prepare survival charts
            hc <- highchart() %>%
                hc_title(text=NULL) %>%
                hc_legend(enabled=FALSE) %>%
                hc_xAxis(title=list(text=""), showLastLabel=TRUE, visible=FALSE,
                         crosshair=FALSE) %>%
                hc_yAxis(title=list(text=""), endOnTick=FALSE, crosshair=FALSE,
                         startOnTick=FALSE, visible=FALSE)  %>%
                hc_tooltip(
                    headerFormat=paste(
                        tags$small("{point.x}", scale <- "days"), br(),
                        span(style="color:{point.color}", "\u25CF "),
                        tags$b("{series.name}"), br()),
                    pointFormat=paste(
                        "Survival proportion: {point.y:.3f}", br(),
                        "Records: {series.options.records}", br(),
                        "Events: {series.options.events}", br(),
                        "Median: {series.options.median}")) %>%
                hc_chart(zoomType=NULL, width=120, height=20, 
                         backgroundColor="", margin=c(2, 0, 2, 0), 
                         style=list(overflow='visible')) %>%
                hc_plotOptions(series=list(stickyTracking=FALSE, cursor="non",
                                           animation=FALSE, fillOpacity=0.25,
                                           marker=list(radius=1)))
            data <- as.character(df[ , 3])
            optimSurv[rownames(df), 3] <- createSparklines(
                hc, data, rownames(df), groups=names(samples), "showSurvCutoff")
            setDifferentialSplicingResetPaging(FALSE)
            setDifferentialSplicingSurvival(optimSurv)
        }
        
        # Make survival columns visible
        visibleCols <- sprintf(
            "var table=$(\"#%s table\")
            table.DataTable().columns([6, 7, 8]).visible(true);",
            ns("statsTable"))
        runjs(visibleCols)
        
        # Scroll to survival column
        scroll <- sprintf("var col=$(\"#%s th[aria-label*='Survival']\")
                          $('body').animate({
                          scrollTop: col.offset().top - 50,
                          scrollLeft: col.offset().left - 300
                          }, 1000);", ns("statsTable"))
        runjs(scroll)
        
        endProcess("survival", time)
    })
    }

#' Set of functions to perform differential analyses
#' 
#' @importFrom shinyBS updateCollapse
#' 
#' @inherit diffSplicingTableServer
diffSplicingSet <- function(session, input, output) {
    ns <- session$ns
    
    observe({
        psi <- getInclusionLevels()
        if (is.null(psi)) {
            show("missingIncLevels")
            hide("statAnalysesOptions")
        } else {
            hide("missingIncLevels")
            show("statAnalysesOptions")
        }
    })
    
    performDiffSplicing <- reactive({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        statsChoices <- input$statsChoices
        pvalueAdjust <- input$pvalueAdjust
        
        totalTime <- startProcess("startAnalyses")
        
        # Prepare groups of samples to analyse
        groups <- getSelectedGroups(input, "diffGroups", "Samples",
                                    filter=colnames(psi))
        if ( !is.null(groups) ) {
            colour     <- attr(groups, "Colour")
            attrGroups <- groups
            psi        <- psi[ , unlist(groups), drop=FALSE]
            groups     <- rep(names(groups), sapply(groups, length))
            attr(groups, "Colour") <- colour
        } else {
            attrGroups <- "All samples"
            groups <- rep(attrGroups, ncol(psi))
        }
        
        # Prepare splicing events to analyse
        ASevents <- getSelectedGroups(input, "diffASevents", "ASevents",
                                      filter=rownames(psi))
        if (!is.null(ASevents) ) 
            psi <- psi[unique(unlist(ASevents)), , drop=FALSE]
        
        stats <- diffAnalyses(psi, groups, statsChoices,
                              pvalueAdjust=pvalueAdjust)
        attr(stats, "groups") <- attrGroups
        setDifferentialSplicing(stats)
        setDifferentialSplicingSurvival(NULL)
        updateCollapse(session, "diffSplicingCollapse", "plotEvents")
        endProcess("startAnalyses", totalTime)
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        psi <- isolate(getInclusionLevels())
        diffSplicing <- isolate(getDifferentialSplicing())
        if ( is.null(psi) ) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if ( !is.null(diffSplicing) ) {
            warningModal(session, "Differential splicing already performed",
                         "Do you wish to discard the current results?",
                         footer=actionButton(ns("replace"), "Discard",
                                             class="btn-warning",
                                             "data-dismiss"="modal"),
                         caller="Differential splicing analysis")
        } else {
            performDiffSplicing()
        }
    })
    
    # Replace previously performed differential analyses
    observeEvent(input$replace, {
        performDiffSplicing()
        # Reset previous results from differential analyses
        setDifferentialSplicingFiltered(NULL)
        setZoom("psi-volcano", NULL)
        setSelectedPoints("psi-volcano", NULL)
        setHighlightedPoints("psi-volcano", NULL)
        setDifferentialSplicingSurvival(NULL)
        setLabelledPoints("psi-volcano", NULL)
    })
}

#' @rdname appServer
diffSplicingTableServer <- function(input, output, session) {
    selectGroupsServer(session, "diffGroups", "Samples")
    selectGroupsServer(session, "diffASevents", "ASevents")
    selectGroupsServer(session, "sampleFiltering", "Samples",
                       # Prefer TCGA tumour samples
                       preference="Primary solid Tumor")
    
    observeEvent(input$loadClinical, 
                 missingDataGuide("Clinical data"))
    observeEvent(input$loadIncLevels, 
                 missingDataGuide("Inclusion levels"))
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    diffSplicingSet(session, input, output)
    analysesPlotSet(
        session, input, output, "PSI", "psi-volcano", getDifferentialSplicing,
        getDifferentialSplicingFiltered, getDifferentialSplicingSurvival)
    analysesTableSet(
        session, input, output, "PSI", "psi-volcano", getDifferentialSplicing, 
        getDifferentialSplicingFiltered, setDifferentialSplicingFiltered,
        getDifferentialSplicingSurvival, getDifferentialSplicingColumns, 
        setDifferentialSplicingColumns, getDifferentialSplicingResetPaging, 
        setDifferentialSplicingResetPaging)
    
    # Optimal survival difference given an inclusion level cutoff for a 
    # specific alternative splicing event
    optimSurvDiffSet(session, input, output)
}

attr(diffSplicingTableUI, "loader") <- "diffSplicing"
attr(diffSplicingTableUI, "name") <- "Exploratory (multiple splicing events)"
attr(diffSplicingTableUI, "selectEvent") <- FALSE
attr(diffSplicingTableServer, "loader") <- "diffSplicing"