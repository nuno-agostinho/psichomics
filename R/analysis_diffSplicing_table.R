#' Interface for differential analyses on all splicing events
#' 
#' @param id Character: identifier
#' 
#' @importFrom shinyjs disabled
#' @importFrom shiny downloadLink selectizeInput uiOutput actionButton tags
#' checkboxGroupInput helpText tagList sidebarLayout mainPanel
#' @importFrom DT dataTableOutput
#' @importFrom highcharter highchartOutput
#' 
#' @return HTML elements
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
    
    sidebar <- sidebarPanel(
        radioButtons(ns("groupsSelect"), 
                     "Clinical groups on which to perform the analyses",
                     choices=c("Sample types"="samples",
                               "Patients' clinical groups"="patients")),
        conditionalPanel(
            sprintf("input[id='%s'] == 'patients'", ns("groupsSelect")),
            selectGroupsUI(ns("diffGroups"), label=NULL)),
        conditionalPanel(
            sprintf("input[id='%s'] == 'samples'", ns("groupsSelect")),
            uiOutput(ns("groupsInfo"))), 
        hr(), checkboxGroupInput(
            ns("statsChoices"),
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
        # Disable checkbox of basic statistics
        tags$script('$("[value=basicStats]").attr("disabled", true);'),
        helpText("For each alternative splicing event, groups with one or less",
                 "non-missing values are discarded."), hr(),
        selectizeInput(ns("pvalueAdjust"), selected="BH",
                       "Adjust p-values of statistical tests", pvalueAdjust),
        disabled(div(id=ns("downloadStats"), class="btn-group",
                     tags$button(class="btn btn-default dropdown-toggle",
                                 type="button", "data-toggle"="dropdown",
                                 "aria-haspopup"="true",
                                 "aria-expanded"="false", 
                                 icon("download"), 
                                 "Save table", tags$span(class="caret")),
                     tags$ul(class="dropdown-menu", 
                             tags$li(downloadLink(ns("downloadAll"), 
                                                  "All data")),
                             tags$li(downloadLink(ns("downloadSubset"), 
                                                  "Filtered data"))))),
        processButton(ns("startAnalyses"), "Perform analyses"),
        uiOutput(ns("survivalOptions")))
    
    td <- lapply(seq(10), function(i) {
        chart <- highchartOutput(paste0(ns("curves"), i), height="100px")
        chart[[1]] <- tagAppendAttributes(
            chart[[1]], style="min-width: 150px;")
        
        tags$td(style="word-wrap: break-word;", style="white-space: normal;", 
                style="padding: 5px;",
                uiOutput(paste0(ns("eventText"), i)), chart,
                uiOutput(paste0(ns("eventSurvStats"), i)))
    })
    
    tagList(uiOutput(ns("modal")),
            sidebarLayout(
                sidebar, 
                mainPanel(
                    dataTableOutput(ns("statsTable")),
                    highchartOutput(ns("densitySparklines"), 0, 0),
                    hidden(
                        tags$table(id=ns("survTable"), class="table", 
                                   class="table-bordered",
                                   class="text-center",
                                   style="margin-top: 20px;",
                                   tags$tbody(
                                       do.call(tags$tr, td[1:5]),
                                       do.call(tags$tr, td[6:10])))))))
}

#' Interface for calculating optimal cut-off and p-value for survival curves
#' differences
#' @param ns Namespace function
#' @return HTML elements to calculate optimal survival difference
optimSurvDiffOptions <- function(ns) {
    tagList(
        hr(), h3("Survival analyses by splicing quantification cut-off"),
        helpText("For each splicing event, find the PSI cut-off that maximizes",
                 "differences in survival."),
        radioButtons(ns("censoring"), "Data censoring", selected="right",
                     inline=TRUE, choices=c(Left="left", Right="right",
                                            Interval="interval", 
                                            "Interval 2"="interval2")),
        selectizeInput(ns("timeStart"), choices = NULL, "Follow up time"),
        # If the chosen censoring contains the word 'interval', show this input
        conditionalPanel(
            sprintf("input[id='%s'].indexOf('interval') > -1", ns("censoring")),
            selectizeInput(ns("timeStop"), choices=NULL, "Ending time")),
        helpText("For patients for which there is no event reported, time",
                 "to last follow up is used instead."),
        selectizeInput(ns("event"), choices = NULL, 
                       "Event of interest"),
        radioButtons(
            ns("selected"), "Perform survival analysis based on",
            choices=c(
                "Splicing events shown in the screen"="shown",
                "Filtered splicing events (may be a slow process)"="filtered",
                "All splicing events (slow process)"="all")),
        processButton(ns("survival"), "Plot survival curves")
    )
}

#' Optimal survival difference given an inclusion level cut-off for a specific
#' alternative splicing event
#' 
#' @importFrom shinyjs runjs
#' @importFrom shiny renderText
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' @return NULL (this function is used to modify the Shiny session's state)
optimSurvDiff <- function(session, input, output) {
    ns <- session$ns
    
    # Interface of survival analyses
    output$survivalOptions <- renderUI({
        if (!is.null(getDifferentialAnalyses()) && !is.null(getClinicalData()))
            optimSurvDiffOptions(ns)
    })
    
    # Update clinical parameters
    observe({
        if (!is.null(getDifferentialAnalyses()) && !is.null(getClinicalData()))
            updateClinicalParams(session)
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        if (!is.null(getDifferentialAnalyses()) &&
            !is.null(getClinicalData()) && !is.null(input$censoring)) {
            label <- "Follow up time"
            if (grepl("interval", input$censoring, fixed=TRUE))
                label <- "Starting time"
            updateSelectizeInput(session, "timeStart", label=label)
        }
    })
    
    #' Calculate optimal survival cut-off for the inclusion levels of a given
    #' alternative splicing event
    observeEvent(input$survival, {
        time <- startProcess("survival")
        isolate({
            # Get tumour sample IDs (normal and control samples are not
            # interesting for survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            # Group samples by the inclusion levels cut-off
            clinical <- getClinicalData()
            clinicalIDs <- nrow(clinical)
            groups <- rep(NA, clinicalIDs)
            
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            psi       <- getInclusionLevels()
            statsTable <- getDifferentialAnalyses()
            optimSurv <- getDifferentialAnalysesSurvival()
            display   <- input$statsTable_rows_current
            filtered  <- input$statsTable_rows_all
            selected  <- input$selected
        })
        
        if ("shown" %in% selected) {
            if (!is.null(display)) {
                events <- rownames(statsTable)[display]
                subset <- psi[events, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get events",
                           "shown in the table.")
                endProcess("survival")
                return(NULL)
            }
        } else if ("filtered" %in% selected) {
            if (!is.null(filtered)) {
                events <- rownames(statsTable)[filtered]
                subset <- psi[events, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get the events",
                           "from the table.")
                endProcess("survival")
                return(NULL)
            }
        } else if ("all" %in% selected) {
            subset <- psi
        }
        startProgress("Performing survival analysis", nrow(subset))
        
        opt <- apply(subset, 1, function(vector) {
            # Retrieve numeric PSIs from tumour samples
            eventPSI <- rep(NA, nrow(clinical))
            eventPSI[tumour] <- as.numeric(vector[toupper(names(tumour))])
            
            opt <- suppressWarnings(
                optim(0, testSurvivalCutoff, data=eventPSI, filter=tumour, 
                      group=groups, clinical=clinical, censoring=censoring,
                      timeStart=timeStart, timeStop=timeStop, event=event,
                      # Method and parameters interval
                      method="Brent", lower=0, upper=1))
            
            updateProgress("Survival analysis", console=FALSE)
            return(c("Optimal survival PSI cut-off"=opt$par,
                     "Minimal survival p-value"=opt$value))
        })
        
        if (length(opt) == 0) {
            errorModal(session, "No survival analyses",
                       "Optimal PSI cut-off for the selected alternative",
                       "splicing events returned no survival analyses.")
        } else {
            df <- data.frame(t(opt))
            if (is.null(optimSurv)) {
                # Prepare survival table
                nas <- rep(NA, nrow(statsTable))
                optimSurv <- data.frame(nas, nas)
                rownames(optimSurv) <- rownames(statsTable)
                colnames(optimSurv) <- colnames(df)
            }
            for (col in names(df)) optimSurv[rownames(df), col] <- df[ , col]
            setDifferentialAnalysesSurvival(optimSurv)
        }
        endProcess("survival", time)
    })
    
    # Hide survival curves if there's no data
    observe({
        if (is.null(getDifferentialAnalysesSurvival()))
            hide("survTable", anim=TRUE, animType="fade")
        else
            show("survTable", anim=TRUE, animType="fade")
    })
    
    observe({
        optimSurv <- getDifferentialAnalysesSurvival()
        if (is.null(optimSurv)) {
            lapply(seq(10), function(i) {
                output[[paste0("eventText", i)]] <- renderText(NULL)
                output[[paste0("eventSurvStats", i)]] <- renderUI(NULL)
            })
            return(NULL)
        }
        
        lapply(seq(10), function(i) {
            if (!is.null(row)) {
                row <- input$statsTable_rows_current[i]
                cutoff <- optimSurv[row, 1]
                pvalue <- optimSurv[row, 2]
            }
            
            if (is.null(row)) {
                stat <- NULL
                splicingEvent <- NULL
            } else if (is.na(cutoff)) {
                splicingEvent <- rownames(optimSurv)[row]
                stat <- div(class="panel panel-default", div(
                    class="panel-heading",
                    "Survival analysis not yet requested for this event"))
            } else if (is.na(pvalue)) {
                splicingEvent <- rownames(optimSurv)[row]
                stat <- tagList(
                    "No optimal PSI cut-off", tags$br(), "found for this event")
            } else {
                splicingEvent  <- rownames(optimSurv)[row]
                stat <- tags$table(class="table", class="table-condensed", 
                                   style="margin-bottom: 0",
                                   tags$tbody(
                                       tags$tr(tags$td(tags$b("PSI cut-off")),
                                               tags$td(roundDigits(cutoff))),
                                       tags$tr(tags$td(tags$b("p-value")),
                                               tags$td(pvalue))))
            }
            
            if (is.null(splicingEvent) || is.na(splicingEvent)) stat <- NULL
            output[[paste0("eventText", i)]] <- renderUI(
                tags$a(
                    gsub("_", " ", splicingEvent),
                    class="label label-default", style="display: inline-block;",
                    style="white-space: normal;", 
                    onclick=paste0("showSurvCutoff('", splicingEvent, "')")))
            output[[paste0("eventSurvStats", i)]] <- renderUI(tags$small(stat))
        })
    })
    
    observe({
        optimSurv <- getDifferentialAnalysesSurvival()
        if (is.null(optimSurv)) {
            lapply(seq(10), function(i) 
                output[[paste0("curves", i)]] <- renderHighchart(NULL))
            return(NULL)
        }
        
        isolate({
            # Get tumour sample IDs (normal and control samples are not
            # interesting for survival analysis)
            match <- getClinicalMatchFrom("Inclusion levels")
            types <- parseSampleGroups(names(match))
            tumour <- match[!grepl("Normal|Control", types)]
            
            clinical <- getClinicalData()
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            psi       <- getInclusionLevels()
        })
        
        # Interface for the survival curves of 10 splicing events
        lapply(seq(10), function(i) {
            hc <- plotMiniSurvivalCurves(i, input, optimSurv, clinical, tumour,
                                         psi, censoring, event, timeStart, 
                                         timeStop)
            output[[paste0("curves", i)]] <- renderHighchart(hc)
        })
    })
}

#' Perform and plot survival curves
#' 
#' @inheritParams processSurvTerms
#' @param survParams List of parameters to plot survival curves
#' @param i Numeric: index of the survival curves plot of interest
#' @param input Shiny input
#' @param filter Numeric or character: filtered samples
#' @param psi Data frame or matrix: alternative splicing quantification
#' 
#' @importFrom highcharter hc_legend hc_xAxis hc_yAxis hc_chart hc_plotOptions
#' 
#' @return A \code{"highchart"} object to plot
plotMiniSurvivalCurves <- function(i, input, survParams, clinical, filter, psi, 
                                   censoring, event, timeStart, timeStop) {
    row <- input$statsTable_rows_current[i]
    cutoff <- survParams[row, 1]
    
    if (!is.null(row) && !is.na(cutoff)) {
        show(paste0("curves", i), anim=TRUE)
        splicingEvent  <- rownames(survParams)[row]
        
        # Retrieve numeric PSIs from filtered samples
        eventPSI <- rep(NA, nrow(clinical))
        eventPSI[filter] <- as.numeric(
            psi[splicingEvent, toupper(names(filter))])
        group <- eventPSI >= cutoff
        
        group[group == "TRUE"]  <- paste(">=", cutoff)
        group[group == "FALSE"] <- paste("<", cutoff)
        
        survTerms <- processSurvTerms(clinical, censoring, event, timeStart, 
                                      timeStop, group)
        surv <- survfit(survTerms)
        
        hc <- plotSurvivalCurves(surv, mark = FALSE, title=NULL) %>%
            hc_legend(enabled=FALSE) %>%
            hc_xAxis(title=list(text=""), showLastLabel=TRUE) %>%
            hc_yAxis(title=list(text="")) %>%
            hc_chart(zoomType=NULL) %>%
            hc_chart(events=list(click=JS(paste0(
                "function(e) { showSurvCutoff('", splicingEvent, "') }"))))
    } else {
        hc <- NULL
        hide(paste0("curves", i), anim=TRUE)
    }
    return(hc)
}

#' Server logic of the exploratory differential analyses
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny checkboxGroupInput
#' @importFrom shinyjs toggleState disable
#' @importFrom DT replaceData dataTableProxy
#' @importFrom utils write.table
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
diffSplicingTableServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "diffGroups", "Clinical data")
    
    # Information on the data groups from TCGA
    output$groupsInfo <- renderUI({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        
        if (is.null(psi)) return(tagList(
            helpText(icon("exclamation-circle"), 
                     "No alternative splicing quantification loaded.",
                     tags$a("Load or calculate it.",
                            onclick=loadRequiredData("Inclusion levels")))))
        
        # Separate samples by their type
        ids <- names(psi)
        type <- parseSampleGroups(ids)
        
        groups <- unique(type)
        checkboxGroupInput(ns("sampleTypes"),
                           "Choose sample types for comparison:",
                           groups, selected=groups)
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    performDiffAnalyses <- reactive({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        select <- input$groupsSelect
        clinicalGroups <- input$diffGroups
        statsChoices <- input$statsChoices
        pvalueAdjust <- input$pvalueAdjust
        col <- getDiffSplicingGroups()
        if (is.null(col) || col=="") return(NULL)
        
        totalTime <- startProcess("startAnalyses")
        groups <- prepareGroupsDiffSplicing(psi, col)
        psi <- groups$psi
        groups <- groups$groups
        
        stats <- diffAnalyses(psi, groups, statsChoices, 
                              pvalueAdjust=pvalueAdjust,
                              progress=updateProgress)
        attr(stats, "groups") <- getDiffSplicingGroups()
        setDifferentialAnalyses(stats)
        setDifferentialAnalysesSurvival(NULL)
        endProcess("startAnalyses", totalTime)
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            select <- input$groupsSelect
            clinicalGroups <- input$diffGroups
            statsChoices <- input$statsChoices
            diffSplicing <- getDifferentialAnalyses()
        })
        
        emptyGroupsError <- function(session) {
            errorModal(session, "Select groups",
                       "The groups on which to perform statistical analysis",
                       "cannot be empty.")
        }
        
        carryOn <- FALSE
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if (select == "patients") {
            if (is.null(clinicalGroups)) {
                emptyGroupsError(session)
            } else {
                clinicalGroups <- isolate(
                    getGroupsFrom("Clinical data")[clinicalGroups])
                intersection <- unlist(clinicalGroups)
                names(intersection) <- rep(names(clinicalGroups), 
                                           sapply(clinicalGroups, length))
                dup <- duplicated(intersection) | duplicated(intersection, 
                                                             fromLast=TRUE)
                if (any(dup)) {
                    intersected <- unique(names(intersection[dup]))
                    errorModal(session, "Group intersection", 
                               "Differential analysis cannot be carried with",
                               "samples belonging to two groups. The following",
                               "groups share samples:",
                               tags$kbd(paste(intersected, collapse = ", ")))
                } else {
                    carryOn <- TRUE
                }
            }
        } else if (select == "samples") {
            if (is.null(input$sampleTypes))
                emptyGroupsError(session)
            else
                carryOn <- TRUE
        }
        
        if (!carryOn) {
            return(NULL)
        } else if (!is.null(diffSplicing)) {
            warningModal(session, "Differential analyses already performed",
                         "Do you wish to replace the loaded analyses?",
                         footer=actionButton(ns("replace"), "Replace",
                                             class="btn-warning",
                                             "data-dismiss"="modal"))
        } else {
            performDiffAnalyses()
        }
    })
    
    observeEvent(input$replace, performDiffAnalyses())
    
    output$statsTable <- renderDataTableSparklines({
        stats <- getDifferentialAnalyses()
        if (!is.null(stats)) {
            # Remove columns of no interest
            return(stats[, !grepl("method|data.name", colnames(stats))])
        }
    }, style="bootstrap", selection="none", filter="top", server=TRUE,
    extensions="Buttons", options=list(
        pageLength=10, rowCallback=JS("createDiffSplicingLinks"), 
        dom="Bfrtip", buttons=I("colvis"),
        columnDefs=list(list(targets=5, searchable=FALSE))))
    
    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategory(), "diff. splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- NULL
            if (!is.null(stats)) densityCol <- match("PSI distribution", 
                                                     colnames(stats))
            
            write.table(stats[-densityCol], file, quote=FALSE, sep="\t",
                        row.names=FALSE)
        }
    )
    
    # Download filtered table
    output$downloadSubset <- downloadHandler(
        filename=paste(getCategory(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- NULL
            if (!is.null(stats)) densityCol <- match("PSI distribution", 
                                                     colnames(stats))
            
            write.table(stats[input$statsTable_rows_all, -densityCol], 
                        file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
    
    # Optimal survival difference given an inclusion level cut-off for a 
    # specific alternative splicing event
    optimSurvDiff(session, input, output)
    
    # Disable download button if statistical table is NULL
    observe({
        if (is.null(getDifferentialAnalyses()))
            disable("downloadStats")
        else
            enable("downloadStats")
    })
    
    observe({
        groups <- NULL
        if (!is.null(input$groupsSelect) && input$groupsSelect != "samples") {
            groups <- input$diffGroups
        } else if (!is.null(input$sampleTypes)) {
            groups <- input$sampleTypes
            attr(groups, "samples") <- TRUE
        }
        setDiffSplicingGroups(groups)
    })
}

attr(diffSplicingTableUI, "loader") <- "diffSplicing"
attr(diffSplicingTableUI, "name") <- "All events (table)"
attr(diffSplicingTableUI, "selectEvent") <- FALSE
attr(diffSplicingTableServer, "loader") <- "diffSplicing"