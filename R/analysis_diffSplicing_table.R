#' Interface for differential analyses on all splicing events
#' 
#' @param id Character: identifier
#' 
#' @importFrom shinyjs disabled
#' @importFrom shiny downloadLink selectizeInput uiOutput actionButton tags
#' checkboxGroupInput helpText tagList sidebarLayout mainPanel textOutput
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
        selectizeInput(ns("groupsCol"), choices=NULL,
                       "Clinical groups on which to perform the analyses"), 
        uiOutput(ns("groupsInfo")), hr(),
        checkboxGroupInput(
            ns("statsChoices"),
            "Choose statistical analyses to perform:",
            # Basic stats is on and disabled by JavaScript
            c("Variance and median"="basicStats",
              "Wilcoxon signed rank test (1 group)"="wilcoxSignedRank",
              "Wilcoxon rank sum test (2 groups)"="wilcoxRankSum",
              "Kruskal-Wallis rank sum test (2 or more groups)"="kruskal", 
              "Levene's test (2 or more groups)"="levene",
              "Alternative splicing quantification density"="density"),
            selected=c("basicStats", "kruskal", "levene", "density",
                       "wilcoxSignedRank", "wilcoxRankSum")),
        # Disable checkbox of basic statistics
        tags$script('$("[value=basicStats]").attr("disabled", true);'),
        helpText("For patients for which there is no event reported, time",
                 "to last follow up is used instead."), hr(),
        selectizeInput(ns("pvalueAdjust"), selected="BH",
                       "Adjust p-values of statistical tests", pvalueAdjust),
        disabled(div(id=ns("downloadStats"), class="btn-group",
                     tags$button(class="btn btn-default dropdown-toggle",
                                 type="button", "data-toggle"="dropdown",
                                 "aria-haspopup"="true",
                                 "aria-expanded"="false", 
                                 icon("download"), 
                                 "Download table", tags$span(class="caret")),
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
                textOutput(paste0(ns("eventText"), i)), chart,
                uiOutput(paste0(ns("eventSurvStats"), i)))
    })
    
    tagList(uiOutput(ns("modal")),
            sidebarLayout(sidebar, 
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
optimSurvDiffUI <- function(ns) {
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
        conditionalPanel(paste0(
            "input[id='", ns("censoring"), "'].indexOf('interval') > -1"),
            selectizeInput(ns("timeStop"), choices=NULL, "Ending time")),
        helpText("In case there's no record for a patient, the",
                 "days to last follow up will be used instead."),
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
optimSurvDiff <- function(session, input, output) {
    ns <- session$ns
    
    # Interface of survival analyses
    output$survivalOptions <- renderUI({
        if (!is.null(getDifferentialAnalyses()) && !is.null(getClinicalData()))
            optimSurvDiffUI(ns)
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
            optimSurv <- getDifferentialAnalysesSurvival()
            display   <- input$statsTable_rows_current
            filtered  <- input$statsTable_rows_all
            selected  <- input$selected
        })
        
        if ("shown" %in% selected) {
            if (!is.null(display)) {
                subset <- psi[display, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get events",
                           "shown in the table.")
                endProcess("survival")
                return(NULL)
            }
        } else if ("filtered" %in% selected) {
            if (!is.null(filtered)) {
                subset <- psi[filtered, ]
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
            v <- as.numeric(vector[toupper(names(tumour))])
            
            opt <- suppressWarnings(
                optim(0, testSurvivalCutoff, data=v, filter=tumour,
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
                nas <- rep(NA, nrow(psi))
                optimSurv <- data.frame(nas, nas)
                rownames(optimSurv) <- rownames(psi)
                colnames(optimSurv) <- colnames(df)
            }
            for (col in names(df)) optimSurv[rownames(df), col] <- df[ , col]
            setDifferentialAnalysesSurvival(optimSurv)
        }
        endProcess("survival", time)
    })
    
    # Update groups used for differential splicing analysis
    observe( setDiffSplicingGroups(input$groupsCol) )
    
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
                stat <- "No optimal PSI cut-off found for this event"
            } else {
                splicingEvent  <- rownames(optimSurv)[row]
                stat <- tagList(
                    tags$b("PSI cut-off:"), roundDigits(cutoff), tags$br(),
                    tags$b("p-value:"), pvalue)
            }
            
            if (is.na(splicingEvent)) stat <- NULL
            output[[paste0("eventText", i)]] <- renderText(
                gsub("_", " ", splicingEvent))
            output[[paste0("eventSurvStats", i)]] <- renderUI(stat)
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
            
            # Group samples by the inclusion levels cut-off
            clinical <- getClinicalData()
            clinicalIDs <- nrow(clinical)
            groups <- rep(NA, clinicalIDs)
            
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            psi       <- getInclusionLevels()
            display   <- input$statsTable_rows_current
            filtered  <- input$statsTable_rows_all
            selected  <- input$selected
        })
        
        lapply(seq(10), function(i) {
            row <- input$statsTable_rows_current[i]
            cutoff <- optimSurv[row, 1]
            
            if (!is.null(row) && !is.na(cutoff)) {
                show(paste0("curves", i), anim=TRUE)
                splicingEvent  <- rownames(optimSurv)[row]
                eventPSI <- as.numeric(psi[splicingEvent, 
                                           toupper(names(tumour))])
                
                group  <- rep(NA, nrow(clinical))
                group[tumour] <- eventPSI >= cutoff
                group[group == "TRUE"]  <- paste(">=", cutoff)
                group[group == "FALSE"] <- paste("<", cutoff)
                
                survTerms <- processSurvTerms(
                    clinical, censoring=censoring, event=event,
                    timeStart=timeStart, timeStop=timeStop, group=group)
                surv <- survfit(survTerms$form, data=survTerms$survTime)
                
                hc <- plotSurvivalCurves(surv, mark = FALSE, title=NULL) %>%
                    hc_legend(enabled=FALSE) %>%
                    hc_xAxis(title=list(text="")) %>%
                    hc_yAxis(title=list(text="")) %>%
                    hc_chart(zoomType=NULL)
            } else {
                hc <- NULL
                hide(paste0("curves", i), anim=TRUE)
            }
            
            output[[paste0("curves", i)]] <- renderHighchart(hc)
        })
    })
}

#' Server logic of the exploratory differential analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shinyjs toggleState disable
#' @importFrom DT replaceData dataTableProxy
diffSplicingTableServer <- function(input, output, session) {
    ns <- session$ns
    
    # Information on the data groups from TCGA
    output$groupsInfo <- renderUI({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        
        if (is.null(psi)) return(tagList(
            helpText(icon("exclamation-circle"), 
                     "No alternative splicing quantification loaded.",
                     tags$a(href="#", "Load or calculate it.",
                            onclick=loadRequiredData("Inclusion levels")))))
        
        # Separate samples by their type
        ids <- names(psi)
        type <- parseSampleGroups(ids)
        
        bullet <- "\u2022"
        groups <- NULL
        for (each in unique(type))
            groups <- tagList(groups, br(), bullet, each)
        
        return(tagList(
            helpText("Sample types available in loaded data:", groups)))
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    performStatsAnalyses <- reactive({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        col <- input$groupsCol
        statsChoices <- input$statsChoices
        pvalueAdjust <- input$pvalueAdjust
        
        totalTime <- startProcess("startAnalyses")
        if (col == "Sample types") {
            # Separate samples by their groups
            ids <- names(psi)
            groups <- parseSampleGroups(ids)
        } else {
            # Get groups from column of interest
            clinical <- getClinicalData()
            col <- clinical[[col]]
            
            # Match groups from patients with respective samples
            matches <- getClinicalMatchFrom("Inclusion levels")
            groups <- rep(NA, ncol(psi))
            names(groups) <- colnames(psi)
            samples <- toupper(names(matches))
            groups[samples] <- as.character(col[matches])
            
            # Remove samples with no groups
            nasGroups <- !is.na(groups)
            psi       <- psi[nasGroups]
            groups    <- groups[nasGroups]
        }
        
        stats <- statsAnalyses(psi, groups, statsChoices, 
                               pvalueAdjust=pvalueAdjust,
                               progress=updateProgress)
        attr(stats, "groups") <- getDiffSplicingGroups()
        setDifferentialAnalyses(stats)
        endProcess("startAnalyses", totalTime)
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            col <- input$groupsCol
            statsChoices <- input$statsChoices
            diffSplicing <- getDifferentialAnalyses()
        })
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if (is.null(col)) {
            errorModal(session, "Select groups",
                       "The groups on which to perform statistical analysis",
                       "cannot be empty.")
        } else if (!is.null(diffSplicing)) {
            warningModal(session, "Differential analyses already performed",
                         "Do you wish to replace the loaded analyses?",
                         footer=actionButton(ns("replace"), "Replace",
                                             class="btn-warning",
                                             "data-dismiss"="modal"))
        } else {
            performStatsAnalyses()
        }
    })
    
    observeEvent(input$replace, performStatsAnalyses())
    
    output$statsTable <- renderDataTableSparklines({
        stats <- getDifferentialAnalyses()
        if (!is.null(stats)) {
            # Remove columns of no interest
            return(stats[, !grepl("method|data.name", colnames(stats))])
        }
    }, style="bootstrap", selection="none", filter='top', server=TRUE,
    extensions="Buttons", options=list(
        pageLength=10, rowCallback=JS("createDiffSplicingLinks"), 
        dom='Bfrtip', buttons=I('colvis'),
        columnDefs=list(list(targets=1, searchable=FALSE))))
    
    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategory(), "diff. splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- NULL
            if (!is.null(stats)) densityCol <- match("Density", colnames(stats))
            
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
            if (!is.null(stats)) densityCol <- match("Density", colnames(stats))
            
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
    
    # Update groups columns
    observe({
        clinical <- getClinicalData()
        psi <- getInclusionLevels()
        
        if (!is.null(clinical)) {
            updateSelectizeInput(
                session, "groupsCol", choices=list(
                    "Clinical groups for samples"=c("Sample types"="Sample types"),
                    "Clinical groups for patients"=names(clinical),
                    "d"=c("Start typing to search for clinical groups"="")))
        } else if (!is.null(psi)) {
            updateSelectizeInput(
                session, "groupsCol", choices=list(
                    "Clinical groups for samples"=c("Sample types"="Sample types"),
                    "d"=c("No clinical data loaded"="")))
        } else {
            updateSelectizeInput(session, "groupsCol", 
                                 choices=c("No clinical data loaded"=""))
        }
    })
}

attr(diffSplicingTableUI, "loader") <- "diffSplicing"
attr(diffSplicingTableUI, "name") <- "All events (table)"
attr(diffSplicingTableUI, "selectEvent") <- FALSE
attr(diffSplicingTableServer, "loader") <- "diffSplicing"