diffAnalysisTableUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebarPanel(
                uiOutput(ns("groupsInfo")),
                checkboxGroupInput(
                    ns("statsChoices"),
                    "Choose statistical analyses to perform:",
                    # Basic stats is on and disabled by JavaScript
                    c("Variance and median"="basicStats",
                      "Wilcoxon Test (1 or 2 groups)"="wilcox",
                      "Kruskal-Wallis Rank Sum Test (2 or more groups)"="kruskal", 
                      "Levene's test (2 or more groups)"="levene"),
                    selected = c("basicStats", "wilcox", "kruskal", "levene")),
                # Disable checkbox of basic statistics
                HTML("<script>",
                     '$("[value=basicStats]").attr("disabled", true);', 
                     "</script>"),
                actionButton(ns("startAnalyses"), class = "btn-primary", 
                             "Perform analyses")
            ), mainPanel(
                uiOutput(ns("showColumns")),
                dataTableOutput(ns("statsTable"))
            )
        )
    )
}

#' @importFrom lawstat levene.test
#' @importFrom stats kruskal.test median wilcox.test var
diffAnalysisTableServer <- function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            statsChoices <- input$statsChoices
        })
        if (is.null(psi)) {
            errorModal(session, "No AS event quantification",
                       "Insert or quantify alternative splicing events by",
                       "going to the", icon("table"), tags$b("Data"),
                       "tab and opening",  icon("calculator"),
                       tags$b("Quantify alternative splicing events"), ".")
            return(NULL)
        }
        ids <- names(psi)
        
        # Separate samples by their type
        type <- getSampleTypes(ids)
        
        # Information on the data groups from TCGA
        output$groupsInfo <- renderUI({
            bullet  <- "\u2022"
            groups <- NULL
            for (each in unique(type))
                groups <- tagList(groups, br(), bullet, each)
            
            return(tagList(
                helpText("The data contains the following groups:", groups),
                hr()))
        })
        
        # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
        step <- 100 # Avoid updating after analysing each event
        startProgress("Performing statistical analysis", 
                      divisions=1+round(nrow(psi)/step))
        time <- Sys.time()
        
        count <- 0
        stats <- apply(psi, 1, function(row, type) {
            count <<- count + 1
            if (count %% step == 0)
                updateProgress("Performing statistical analysis", console=FALSE)
            
            # Filter groups with less data points than required
            threshold <- 1
            names(row) <- type
            
            #' Filter row by a given threshold
            filterByThreshold <- function(thisType, allTypes, row, threshold) {
                row <- row[thisType == allTypes]
                if ( sum(!is.na(row)) >= threshold )
                    return(row)
            }
            row <- lapply(unique(type), filterByThreshold, type, row, threshold)
            
            row  <- unlist(row)
            type <- names(row)
            row  <- as.numeric(row)
            len  <- length(unique(type))
            
            # Wilcoxon tests
            wilcox <- NULL
            if ("wilcox" %in% statsChoices) {
                if (len == 2) {
                    typeOne <- type == unique(type)[1]
                    wilcox  <- suppressWarnings(wilcox.test(row[typeOne], 
                                                            row[!typeOne]))
                } else if (len == 1) {
                    wilcox <- suppressWarnings(wilcox.test(row))
                }
            }
            
            # Kruskal-Wallis test
            kruskal <- NULL
            if ("kruskal" %in% statsChoices && len >= 2) {
                kruskal <- tryCatch(kruskal.test(row, factor(type)),
                                    error=return)
                if ("error" %in% class(kruskal)) kruskal <- NULL
            }
            
            # Levene's test
            levene <- NULL
            if ("levene" %in% statsChoices && len >= 2) {
                nas <- is.na(row)
                levene <- tryCatch(levene.test(row[!nas], factor(type[!nas])),
                                   error=return)
                if ("error" %in% class(levene)) levene <- NULL
            }
            
            # Variance and median
            group <- split(row, type)
            samples <- lapply(group,
                              function(i) sum(!is.na(i))) # Number of samples
            med <- lapply(group, median, na.rm=TRUE) # Median
            var <- lapply(group, var, na.rm=TRUE) # Variance
            
            row <- c(Samples=samples, Wilcox=wilcox, Kruskal=kruskal, 
                     Levene=levene, Variance=var, Median=med)
            row <- row[!vapply(row, is.null, logical(1))] # Remove NULL
            row <- data.frame(row, stringsAsFactors=FALSE)
            return(row)
        }, factor(type))
        
        print(Sys.time() - time)
        
        browser()
        updateProgress("Performing statistical analysis", console=FALSE)
        # Convert to data frame
        df <- do.call(rbind.fill, stats)
        df <- df[, !grepl("method|data.name", colnames(df))]
        
        # Calculate delta variance and delta median if there are only 2 groups
        deltaVar <- df[, grepl("Variance", colnames(df))]
        if (ncol(deltaVar) == 2) {
            deltaVar <- deltaVar[, 2] - deltaVar[, 1]
            deltaMed <- df[, grepl("Median", colnames(df))]
            deltaMed <- deltaMed[, 2] - deltaMed[, 1]
            df <- cbind(df, deltaVar, deltaMed)
        }
        
        stats <- cbind(Event = rownames(psi), df)
        setDifferentialAnalyses(stats)
        
        # parallel::stopCluster(cl)
        print(Sys.time() - time)
        closeProgress()
    })
    
    # Show table and respective interface when statistical table is available
    observe({
        stats <- getDifferentialAnalyses()
        if (is.null(stats)) return(NULL)
        
        # Columns to show in statistical table
        output$showColumns <- renderUI({
            tagList(
                downloadButton(ns("download"), "Download table"),
                selectizeInput(ns("columns"), "Show columns", multiple=TRUE,
                               choices=colnames(stats), width="auto",
                               selected=colnames(stats),
                               options=list(plugins=list('remove_button', 
                                                         'drag_drop'))),
                hr())
        })
        
        # Render statistical table with the selected columns
        output$statsTable <- renderDataTable(
            stats[, input$columns], options=list(pageLength=10, scrollX=TRUE))
        
        # Prepare table to be downloaded
        output$download <- downloadHandler(
            filename=paste(getCategories(), "Differential splicing analyses"),
            content=function(file)
                write.table(stats, file, quote=FALSE, row.names=FALSE, sep="\t")
        )
    })
}

attr(diffAnalysisTableUI, "loader") <- "analysis"
attr(diffAnalysisTableUI, "name") <- "Differential analysis (exploratory)"
attr(diffAnalysisTableServer, "loader") <- "analysis"