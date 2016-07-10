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
                    selected = c("basicStats", "kruskal", "levene")),
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
        
        group <- unique(type)
        len <- length(group)
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
        
        time <- Sys.time()
        if (len > 1) {
            print("Performing statistical analyses...")
            stats <- apply(psi, 1, function(row, type) {
                # Kruskal-Wallis test
                kruskal <- NULL
                if ("kruskal" %in% statsChoices) {
                    kruskal <- tryCatch(kruskal.test(row, factor(type)),
                                        error=return)
                    if ("error" %in% class(kruskal)) kruskal <- NA
                }
                # Levene's test
                levene <- NULL
                if ("levene" %in% statsChoices) {
                    nas <- is.na(row)
                    levene <- tryCatch(levene.test(row[!nas],
                                                   factor(type[!nas])),
                                       error=return)
                    if ("error" %in% class(levene)) levene <- NA
                }
                group <- split(row, type)
                samples <- lapply(group, function(i) sum(!is.na(i))) # Number of samples
                med <- lapply(group, median, na.rm=TRUE) # Median
                var <- lapply(group, var, na.rm=TRUE) # Variance
                return(c(Samples=samples, Kruskal=kruskal, Levene=levene, 
                         Variance=var, Median=med))
            }, factor(type))
            
            # Convert to data frame
            df <- do.call(rbind, stats)
            df <- df[, !grepl("method|data.name", colnames(df))]
            
            # Convert to numeric
            df2 <- data.matrix(matrix(ncol=ncol(df), nrow=nrow(df)))
            for (i in seq(ncol(df))) df2[, i] <- as.numeric(unlist(df[ , i]))
            rownames(df2) <- rownames(df)
            colnames(df2) <- colnames(df)
            
            # Show data frame with not a single NA
            ## TODO(NunoA): we shouldn't discard rows with a single NA...
            df2 <- df2[rowSums(is.na(df2)) == 0, ]
            deltaVar <- df2[, grepl("Variance", colnames(df))]
            deltaVar <- deltaVar[, 2] - deltaVar[, 1]
            deltaMed <- df2[, grepl("Median", colnames(df))]
            deltaMed <- deltaMed[, 2] - deltaMed[, 1]
            df3 <- cbind(df2, deltaVar, deltaMed)
            df4 <- data.frame(data.matrix(df3))
            stats <- cbind(Event = rownames(df4), df4)
            setDifferentialAnalyses(stats)
        }
        print(Sys.time() - time)
    })
    
    observe({
        stats <- getDifferentialAnalyses()
        if (is.null(stats)) return(NULL)
        
        # Columns to show in statistical table
        output$showColumns <- renderUI({
            tagList(
                hr(),
                selectizeInput(ns("columns"), "Show columns", multiple=TRUE,
                               choices=colnames(stats), 
                               selected=colnames(stats),
                               options=list(plugins=list('remove_button', 
                                                         'drag_drop'))),
                downloadButton(ns("download"), "Download"))
        })
        
        # Render statistical table with the selected columns
        output$statsTable <- renderDataTable({
            cols <- colnames(stats) %in% input$columns
            stats[, cols]
        }, options=list(pageLength=10, scrollX=TRUE))
        
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