# The name used for the plot must be unique
plot <- "Differential analysis table"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput(id("statsChoices"),
                               "Choose statistical analyses to perform:",
                               c("Wilcoxon Test"="wilcox",
                                 "Kruskal-Wallis Rank Sum Test"="kruskal", 
                                 "Levene's test"="levene"),
                               selected = c("kruskal", "levene")),
            downloadButton(id("download"), "Download"),
            actionButton(id("startAnalyses"), class = "btn-primary", 
                         "Perform selected tests")
        ), mainPanel(
            dataTableOutput(id("statsTable"))
        )
    )
)

server <- function(input, output, session) {
    observeEvent(input[[id("startAnalyses")]], {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            statsChoices <- input[[id("statsChoices")]]
        })
        if (is.null(psi)) {
            errorModal(session, "No inclusion levels data",
                       "Insert inclusion levels data first.")
            return(NULL)
        }
        ids <- names(psi)
        
        # Separate samples by their type
        typeList <- readRDS("data/TCGAsampleType.RDS")
        type <- gsub(".*?-([0-9]{2}).-.*", "\\1", ids, perl = TRUE)
        type <- typeList[type]
        
        group <- unique(type)
        len <- length(group)
        
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
                    levene <- tryCatch(lawstat::levene.test(row[!nas],
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
            
            output[[id("statsTable")]] <- renderDataTable(
                stats, options=list(pageLength=10, scrollX=TRUE))
            
            output[[id("download")]] <- downloadHandler(
                filename = paste(getCategories(),  
                                 "Differential splicing analyses"),
                content = function(file)
                    write.table(stats, file, quote=FALSE, row.names=FALSE, 
                                sep="\t")
            )
        }
        print(Sys.time() - time)
    })
}