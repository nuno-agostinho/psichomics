## TODO(NunoA): plot using boxplots

# The name used for the plot must be unique
plot <- "Differential analysis table"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput(id("statsChoices"),
                               "Choose statistical analyses to be performed:",
                               c("Wilcoxin Test"="wilcox",
                                 "Kruskal-Wallis Rank Sum Test"="kruskal", 
                                 "Levene's test"="levene")),
            actionButton(id("startAnalyses"), "Perform selected tests")
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
                k <- NULL
                if ("kruskal" %in% statsChoices) {
                    k <- tryCatch(kruskal.test(row, factor(type)),
                                  error=return)
                    if ("error" %in% class(k))
                        k <- NA
                }
                # Levene's test
                l <- NULL
                if ("levene" %in% statsChoices) {
                    nas <- is.na(row)
                    l <- tryCatch(lawstat::levene.test(row[!nas],
                                                       factor(type[!nas])),
                                  error=return)
                    if ("error" %in% class(l)) l <- NA
                }
                return(list(unlist(k), unlist(l)))
            }, factor(type))

            assign("stats", stats, .GlobalEnv)
            
            # Convert to data frame
            df <- do.call(rbind, stats)
            df <- lapply(seq(ncol(df)), function(i) do.call(rbind, df[ , i]))
            
            ns <- lapply(df, colnames)
            ns <- as.character(mapply(function(i, e) paste(i, e), 
                                      c("kruskal", "levene"), ns))
            df <- do.call(cbind, df)
            colnames(df) <- ns
            df <- df[, -c(4, 5, 8, 9)]
            
            # Convert to numeric
            df2 <- data.matrix(matrix(ncol=ncol(df), nrow=nrow(df)))
            for (i in seq(ncol(df))) df2[, i] <- as.numeric(unlist(df[ , i]))
            rownames(df2) <- rownames(df)
            colnames(df2) <- colnames(df)
            
            output[[id("statsTable")]] <- renderDataTable(
                cbind(Event = rownames(df), df2),
                options=list(pageLength=10, scrollX=TRUE))
        }
        print(Sys.time() - time)
    })
}