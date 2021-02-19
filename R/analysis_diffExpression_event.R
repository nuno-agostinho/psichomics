#' @rdname appUI
#'
#' @importFrom highcharter highchartOutput
#' @importFrom shiny tagList uiOutput NS sidebarLayout numericInput h3 mainPanel
#' actionButton sidebarPanel
diffExpressionEventUI <- function(id) {
    ns <- NS(id)
    return(diffEventUI(id, ns, psi=FALSE))
}

#' @rdname appServer
#'
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide
diffExpressionEventServer <- function(input, output, session) {
    ns <- session$ns

    selectGroupsServer(session, "diffGroups", "Samples")

    observe({
        geneExpr <- getGeneExpression(input$geneExpr)
        gene <- input$gene
        if (is.null(geneExpr) || is.null(gene) || gene == "") return(NULL)
        data <- as.numeric(geneExpr[gene, ])
        updateCheckboxInput(session, "rug", value=length(data) < 500)
    })

    observeEvent(input$analyse, {
        geneExpr <- getGeneExpression(input$geneExpr)
        if (is.null(geneExpr)) {
            missingDataModal(session, "Gene expression", ns("missingGeneExpr"))
            return(NULL)
        }

        gene <- input$gene
        if (is.null(gene) || gene == "") {
            errorModal(session, "No gene selected", "Please select a gene.",
                       caller="Differential expression analysis")
            return(NULL)
        }

        # Prepare groups of samples to analyse
        groups <- getSelectedGroups(input, "diffGroups", "Samples",
                                    filter=colnames(geneExpr))
        colour <- attr(groups, "Colour")
        if ( !is.null(groups) ) {
            attrGroups <- groups
            geneExpr <- geneExpr[ , unlist(groups), drop=FALSE]
            groups <- rep(names(groups), sapply(groups, length))
        } else {
            attrGroups <- "All samples"
            groups <- rep(attrGroups, ncol(geneExpr))
        }

        # Check if analyses were already performed
        stats <- getDifferentialExpression()
        if (!is.null(stats) && identical(attrGroups, attr(stats, "groups"))) {
            stat  <- stats[gene, ]
            adjustedPvalue <- grep("p-value (", names(stats), fixed=TRUE,
                                   value=TRUE)

            logFC      <- signifDigits(stat$"log2 Fold-Change")
            logFCconf1 <- signifDigits(stat$"conf. int1")
            logFCconf2 <- signifDigits(stat$"conf. int2")
            modTstats  <- signifDigits(stat$"moderated t-statistics")
            modTpvalue <- signifDigits(stat$"p-value")
            modTadjustedPvalue <- signifDigits(stat[[adjustedPvalue]])
            bStats     <- signifDigits(stat$"B-statistics")

            diffStats <- tagList(
                hr(), h4("Differential expression summary"),
                tags$b("log2 Fold-Change: "),
                sprintf("%s (%s to %s)", logFC, logFCconf1, logFCconf2), br(),
                tags$b("Moderated t-statistics: "), modTstats, br(),
                tags$b("p-value: "), modTpvalue, br(),
                tags$b(paste0(adjustedPvalue, ": ")), modTadjustedPvalue, br(),
                tags$b("B-statistics: "), bStats)
        } else {
            stat  <- NULL
            diffStats <- NULL
        }

        # Separate samples by their groups
        eventGE <- as.numeric(geneExpr[gene, ])
        names(eventGE) <- colnames(geneExpr)
        eventGE <- filterGroups(eventGE, groups, 2)
        groups  <- attr(eventGE, "Groups")
        attr(groups, "Colour") <- colour

        plot <- plotDistribution(eventGE, groups, psi=FALSE,
                                 title=paste(gene, "gene expression"),
                                 type=isolate(input$plotType), rug=input$rug)
        output$density <- renderHighchart(plot)

        output$basicStats <- renderUI(basicStats(eventGE, groups))
        output$diffStats  <- renderUI(diffStats)
        output$ttest      <- renderUI(ttest(eventGE, groups, stat))
        output$wilcox     <- renderUI(wilcox(eventGE, groups, stat))
        output$kruskal    <- renderUI(kruskal(eventGE, groups, stat))
        output$levene     <- renderUI(levene(eventGE, groups, stat))
        output$fligner    <- renderUI(fligner(eventGE, groups, stat))
        # output$fisher   <- renderUI(fisher(eventGE, groups))
        # output$spearman <- renderUI(spearman(eventGE, groups))

        show("survivalButton")
        show("singleEventInfo")
    })

    observeEvent(input$missingGeneExpr,
                 missingDataGuide("Gene expression"))
    observeEvent(input$missingGeneExprButton,
                 missingDataGuide("Gene expression"))

    # Update available gene choices depending on gene expression data loaded
    # Reactive avoids updating if the input remains the same
    updateGeneChoices <- reactive({
        geneExpr <- getGeneExpression(input$geneExpr)
        genes <- rownames(geneExpr)
        updateSelectizeInput(session, "gene", choices=genes, server=TRUE)
    })

    observe({
        geneExpr <- getGeneExpression()
        if (is.null(geneExpr)) {
            show("missingGeneExpr")
        } else {
            updateSelectizeInput(session, "geneExpr",
                                 choices=rev(names(geneExpr)))
            hide("missingGeneExpr")
        }
    })

    # Show options if gene expression data is available, update available gene
    # expression data choices and update available genes for selection
    observe({
        geneExpr <- getGeneExpression(input$geneExpr)
        if (is.null(geneExpr)) {
            hide("singleEventOptions")
            hide("survivalButton")
            hide("singleEventInfo")
            # Gene-related tasks
            hide("gene")
        } else {
            show("singleEventOptions")
            # Gene-related tasks
            hide("gene")
            updateGeneChoices()
            show("gene")
        }
    })
}

attr(diffExpressionEventUI, "loader") <- "diffExpression"
attr(diffExpressionEventUI, "name") <- "Individual gene"
attr(diffExpressionEventServer, "loader") <- "diffExpression"
