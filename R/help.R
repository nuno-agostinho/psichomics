#' @rdname appUI
#' @importFrom shiny textOutput
helpUI <- function(id, tab) {
    ns <- NS(id)

    # if (requireNamespace("parallel", quietly = TRUE)) {
    #     cores <- parallel::detectCores()
    #     coresInput <- tagList(
    #         sliderInput(ns("cores"), h4("Number of cores"), value=1, min=1,
    #                     step=1, max=cores, width="auto", post=" core(s)"),
    #         helpText("A total of", cores, "cores were detected.")
    #     )
    # } else {
    #     coresInput <- numericInput(ns("cores"), h4("Number of cores"), value=1,
    #                                min=1, step=1, width="auto")
    # }

    link <- function(href, ..., target="_blank") {
        a <- tags$a(href=href, target=target, ...)
    }

    linkItem <- function(href, ..., target="_blank") {
        a <- link(href, ..., target=target)
        a <- tags$li(class="list-group-item", a)
        return(a)
    }

    guiLink <- paste0("https://nuno-agostinho.github.io/psichomics/",
                      "articles/GUI_tutorial.html")
    cliLink <- paste0("https://nuno-agostinho.github.io/psichomics/",
                      "articles/CLI_tutorial.html")
    customDataLink <- paste0("https://nuno-agostinho.github.io/psichomics/",
                             "articles/custom_data.html")
    customAnnotLink <- paste0("https://nuno-agostinho.github.io/psichomics/",
                              "articles/AS_events_preparation.html")
    tutorials <- div(
        class="panel", class="panel-default",
        div(class="panel-heading", icon("file-alt"), tags$b("Tutorials")),
        tags$ul(
            class="list-group",
            linkItem(guiLink, "Visual interface tutorial"),
            linkItem(cliLink, "Command-line interface (CLI) tutorial"),
            linkItem(customDataLink, "Loading user-provided data"),
            linkItem(customAnnotLink,
                     "Preparing custom alternative splicing annotations")))

    supportLink  <- "https://support.bioconductor.org/t/psichomics/"
    githubIssues <- "https://github.com/nuno-agostinho/psichomics/issues/new"
    feedback <- div(
        class="panel", class="panel-default",
        div(class="panel-heading", icon("comments"), tags$b("Feedback")),
        tags$ul(class="list-group",
                linkItem(supportLink, "Questions and general support"),
                linkItem(githubIssues, "Suggestions and bug reports")))

    groupSite <- "http://imm.medicina.ulisboa.pt/group/distrans/"
    immSite   <- "http://imm.medicina.ulisboa.pt"

    copyright <- sprintf("psichomics %s, 2015-2024",
                         packageVersion("psichomics"))

    about <- div(
        class="panel", class="panel-default",
        div(class="panel-heading",
            icon("info-circle"), tags$b("About")),
        div(class="panel-body",
            "psichomics is an interactive R package for integrative analyses",
            "of alternative splicing and gene expression from large",
            "transcriptomic datasets including those from",
            link("https://cancergenome.nih.gov",
                 "The Cancer Genome Atlas (TCGA)"),
            "and from the",
            link("https://www.gtexportal.org/home/",
                 "Genotype-Tissue Expression (GTEx)"),
            "project, as well as user-provided data."),
        tags$ul(class="list-group",
                tags$li(class="list-group-item",
                        tags$b("Developer:"),
                        link("mailto:nunodanielagostinho@gmail.com",
                             "Nuno Saraiva-Agostinho", icon("envelope"))),
                tags$li(class="list-group-item",
                        tags$b("Supervisor:"), "Nuno L. Barbosa-Morais"),
                linkItem(groupSite, "Disease Transcriptomics lab"),
                linkItem(immSite, "Instituto de Medicina Molecular"),
                tags$li(class="list-group-item",
                        tags$small(class="help-block", style="margin: 0;",
                                   style="text-align: right;", copyright))))

    credits <- c("Lina Gallego", "Marie Bordone", "Mariana Ferreira",
                 "Teresa Maia", "Carolina Leote", "Juan Carlos Verjan",
                 "Bernardo de Almeida")
    credits <- lapply(credits, tags$li, class="list-group-item")
    acknowledgments <- div(
        class="panel", class="panel-default",
        div(class="panel-heading",
            icon("life-ring"), tags$b("Acknowledgments")),
        div(class="panel-body",
            "This work would not be possible without the support of current",
            "and former members of Nuno Morais lab. Thank you for your help."),
        do.call(tags$ul, c(credits, list(class="list-group"))))

    tab(title="Help", icon="question", linkToArticles(),
        h2("Settings", style="margin-top: 0;"), fluidRow(
            # column(4, coresInput),
            column(4,
                   sliderInput(ns("precision"), h4("Numeric precision"),
                               value=3, min=0, max=10, step=1, width="auto",
                               post=" decimal(s)"),
                   textOutput(ns("precisionExample")),
                   helpText("Only applies to new calculations.")),
            column(4,
                   sliderInput(ns("significant"), h4("Significant digits"),
                               value=3, min=0, max=10, step=1, width="auto",
                               post=" digit(s)"),
                   textOutput(ns("significantExample")),
                   helpText("Only applies to new calculations."))),
        h2("Support"),
        fluidRow(column(4, tutorials, feedback),
                 column(4, about),
                 column(4, acknowledgments)))
}

#' @rdname appServer
#' @importFrom shiny observe renderText
helpServer <- function(input, output, session) {
    # observe(setCores(input$cores))
    observe({
        setPrecision(input$precision)
        output$precisionExample <- renderText(
            paste("Example:",
                  formatC(283.5837243243,
                          digits=getPrecision(), format="f")))
    })

    observe({
        setSignificant(input$significant)
        output$significantExample <- renderText(
            paste("Example:",
                  formatC(5.849371935e-06,
                          getSignificant(), format="g")))
    })
}

attr(helpUI, "loader") <- "app"
attr(helpServer, "loader") <- "app"
