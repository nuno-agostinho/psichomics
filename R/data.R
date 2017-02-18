## TODO(NunoA): should default columns be a perfect match or just a partial
## match? If only a partial match... that would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

#' Get data types available from Firehose
#' 
#' @importFrom R.utils capitalize
#' 
#' @return Named character vector
#' @export
#' 
#' @examples
#' getFirehoseDataTypes()
getFirehoseDataTypes <- function() {
    choices <- list("RNA sequencing"=c(
        paste0(c("junction", "exon"), "_quantification"), "Preprocess",
        paste0("RSEM_", c("isoforms", "genes")),
        paste0(c("junction", "gene", "exon"),
               "_expression"), "genes_normalized"))
    names(choices[[1]]) <- capitalize(gsub("_", " ", choices[[1]]))
    return(choices)
}

#' Create a modal warning the user of already loaded data
#' @param modalId Character: identifier of the modal
#' @param replaceButtonId Character: identifier of the button to replace data
#' @param keepButtonId Character: identifier of the button to append data
#' @param session Shiny session
#' @return HTML elements for a warning modal reminding data is loaded
loadedDataModal <- function(session, modalId, replaceButtonId, keepButtonId) {
    ns <- session$ns
    warningModal(session, "Data already loaded",
                 "Would you like to", tags$b("replace"), "the loaded data or",
                 tags$b("keep"), "both the previous and new data?",
                 footer = tagList(
                     actionButton(ns(keepButtonId), "data-dismiss"="modal",
                                  label="Keep both"),
                     actionButton(ns(replaceButtonId), class="btn-warning",
                                  "data-dismiss"="modal", label="Replace")),
                 modalId=modalId)
}

#' Process dataset names
#' 
#' @details Avoid duplicated names and append the technology used for junction
#' quantification
#' 
#' @param data List of lists of data frames
#' 
#' @return Processed list of lists of data frames
processDatasetNames <- function(data) {
    newData <- data
    # Avoid duplicate names in categories
    names(newData) <- renameDuplicated(names(data), 
                                       names(data)[duplicated(names(data))])
    
    ns <- lapply(newData, names)
    for (each in names(ns)) {
        nse <- names(newData[[each]])
        
        # For junction quantification, add the respective sequencing technology    
        index <- nse == "Junction quantification"
        for (k in seq_along(nse)) {
            if (index[[k]]) {
                file <- attr(newData[[each]][[k]], "filename")
                if (grepl("illuminahiseq", file, fixed=TRUE))
                    names(newData[[each]])[[k]] <- paste(
                        names(newData[[each]])[[k]], "(Illumina HiSeq)")
                else if (grepl("illuminaga", file, fixed=TRUE)) 
                    names(newData[[each]])[[k]] <- paste(
                        names(newData[[each]])[[k]], "(Illumina GA)")
            }
        }
        
        # Avoid duplicate names in datasets from the same category
        nse <- names(newData[[each]])
        names(newData[[each]]) <- renameDuplicated(nse, nse[duplicated(nse)])
    }
    return(newData)
}

#' User interface of the data module
#' @param id Character: identifier
#' @param tab Function to create tab
#' @return HTML elements
dataUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "data", bsCollapsePanel,
                             priority="localDataUI")
    
    tab(title="Data", icon="table",
        sidebarLayout(
            sidebarPanel( do.call(bsCollapse, c(id=ns("accordion"), uiList)) ),
            mainPanel( uiOutput(ns("tablesOrAbout")) ) ))
}

#' Creates a tabPanel template for a datatable with a title and description
#'
#' @param ns Namespace function
#' @param title Character: tab title
#' @param tableId Character: id of the datatable
#' @param description Character: description of the table (optional)
#' @param columns Character: column names of the datatable
#' @param visCols Boolean: visible columns
#' @param data Data frame: dataset of interest
#'
#' @importFrom shinyBS bsTooltip
#' @importFrom DT dataTableOutput
#' @importFrom shiny hr br tabPanel selectizeInput column fluidRow p mainPanel
#' downloadButton
#'
#' @return The HTML code for a tabPanel template
tabDataset <- function(ns, title, tableId, columns, visCols, data,
                       description=NULL) {
    tablename <- ns(paste("table", tableId, sep="-"))
    
    downloadId <- paste(tablename, "download", sep="-")
    download <- downloadButton(downloadId, "Save table",
                               class="pull-right btn-info")
    
    if(!is.null(description)) {
        description <- p(tags$strong("Table description:"), description)
        download <- fluidRow(column(10, description), column(2, download))
    }
    
    # Get class of each column
    colType <- sapply(1:ncol(data), function(i) class(data[[i]]))
    colType[colType == "character"] <- "string"
    
    # Show class of each column
    choices <- columns
    names(choices) <- sprintf("%s (%s class)", columns, colType)
    
    visibleColumns <- selectizeInput(
        paste(tablename, "columns", sep="-"), label="Visible columns", 
        choices=choices, selected=visCols, multiple=TRUE, width="auto", 
        options=list(plugins=list('remove_button', 'drag_drop'), render=I(
            "{ item: function(item, escape) {
                return '<div>' + escape(item.value) + '</div>'; } }")))
    tabPanel(title, br(), download, visibleColumns, hr(),
             dataTableOutput(tablename))
}

#' Render a specific data tab (including data table and related interface)
#' 
#' @param index Integer: index of the data to load
#' @param data Data frame: data with everything to load
#' @param name Character: name of the dataset
#' @param input Shiny session input
#' @param output Shiny session output
#' 
#' @importFrom DT renderDataTable
#' @importFrom shiny downloadHandler br
#' @importFrom utils write.table
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
createDataTab <- function(index, data, name, input, output) {
    tablename <- paste("table", name, index, sep="-")
    
    table <- data[[index]]
    # Only show default columns if they are defined (don't cause problems)
    subsetToShow <- table
    
    visCols <- input[[paste(tablename, "columns", sep="-")]]
    if (!is.null(visCols)) {
        match <- visCols %in% colnames(table)
        subsetToShow <- subset(table, select=visCols[match])
    }
    
    output[[tablename]] <- renderDataTable(
        subsetToShow, style="bootstrap", selection='none', filter="top",
        options=list(pageLength=10))
    
    downloadId <- paste(tablename, "download", sep="-")
    output[[downloadId]] <- downloadHandler(
        filename = paste(name, attr(table, "tablename")),
        content = function(file) {
            res <- cbind(rownames(table), table)
            names(res)[1] <- attr(table, "dataType")
            write.table(res, file, quote=FALSE, row.names=FALSE, sep="\t")
        })
}

#' Server logic of the data module
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#'
#' @importFrom shiny selectInput tabsetPanel tags h1 h2 HTML fluidRow column
#' tagList
#'
#' @return Part of the server logic related to this tab
dataServer <- function(input, output, session) {
    ns <- session$ns
    
    analysesDescription <- tagList(
        fluidRow(
            column(3, style="padding: 5px !important;",
                   h4("Differential splicing analysis"),
                   "Analyse alternative splicing quantification based on variance",
                   "and median statistical tests. The groups available for",
                   "differential analysis comprise sample types (e.g. normal",
                   "versus tumour) and clinical attributes of patients (e.g.",
                   "tumour stage)."),
            column(3, style="padding: 5px !important;",
                   h4("Gene, transcript and protein information"),
                   "For a given splicing event, examine its gene's annotation",
                   "and corresponding transcripts and proteins. Related",
                   "research articles are also available."),
            column(3, style="padding: 5px !important;",
                   h4("Principal component analysis (PCA)"),
                   "Explore alternative splicing quantification groups using",
                   "associated clinical attributes."),
            column(3, style="padding: 5px !important;",
                   h4("Survival analysis"),
                   "Analyse survival based on clinical attributes (e.g. tumour",
                   "stage, gender and race). Additionally, study the impact of",
                   "the quantification of a single alternative splicing event",
                   "on patient survivability.")))
    
    tcga <- tags$abbr(title="The Cancer Genome Atlas", "TCGA")
    gtex <- tags$abbr(title="Genotype-Tissue Expression", "GTEx")
    
    welcome <- tagList(
        h1("Welcome"), HTML(paste0(
            "Analyse alternative splicing based on transcriptomic and ",
            "clinical data from The Cancer Genome Atlas (", tcga, ") or the ",
            "Genotype-Tissue Expression (", gtex, ") project.")),
        tags$br(), tags$br(), tags$ol(
            id="list",
            tags$li("Load clinical data and alternative splicing",
                    "junction quantification from", tcga, "or", gtex, ".", 
                    tags$br(), tags$small(
                        style="color: gray;",
                        "More data types will soon be supported.")),
            tags$li("Quantify alternative splicing events based on the",
                    "values from the percentage splicing index (PSI)",
                    "metric.",
                    # "The following event types are available:",
                    # "skipped exon (SE), mutually exclusive exon",
                    # "(MXE), alternative 3' and 5' splice site (A3SS",
                    # "and A5SS) and alternative first and last exon",
                    # "(AFE and ALE).", tags$br(),
                    tags$br(), tags$small(
                        style="color: gray;",
                        "Note: as", tcga, "does not include exon-intron",
                        "junction quantification, intron retention",
                        "events are not measurable.")),
            tags$li("Explore statistically significant events or",
                    "individual events of interest through the",
                    "following analyses:")), 
        analysesDescription, br(), br(),
        p(style="text-align:right",
          tags$a(href="http://imm.medicina.ulisboa.pt/group/compbio/",
                 target="_blank", "Nuno Morais Lab, iMM"), 
          "(", tags$a(href="mailto:nunodanielagostinho@gmail.com", 
                      "Nuno Agostinho", icon("envelope-o")), ", 2015-2016)", 
          br(), "Special thanks to my lab colleagues for their work-related",
          br(), "support and supporting chatter."))
    
    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(getData()))
            welcome
        else
            uiOutput(ns("datatabs"))
    })
    
    # Render tables when data changes
    observe({
        data <- getData()
        if (!is.null(data)) {
            for (category in names(data)) {
                categoryData <- data[[category]]
                # Create data tab for each dataset in a data category
                lapply(seq_along(categoryData), createDataTab,
                       data=categoryData, category, input, output)
            }
        }
    })
    
    # Render tabs with data tables
    output$datatabs <- renderUI({
        categoryData <- getCategoryData()
        category <- getCategory()
        
        dataTablesUI <- lapply(
            seq_along(categoryData),
            function(i) {
                data <- categoryData[[i]]
                tabDataset(ns, names(categoryData)[i], 
                           paste(category, i, sep="-"), names(data),
                           attr(data, "show"), data,
                           description=attr(data, "description"))
            })
        do.call(tabsetPanel, c(id=ns("datasetTab"), dataTablesUI))
    })
    
    # Change the active dataset
    observe( setActiveDataset(input$datasetTab) )
    
    # Update patient identifiers when clinical data is available
    observe({
        clinical <- getClinicalData()
        if ( !is.null(clinical) )
            setPatientId(rownames(clinical))
        else
            setPatientId(NULL)
    })
    
    observe({
        sampleInfo <- getSampleInfo()
        if ( !is.null(sampleInfo) )
            setSampleId( rownames(sampleInfo) )
        else
            setSampleId(NULL)
    })
    
    # Match clinical data with sample information
    observe({
        patients <- getPatientId()
        samples  <- getSampleId()
        if ( !is.null(patients) && !is.null(samples) ) {
            startProgress("Matching patients with samples...", 1)
            match <- getPatientFromSample(samples, patients)
            setClinicalMatchFrom("Inclusion levels", match)
            closeProgress("Matching process concluded")
        }
    })
    
    # Run server logic from the scripts
    getServerFunctions("data", priority="localDataServer")
}

attr(dataUI, "loader") <- "app"
attr(dataServer, "loader") <- "app"