#' @rdname appUI
#' 
#' @importFrom shiny textInput
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsPopover
#' @importFrom shinyjs hidden
recountDataUI <- function(id, panel) {
    ns <- NS(id)
    
    title <- "Automatically load SRA data"
    panel(style="info", title=list(icon("plus-circle"), title), value=title, 
          uiOutput(ns("recountDataModal")),
          helpText(
              "Gene expression, junction quantification and sample metadata",
              "from select SRA projects are downloaded through the",
                   a(href="https://jhubiostatistics.shinyapps.io/recount/",
                     target="_blank", "recount"), "R package."),
          div(class="alert", class="alert-info", role="alert", 
              "Data from SRA projects not listed below may be manually loaded",
              "after splice-aware alignment.",
              tags$a(
                  href="http://rpubs.com/nuno-agostinho/psichomics-custom-data",
                  class="alert-link", target="_blank", "Learn more...")),
          div(id=ns("loading"), class="progress",
              div(class="progress-bar progress-bar-striped active",
                  role="progressbar", style="width: 100%", "Loading")),
          uiOutput(ns("sraInterface")))
}

#' Download and load SRA projects
#' 
#' @param project Character: SRA project identifiers to download
#' @param outdir Character: directory to store the downloaded files
#' 
#' @importFrom recount download_study
#' @importFrom data.table fread
#' @importFrom SummarizedExperiment assay seqnames start end strand
#' 
#' @return List containing downloaded projects
loadSRAproject <- function(project, outdir=getDownloadsFolder()) {
    data <- list()
    
    # Warn about SRA projects whose data are not available
    project   <- unique(project)
    available <- project %in% recount::recount_abstract$project
    
    if (all(!available)) {
        stop("No data found for the given SRA projects.")
    } else if (any(!available)) {
        warning("No data found for the following SRA projects: ", 
                paste(project[!available], collapse=", "))
    }
    
    project <- project[available]
    for (sra in project) {
        # Check which files that were not yet downloaded
        folder        <- file.path(outdir, sra)
        geneExpr      <- file.path(folder, "rse_gene.Rdata")
        junctionQuant <- file.path(folder, "rse_jx.Rdata")
        sampleInfo    <- file.path(folder, paste0(sra, ".tsv"))
        
        fileType <- c("rse-gene"=geneExpr,
                      "rse-jx"=junctionQuant, 
                      "phenotype"=sampleInfo)
        fileTypeToDownload <- !sapply(fileType, file.exists)
        fileTypeToDownload <- names(fileType)[fileTypeToDownload]
        
        len <- length(fileTypeToDownload)
        updateProgress(paste("Loading", sra), divisions=3 + len)
        if (len > 0) {
            downloadProject <- function(type, sra, folder, ...) {
                detail <- switch(type, "rse-gene"="Gene expression",
                                 "rse-jx"="Junction quantification",
                                 "phenotype"="Sample metadata")
                updateProgress(paste("Downloading", sra), detail=detail)
                download_study(sra, outdir=folder, type=type, ...)
            }
            
            sapply(fileTypeToDownload, downloadProject, sra=sra, folder=folder)
        }
        
        recountEnv  <- new.env()
        data[[sra]] <- list()
        
        # Gene expression
        updateProgress(paste("Loading", sra), detail="Gene expression")
        load(geneExpr, recountEnv)
        geneExpr <- as.data.frame(assay(recountEnv$rse_gene))
        geneExpr <- addObjectAttrs(
            geneExpr,
            "rowNames"=1,
            "tablename"="Gene expression",
            "description"="Gene expression",
            "dataType"="Gene expression",
            "rows"="genes",
            "columns"="samples")
        data[[sra]][["Gene expression"]] <- geneExpr
        
        # Junction quantification
        updateProgress(paste("Loading", sra), detail="Junction quantification")
        load(junctionQuant, recountEnv)
        rse_jx <- recountEnv$rse_jx
        junctionQuant <- as.data.frame(assay(rse_jx))
        
        ## Remove non-canonical chromosomes
        valid  <- as.vector(seqnames(rse_jx)) %in% 
            paste0("chr", c(1:22, "X", "Y", "M"))
        
        chr    <- as.vector(seqnames(rse_jx))[valid]
        start  <- start(rse_jx)[valid] - 1
        end    <- end(rse_jx)[valid]   + 1
        strand <- as.vector(strand(rse_jx))[valid]
        
        junctionQuant           <- junctionQuant[valid, , drop=FALSE]
        rownames(junctionQuant) <- paste(chr, start, end, strand, sep=":")
        junctionQuant <- addObjectAttrs(
            junctionQuant,
            "rowNames"=1,
            "tablename"="Junction quantification",
            "description"="Read counts of splicing junctions",
            "dataType"="Junction quantification",
            "rows"="splice junctions",
            "columns"="samples")
        data[[sra]][["Junction quantification"]] <- junctionQuant
        
        # Sample metadata
        updateProgress(paste("Loading", sra), detail="Sample metadata")
        format <- loadFileFormats()$recountSampleFormat
        data[[sra]][["Sample metadata"]] <- parseValidFile(sampleInfo, format)
        
        closeProgress()
    }
    return(data)
}

#' @rdname appServer
#' @importFrom shiny updateTextInput
#' @importFrom shinyjs hide show
recountDataServer <- function(input, output, session) {
    ns <- session$ns
    
    output$sraInterface <- renderUI({
        data           <- recount::recount_abstract
        choices        <- data$project
        names(choices) <- sprintf("%s (%s samples)", data$project,
                                  data$number_samples)
        
        ui <- tagList(
            selectizeInput(
                ns("project"), "SRA project", NULL, choices=choices,
                width="100%", multiple=TRUE, options=list(
                    placeholder="Select SRA project(s)",
                    plugins=list("remove_button"))),
            fileBrowserInput(
                ns("dataFolder"), "Folder where data is stored",
                value=getDownloadsFolder(), placeholder="No folder selected",
                info=TRUE, infoFUN=bsTooltip,
                infoTitle=paste("Data will be downloaded if not available in",
                                "this folder.")),
            tags$a(href="https://jhubiostatistics.shinyapps.io/recount/",
                   class="btn btn-default", role="button", target="_blank",
                   tags$i(class="fa fa-external-link"),
                   "Check available datasets"),
            processButton(ns("loadRecountData"), "Load data"))
        hide("loading")
        return(ui)
    })
    
    observeEvent(input$loadRecountData, {
        isolate({
            outdir  <- input$dataFolder
            project <- isolate(input$project)
        })
        
        startProcess("loadRecountData")
        
        if (!is.null(getData())) {
            loadedDataModal(session, "recountDataModal", 
                            "recountDataReplace", "recountDataAppend")
        } else {
            data <- loadSRAproject(project, outdir)
            setData(data)
        }
        endProcess("loadRecountData")
    })
    
    observeEvent(input$recountDataReplace, {
        isolate({
            outdir  <- input$dataFolder
            project <- isolate(input$project)
        })
        
        data <- loadSRAproject(project, outdir)
        setData(data)
    })
    
    observeEvent(input$recountDataAppend, {
        isolate({
            outdir  <- input$dataFolder
            project <- isolate(input$project)
        })
        
        data <- loadSRAproject(project, outdir)
        data <- processDatasetNames(c(getData(), data))
        setData(data)
    })
}

attr(recountDataUI, "loader") <- "data"
attr(recountDataServer, "loader") <- "data"