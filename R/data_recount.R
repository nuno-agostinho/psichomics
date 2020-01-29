#' @rdname appUI
#' 
#' @importFrom shiny textInput
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsPopover
#' @importFrom shinyjs hidden
recountDataUI <- function(id, panel) {
    ns <- NS(id)
    
    title <- "SRA data loading"
    panel(style="info", title=list(icon("plus-circle"), title), value=title, 
          uiOutput(ns("recountDataModal")),
          helpText(
              "Gene expression, junction quantification and sample metadata",
              "from select SRA projects are downloaded through the",
              a(href="https://jhubiostatistics.shinyapps.io/recount/",
                target="_blank", "recount"), "R package."),
          div(class="alert", class="alert-info", role="alert", 
              "SRA data unlisted below may be manually aligned and loaded.",
              tags$a(
                  href=paste0("https://nuno-agostinho.github.io/psichomics/",
                              "articles/custom_data.html"),
                  class="alert-link", target="_blank", "Learn more...")),
          div(id=ns("loading"), class="progress",
              div(class="progress-bar progress-bar-striped active",
                  role="progressbar", style="width: 100%", "Loading")),
          uiOutput(ns("sraInterface")))
}

#' Download and load SRA projects via
#' \href{https://jhubiostatistics.shinyapps.io/recount/}{recount2}
#' 
#' @param project Character: SRA project identifiers (check
#' \code{\link[recount]{recount_abstract}})
#' @param outdir Character: directory to store the downloaded files
#' 
#' @importFrom recount download_study
#' @importFrom data.table fread
#' @importFrom SummarizedExperiment assay seqnames start end strand
#'
#' @family functions associated with SRA data retrieval
#' @return List with loaded projects
#' @export
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
    
    downloadRequiredSRAfiles <- function(folder, geneExpr, junctionQuant, 
                                         sampleInfo, sra) {
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
    }
    
    loadGeneExpression <- function(geneExpr, recountEnv) {
        load(geneExpr, recountEnv)
        geneExpr <- as.data.frame(assay(recountEnv$rse_gene))
        geneExpr <- addObjectAttrs(geneExpr,
                                   "rowNames"=1,
                                   "tablename"="Gene expression",
                                   "description"="Gene expression",
                                   "dataType"="Gene expression",
                                   "rows"="genes",
                                   "columns"="samples")
        return(geneExpr)
    }
    
    loadJunctionQuantification <- function(junctionQuant, recountEnv) {
        load(junctionQuant, recountEnv)
        rse_jx <- recountEnv$rse_jx
        junctionQuant <- as.data.frame(assay(rse_jx))
        
        ## Remove non-canonical chromosomes
        valid  <- as.vector(seqnames(rse_jx)) %in% 
            paste0("chr", c(seq(22), "X", "Y", "M"))
        
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
        return(junctionQuant)
    }
    
    project <- project[available]
    for (sra in project) {
        # Download required files if needed
        folder        <- file.path(outdir, sra)
        geneExpr      <- file.path(folder, "rse_gene.Rdata")
        junctionQuant <- file.path(folder, "rse_jx.Rdata")
        sampleInfo    <- file.path(folder, paste0(sra, ".tsv"))
        downloadRequiredSRAfiles(
            folder, geneExpr, junctionQuant, sampleInfo, sra)
        
        recountEnv  <- new.env()
        data[[sra]] <- list()
        
        # Load gene expression
        updateProgress(paste("Loading", sra), detail="Gene expression")
        data[[sra]][["Gene expression"]] <- loadGeneExpression(
            geneExpr, recountEnv)
        
        # Load junction quantification
        updateProgress(paste("Loading", sra), detail="Junction quantification")
        data[[sra]][["Junction quantification"]] <- loadJunctionQuantification(
            junctionQuant, recountEnv)
        
        # Sample metadata
        updateProgress(paste("Loading", sra), detail="Sample metadata")
        format <- loadFileFormats()$recountSampleFormat
        data[[sra]][["Sample metadata"]] <- parseValidFile(sampleInfo, format)
        
        attr(data[[sra]], "source") <- "recount"
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
            project <- input$project
        })
        
        if (is.null(project)) {
            errorModal(session, "No project selected",
                       "Select at least one project to load.",
                       caller="Load recount data", modalId="recountDataModal")
        } else if (!dir.exists(outdir)) {
            errorModal(session, "Folder not found",
                       "The selected folder", 
                       tags$kbd(prepareWordBreak(outdir)), "was not found.",
                       caller="Load recount data", modalId="recountDataModal")
        } else if (!is.null(getData())) {
            loadedDataModal(session, "recountDataModal", 
                            "recountDataReplace", "recountDataAppend")
        } else {
            startProcess("loadRecountData")
            data <- loadSRAproject(project, outdir)
            setData(data)
            endProcess("loadRecountData")
        }
    })
    
    observeEvent(input$recountDataReplace, {
        isolate({
            outdir  <- input$dataFolder
            project <- input$project
        })
        
        startProcess("loadRecountData")
        data <- loadSRAproject(project, outdir)
        setData(data)
        endProcess("loadRecountData")
    })
    
    observeEvent(input$recountDataAppend, {
        isolate({
            outdir  <- input$dataFolder
            project <- input$project
        })
        
        startProcess("loadRecountData")
        data <- loadSRAproject(project, outdir)
        data <- processDatasetNames(c(getData(), data))
        setData(data)
        endProcess("loadRecountData")
    })
}

attr(recountDataUI, "loader") <- "data"
attr(recountDataServer, "loader") <- "data"