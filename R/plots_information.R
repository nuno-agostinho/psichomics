## TODO: check if API is alive before querying data

#' Query the Ensembl REST API
#' 
#' @param path Character: API path
#' @param query Character: API query
#' @param grch37 Boolean: query the Ensembl GRCh37 API? TRUE by default;
#' otherwise, query the most recent API
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' 
#' @return Parsed response
#' @export
#' 
#' @examples 
#' path  <- "overlap/region/human/7:140424943-140624564"
#' query <- list(feature = "gene")
#' queryEnsembl(path, query, grch37 = TRUE)
#' 
#' path  <- "lookup/symbol/human/BRCA2"
#' query <- list(expand=1)
#' queryEnsembl(path, query, grch37 = TRUE)
queryEnsembl <- function(path, query, grch37 = TRUE) {
    url <- paste0("http://", if(grch37) "grch37.", "rest.ensembl.org")
    resp <- GET(url, path=path, query=query)
    warn_for_status(resp)
    r <- content(resp, "text", encoding = "UTF8")
    return(fromJSON(r))
}

#' Query the Uniprot REST API
#' 
#' @param protein Character: protein to query
#' @param format Character: format of the response
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' 
#' @return Parsed response
#' @export
#' 
#' @examples 
#' protein <- "P51587"
#' format <- "xml"
#' queryUniprot(protein, format)
queryUniprot <- function(protein, format="xml") {
    url <- "http://www.uniprot.org"
    path <- paste0("uniprot/", protein, ".", format)
    resp <- GET(url, path=path)
    warn_for_status(resp)
    r <- content(resp, "text", encoding = "UTF8")
    return(r)
}

infoUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("info")),
        plotOutput(ns("plotTranscripts")),
        uiOutput(ns("selectProtein")),
        highchartOutput(ns("plotProtein"))
    )
}

noinfo <- function(output) {
    output$info <- renderUI("No information available. Try another event.")
    output$plotTranscripts <- renderPlot(NULL)
}

#' Plot a Uniprot protein
#' 
#' @param feature List of XML nodes: features of the Uniprot protein
#' @param proteinLength Integer: protein length
#' 
#' @return Highcharter object
proteinHighcharts <- function(feature, proteinLength) {
    hc <- highchart() %>%
        hc_chart(type="area", zoomType="x") %>%
        hc_xAxis(title=list(text="Position (aminoacids)"), min=0,
                 max=proteinLength) %>%
        hc_yAxis(visible=FALSE) %>%
        hc_tooltip(pointFormat="<b>{series.name} {point.id}</b>
                   <br>{point.variant}{point.description}")
    
    # The diverse types of features available
    types <- unique(sapply(sapply(feature, xmlAttrs), "[[", 1))
    
    featureList <- NULL
    # Reverse elements from features so the first ones (smaller Y) are above
    for (feat in rev(feature)) {
        attrs <- xmlAttrs(feat)
        type  <- attrs[[1]]
        description <- attrs[[2]]
        id <- tryCatch(attrs[[3]], error=function(e) NULL)
        location  <- feat[[match("location", names(feat))]]
        
        # Get original and variant aminoacid
        variant <- match("variation", names(feat))
        if (!is.na(variant) && !is.null(variant)) {
            original  <- xmlToList(feat[[match("original", names(feat))]])
            variation <- xmlToList(feat[[variant]])
            variant <- sprintf("%s>%s: ", original, variation)
        } else {
            variant <- NULL
        }
        
        start <- as.numeric(xmlAttrs(location[[1]]))
        # If there's no stop position, simply sum 1 to the start position
        stop  <- tryCatch(as.numeric(xmlAttrs(location[[2]])),
                          error=function(e) {return(start+1)})
        y     <- match(type, types)
        
        # Create a list with two points based on this region
        temp <- list(
            NULL, 
            list(x=start, y=y, description=description, id=id, variant=variant), 
            list(x=stop,  y=y, description=description, id=id, variant=variant))
        
        # Either make a new list or append to existing
        if (is.null(featureList[type])) {
            featureList[[type]] <- temp[2:3]
        } else {
            featureList[[type]] <- c(featureList[[type]], temp)
        }
    }
    for (type in names(featureList))
        hc <- hc %>% hc_add_series(name=type, data=featureList[[type]])
    return(hc)
}

#' @importFrom Sushi plotGenes zoomsregion labelgenome
#' @importFrom highcharter highchart %>%
infoServer <- function(input, output, session) {
    ns <- session$ns
    
    observe({
        event <- getEvent()
        if (is.null(event) || event == "") return(noinfo(output))
        
        event <- strsplit(event, "_")[[1]]
        gene <- event[length(event)]
        if (gene == "NA") return(noinfo(output))
        
        eventPosition <- event[4:(length(event)-1)]
        eventPosition <- range(as.numeric(eventPosition))
        
        species <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37 <- assembly == "hg19"
        
        path <- paste0("lookup/symbol/", species, "/", gene)
        info <- queryEnsembl(path, list(expand=1), grch37=grch37)
        
        output$plotTranscripts <- renderPlot({
            transcripts <- data.frame()
            
            for (i in 1:nrow(info$Transcript)) {
                transcriptId <- info$Transcript[i, "id"]
                nn <- info$Transcript$Exon[[i]]
                transcripts <- rbind(
                    transcripts,
                    data.frame(chrom=nn$seq_region_name, start=nn$start, 
                               stop=nn$end, gene=transcriptId, score=0, 
                               strand=nn$strand))
            }
            
            chrom <- paste0("chr", transcripts$chrom)
            min <- min(transcripts$start)
            max <- max(transcripts$stop)
            
            if (!is.null(input$zoom) && input$zoom == "all") {
                plotGenes(transcripts, chrom, min, max, fontsize=1.5)
                zoomsregion(eventPosition, highlight=TRUE)
                labelgenome(chrom, min, max, scale="Mb")
            } else if (!is.null(input$zoom) && input$zoom == "event") {
                min <- eventPosition[1]
                max <- eventPosition[2]
                plotGenes(transcripts, chrom, min, max, fontsize=1.5)
                labelgenome(chrom, min, max, scale="Mb")
            }
        })
        
        output$info <- renderUI({
            start <- as.numeric(info$start)
            end   <- as.numeric(info$end)
            
            tagList(
                h2(info$display_name, tags$small(info$id)),
                sprintf("Species: %s (assembly %s)", species, assembly),
                br(), sprintf("Chromosome %s: %s-%s (%s strand)",
                              info$seq_region_name, 
                              format(start, big.mark=",", scientific=FALSE),
                              format(end, big.mark=",", scientific=FALSE),
                              ifelse(info$strand == -1,"reverse", "forward")),
                br(), sprintf("%s (%s)", info$description,
                              info$biotype), br(),
                tags$a("Ensembl", icon("external-link"), target="_blank",
                       class="btn btn-link",
                       href=paste0("http://grch37.ensembl.org/", species, "/",
                                   "Gene/Summary?g=", info$id)),
                tags$a("UCSC", icon("external-link"), target="_blank",
                       class="btn btn-link",
                       href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks",
                                   "?db=hg19&position=chr", 
                                   info$seq_region_name, ":", start, "-", end)),
                radioButtons(ns("zoom"), "Zoom", inline=TRUE,
                             c("Show all transcripts"="all", 
                               "Zoom to splicing event"="event"))
            )
        })
        
        output$selectProtein <- renderUI({
            proteins <- info$Transcript$Translation$id
            proteins <- proteins[!is.na(proteins)]
            selectizeInput(ns("protein"), label="Select protein",
                           choices=proteins)
        })
        
        observe({
            # Convert from ENSEMBL to Uniprot/SWISSPROT
            ensembl <- input$protein
            
            if (is.null(ensembl)) return(NULL)
            print("Looking for ENSEMBL protein in Uniprot...")
            uniprot <- queryEnsembl(paste0("xrefs/id/", ensembl),
                                    list("content-type"="application/json"), 
                                    grch37=grch37)
            uniprot <- uniprot[grepl("SWISSPROT", uniprot$dbname), ]
            
            if (nrow(uniprot) == 0)
                print("No protein from Uniprot :(")
            else {
                protein <- uniprot$primary_id[1]
                resp <- queryUniprot(protein, "xml")
                
                doc <- xmlTreeParse(resp)
                root <- xmlRoot(doc)[[1]]
                feature <- getNodeSet(root, "//feature")
                proteinLength <- as.numeric(
                    xmlAttrs(getNodeSet(root, "//sequence")[[1]])[[1]])
                
                hc <- proteinHighcharts(feature, proteinLength)
                output$plotProtein <- renderHighchart(hc)
                print("Uniprot protein found")
            }
        })
    })
}

attr(infoUI, "loader") <- "plots"
attr(infoUI, "name") <- "Information"
attr(infoServer, "loader") <- "plots"