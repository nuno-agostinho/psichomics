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
    uiOutput(ns("info"))
}

noinfo <- function(output) {
    output$info <- renderUI("No information available. Try another event.")
}

#' Parse XML from Uniprot's RESTful service
#' 
#' @param XML response from Uniprot
#' 
#' @return List containing protein length and data frame of protein features
parseUniprotXML <- function(xml) {
    doc <- xmlTreeParse(xml)
    root <- xmlRoot(doc)[[1]]
    featureNodes <- getNodeSet(root, "//feature")
    proteinLength <- as.numeric(
        xmlAttrs(getNodeSet(root, "//sequence[@length]")[[1]])[[1]])
    
    # Convert list of XMLNodes to list of characters
    l <- lapply(featureNodes, function(feat) { 
        attrs <- xmlAttrs(feat)
        
        location  <- feat[[match("location", names(feat))]]
        start <- as.numeric(xmlAttrs(location[[1]]))
        # If there's no stop position, simply sum 1 to the start position
        stop  <- tryCatch(as.numeric(xmlAttrs(location[[2]])),
                          error=function(e) start+1)
        
        # Get original and variant aminoacid
        variant <- match("variation", names(feat))
        if (!is.na(variant) && !is.null(variant)) {
            original  <- xmlToList(feat[[match("original", names(feat))]])
            variation <- xmlToList(feat[[variant]])
            variant <- sprintf("%s>%s: ", original, variation)
        } else {
            variant <- NULL
        }
        return(c(attrs, start=start, stop=stop, variant=variant))
    })
    
    # Convert list of characters to data frame of characters
    feature <- plyr::ldply(l, rbind)
    for (col in 1:ncol(feature))
        feature[[col]] <- as.character(feature[[col]])
    
    feature$start <- as.numeric(feature$start)
    feature$stop <- as.numeric(feature$stop)
    return(list(proteinLength=proteinLength, feature=feature))
}

#' Plot protein features
#' 
#' @param feature Data frame: protein features
#' @param length Integer: protein length
#' 
#' @return Highcharter object
proteinHighcharts <- function(feature, length) {
    hc <- highchart() %>%
        hc_chart(type="area", zoomType="x") %>%
        hc_xAxis(title=list(text="Position (aminoacids)"), min=0,
                 max=length, allowDecimals=FALSE) %>%
        hc_yAxis(visible=FALSE) %>%
        hc_tooltip(pointFormat="<b>{series.name} {point.id}</b>
                   <br>{point.variant}{point.description}")
    
    # The diverse types of features available
    types <- unique(feature$type)
    
    featureList <- NULL
    # Reverse elements from features so the first ones (smaller Y) are above
    for (feat in nrow(feature):1) {
        feat <- feature[feat, ]
        
        # If there's no stop position, simply sum 1 to the start position
        stop <- ifelse(!is.na(feat$stop), feat$stop, feat$start + 1)
        y <- match(feat$type, types)
        
        # Create a list with two points based on this region
        temp <- list(NULL, 
                     list(x=feat$start, y=y, description=feat$description, 
                          id=feat$id, variant=feat$variant), 
                     list(x=feat$stop,  y=y, description=feat$description, 
                          id=feat$id, variant=feat$variant))
        
        # Either make a new list or append to existing
        if (is.null(featureList[feat$type])) {
            featureList[[feat$type]] <- temp[2:3]
        } else {
            featureList[[feat$type]] <- c(featureList[[feat$type]], temp)
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
                               "Zoom to splicing event"="event")),
                plotOutput(ns("plotTranscripts"), height="200px"),
                uiOutput(ns("selectizeProtein")),
                highchartOutput(ns("plotProtein"), height="200px")
            )
        })
        
        output$selectizeProtein <- renderUI({
            proteins <- info$Transcript$Translation$id
            proteins <- proteins[!is.na(proteins)]
            tagList(
                selectizeInput(ns("selectedProtein"), label="Select protein",
                               choices=proteins),
                uiOutput(ns("proteinLink"))
            )
        })
    })
    
    observe({
        # Convert from ENSEMBL to Uniprot/SWISSPROT
        ensembl <- input$selectedProtein
        
        if (is.null(ensembl)) return(NULL)
        
        species <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37 <- assembly == "hg19"
        
        print("Looking for ENSEMBL protein in Uniprot...")
        uniprot <- queryEnsembl(paste0("xrefs/id/", ensembl),
                                list("content-type"="application/json"), 
                                grch37=grch37)
        uniprot <- uniprot[grepl("SWISSPROT", uniprot$dbname), ]
        
        ensemblLink <- tags$a("Ensembl", icon("external-link"), target="_blank",
                              class="btn btn-link",
                              href=paste0("http://", if(grch37) { "grch37." }, 
                                          "ensembl.org/", species, "/",
                                          "Search/Results?q=", ensembl))
        
        if (nrow(uniprot) == 0) {
            print("No protein match with Uniprot/SWISSPROT")
            output$proteinLink <- renderUI({
                tagList(ensemblLink,
                        tags$a("No Uniprot match was found",
                               icon("chain-broken"), 
                               class="btn btn-link")
                )
            })
            output$plotProtein <- renderHighchart(NULL)
        } else {
            protein <- uniprot$primary_id[1]
            xml <- queryUniprot(protein, "xml")
            parsed <- parseUniprotXML(xml)
            proteinLength <- parsed$proteinLength
            feature <- parsed$feature
            
            hc <- proteinHighcharts(feature, proteinLength)
            output$proteinLink <- renderUI({
                tagList(ensemblLink,
                        tags$a("Uniprot", icon("external-link"),
                               target="_blank", class="btn btn-link",
                               href=paste0("http://www.uniprot.org/uniprot/",
                                           protein))
                )
            })
            output$plotProtein <- renderHighchart(hc)
            print("Uniprot/SWISSPROT protein found")
        }
    })
}

attr(infoUI, "loader") <- "plots"
attr(infoUI, "name") <- "Information"
attr(infoServer, "loader") <- "plots"