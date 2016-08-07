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
#' @return Parsed response or NULL if there's no response
#'
#' @examples
#' path  <- "overlap/region/human/7:140424943-140624564"
#' query <- list(feature = "gene")
#' psichomics:::queryEnsembl(path, query, grch37 = TRUE)
#'
#' path  <- "lookup/symbol/human/BRCA2"
#' query <- list(expand=1)
#' psichomics:::queryEnsembl(path, query, grch37 = TRUE)
queryEnsembl <- function(path, query, grch37 = TRUE) {
    url <- paste0("http://", if(grch37) "grch37.", "rest.ensembl.org")
    resp <- GET(url, path=path, query=query)
    if (http_error(resp)) return(NULL)
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
#'
#' @examples
#' protein <- "P51587"
#' format <- "xml"
#' psichomics:::queryUniprot(protein, format)
queryUniprot <- function(protein, format="xml") {
    url <- "http://www.uniprot.org"
    path <- paste0("uniprot/", protein, ".", format)
    resp <- GET(url, path=path)
    warn_for_status(resp)
    r <- content(resp, "text", encoding = "UTF8")
    return(r)
}

#' Convert a protein's Ensembl identifier to UniProt identifier
#' 
#' @param protein Character: Ensembl protein identifier
#' 
#' @return UniProt protein identifier
#' @export
ensemblToUniprot <- function(protein) {
    external <- queryEnsembl(paste0("xrefs/id/", protein),
                             list("content-type"="application/json"),
                             grch37=TRUE)
    db <- external[grepl("Uniprot", external$dbname), ]
    uniprot <- db$primary_id
    names(uniprot) <- sprintf("%s (%s)", db$display_id, db$db_display_name)
    return(uniprot)
}

#' Information's user interface
#' @param id Character: identifier
#' @importFrom shiny uiOutput
#' @return HTML elements
infoUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("info"))
}

#' Interface when no information could be retrieved
#' @param output Shiny output
#' @param title Character: title of the message to show to the user
#' @param description Character: description of the message to show to the user
#' @importFrom shiny renderUI h3 br tags
noinfo <- function(output, title="No information available for this event.",
                   description="Select another alternative splicing event.") {
    output$info <- renderUI( h3(title, br(), tags$small(description)) )
}

#' Parse XML from Uniprot's RESTful service
#'
#' @param xml response from Uniprot
#'
#' @importFrom XML xmlTreeParse xmlRoot getNodeSet xmlAttrs xmlToList
#' @importFrom plyr ldply
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
    feature <- ldply(l, rbind)
    for (col in 1:ncol(feature))
        feature[[col]] <- as.character(feature[[col]])
    
    feature$start <- as.numeric(feature$start)
    feature$stop <- as.numeric(feature$stop)
    return(list(proteinLength=proteinLength, feature=feature))
}

#' Plot protein features
#'
#' @param protein Character: UniProt protein identifier
#'
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip
#' hc_add_series
#'
#' @return highchart object
#' @export
plotProtein <- function(protein) {
    xml     <- queryUniprot(protein, "xml")
    parsed  <- parseUniprotXML(xml)
    length  <- parsed$proteinLength
    feature <- parsed$feature
    
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

#' Plot transcripts
#' 
#' @param info Information retrieved from ENSEMBL
#' @param eventPosition Numeric: coordinates of the alternative splicing event
#' 
#' @importFrom shiny renderPlot
#' @export
plotTranscripts <- function(info, eventPosition) {
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
    
    
    plotGenes(transcripts, chrom, min, max, fontsize=1.5)
    zoomsregion(eventPosition, highlight=TRUE)
    labelgenome(chrom, min, max, scale="Mb")
}

#' Render genetic information
#' 
#' @param ns Namespace function
#' @param info Information as retrieved from ENSEMBL
#' @param species Character: species name
#' @param assembly Character: assembly version
#' @param grch37 Boolean: use version GRCh37 of the genome?
#' 
#' @importFrom shiny renderUI h2 h4 plotOutput
renderGeneticInfo <- function(ns, info, species, assembly, grch37) {
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
               href=paste0("http://", if(grch37) { "grch37." },
                           "ensembl.org/", species, "/",
                           "Gene/Summary?g=", info$id)),
        tags$a("UCSC", icon("external-link"), target="_blank",
               class="btn btn-link",
               href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks",
                           if(grch37) { "?db=hg19" }, "&position=chr",
                           info$seq_region_name, ":", start, "-", end)),
        tags$a("GeneCards", icon("external-link"), target="_blank",
               class="btn btn-link",
               href=paste0("http://www.genecards.org/cgi-bin/",
                           "carddisp.pl?gene=", info$id)),
        if (species == "human") 
            tags$a("Human Protein Atlas", icon("external-link"),
                   target="_blank", class="btn btn-link",
                   href=paste0("http://www.proteinatlas.org/", info$id, 
                               "/cancer")),
        h4("Transcripts"),
        plotOutput(ns("plotTranscripts"), height="200px"),
        uiOutput(ns("selectizeProtein")),
        highchartOutput(ns("plotProtein"), height="200px")
    )
}

#' Query information from Ensembl by a given alternative splicing event
#' 
#' @param event Character: alternative splicing event identifier
#' @param ... Arguments to pass to \code{queryEnsemblByGene}
#' 
#' @return Information from Ensembl
#' @export
queryEnsemblByEvent <- function(event, ...) {
    gene <- parseEvent(event)$gene
    if (gene == "NA")
        stop("This event has no gene associated")
    return(queryEnsemblByGene(gene, ...))
}

#' Query information from Ensembl by a given gene
#' 
#' @param gene Character: gene identifier
#' @param species Character: species
#' @param assembly Character: assembly version
#' 
#' @return Information from Ensembl
#' @export
queryEnsemblByGene <- function(gene, species, assembly) {
    grch37 <- assembly == "hg19"
    path   <- paste0("lookup/symbol/", species, "/", gene)
    info   <- queryEnsembl(path, list(expand=1), grch37=grch37)
    return(info)
}

#' Server logic
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom Sushi plotGenes zoomsregion labelgenome
#' @importFrom highcharter highchart %>%
infoServer <- function(input, output, session) {
    ns <- session$ns
    
    observe({
        event <- getEvent()
        if (is.null(getInclusionLevels()))
            return(noinfo(output, "Quantify alternative splicing",
                          paste("To perform this analysis, you need to quantify",
                                "alternative splicing.")))
        else if (is.null(event) || event == "") return(noinfo(output))
        
        parsed <- parseEvent(event)
        if (parsed$gene == "NA") return(noinfo(output))
        
        species  <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37   <- assembly == "hg19"
        if(is.null(species) || is.null(assembly)) stop("NULL")
        
        path <- paste0("lookup/symbol/", species, "/", parsed$gene)
        info <- queryEnsembl(path, list(expand=1), grch37=grch37)
        
        if (is.null(info))
            return(noinfo(output, title="No Ensembl match retrieved.",
                          description="Please, try again."))
        
        output$info <- renderUI(
            renderGeneticInfo(ns, info, species, assembly, grch37))
        
        output$plotTranscripts <- renderPlot(
            plotTranscripts(info, parsed$pos))
        
        # Show NULL so it doesn't show previous results when loading
        output$selectizeProtein <- renderUI("Loading...")
        output$plotProtein <- renderHighchart(NULL)
        
        output$selectizeProtein <- renderUI({
            proteins <- info$Transcript$Translation$id
            names(proteins) <- paste("Transcript:", info$Transcript$id, "/", 
                                     "Protein:", proteins)
            proteins <- proteins[!is.na(proteins)]
            tagList(
                selectizeInput(ns("selectedProtein"), label="Select protein",
                               choices=proteins, width="auto"),
                uiOutput(ns("proteinLink")))
        })
    })
    
    observe({
        ensembl <- input$selectedProtein
        if (is.null(ensembl)) return(NULL)
        
        species  <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37   <- assembly == "hg19"
        if(is.null(species) || is.null(assembly)) return(NULL)
        
        # Convert from ENSEMBL to Uniprot/SWISSPROT
        print("Looking for ENSEMBL protein in Uniprot...")
        uniprot <- queryEnsembl(paste0("xrefs/id/", ensembl),
                                list("content-type"="application/json"),
                                grch37=grch37)
        if (is.null(uniprot)) return(NULL)
        uniprot <- uniprot[grepl("SWISSPROT", uniprot$dbname), ]
        
        href <- paste0("http://", if(grch37) { "grch37." }, "ensembl.org/", 
                       species, "/Search/Results?q=", ensembl)
        ensemblLink <- tags$a("Ensembl", icon("external-link"), target="_blank",
                              class="btn btn-link", href=href)
        
        if (nrow(uniprot) == 0) {
            print("No protein match with Uniprot/SWISSPROT")
            noMatch <- tagList(ensemblLink, tags$a("No Uniprot match", 
                                                   icon("chain-broken"),
                                                   class="btn btn-link"))
            output$proteinLink <- renderUI(noMatch)
            output$plotProtein <- renderHighchart(NULL)
        } else {
            protein <- uniprot$primary_id[1]
            hc <- plotProtein(protein)
            output$proteinLink <- renderUI({
                href <- paste0("http://www.uniprot.org/uniprot/", protein)
                tagList(ensemblLink, 
                        tags$a("Uniprot", icon("external-link"), href=href,
                               target="_blank", class="btn btn-link"))
            })
            output$plotProtein <- renderHighchart(hc)
            print("Uniprot/SWISSPROT protein found")
        }
    })
}

attr(infoUI, "loader") <- "analysis"
attr(infoUI, "name") <- "Gene, transcript and protein information"
attr(infoUI, "selectEvent") <- TRUE
attr(infoServer, "loader") <- "analysis"