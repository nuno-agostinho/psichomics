## TODO: check if API is alive before querying data

#' Query the Ensembl REST API
#'
#' @param path Character: API path
#' @param query Character: API query
#' @param grch37 Boolean: query the Ensembl GRCh37 API? TRUE by default;
#' otherwise, query the most recent API
#'
#' @importFrom httr GET timeout
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
    resp <- tryCatch(GET(url, path=path, query=query, config=timeout(10)), 
                     error=return)
    
    if (is(resp, "error") || http_error(resp) || is.null(resp)) # Time out
        return(NULL) # for instance, time out
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

#' Query the PubMed REST API
#' 
#' @param primary Character: primary search term
#' @param ... Character: other relevant search terms
#' @param top Numeric: number of articles to retrieve (3 by default)
#' @param field Character: field of interest where to look for terms ("abstract"
#'  by default)
#' @param sort Character: sort by a given parameter ("relevance" by default)
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#'
#' @return Parsed response
#'
#' @examples
#' psichomics:::queryPubMed("BRCA1", "cancer", "adrenocortical carcinoma")
queryPubMed <- function(primary, ..., top=3, field="abstract", 
                        sort="relevance") {
    args  <- unlist(list(...))
    for (each in seq_along(args)) {
        args[each] <- paste0("(", paste(primary, args[each], sep=" AND "), ")")
    }
    terms <- paste(args, collapse=" OR ")
    terms <- paste(primary, terms, sep=" OR ")
    
    url <- "https://eutils.ncbi.nlm.nih.gov"
    query <- list(db="pmc", term=terms, retmax=top, tool="psichomics", 
                  field=field, sort=sort,
                  email="nunodanielagostinho@gmail.com", retmode="json")
    resp <- GET(url, path="entrez/eutils/esearch.fcgi", query=query)
    warn_for_status(resp)
    search <- content(resp, "text", encoding = "UTF8")
    search <- fromJSON(search)[[2]]
    
    # Get summary information on the articles
    ids <- paste(search$idlist, collapse="+")
    query <- list(db="pmc", tool="psichomics", id=ids,
                  email="nunodanielagostinho@gmail.com", retmode="json")
    resp <- GET(url, path="entrez/eutils/esummary.fcgi", query=query)
    warn_for_status(resp)
    metadata <- content(resp, "text", encoding = "UTF8")
    metadata <- fromJSON(metadata)[[2]][-1]
    return(c(search=list(search), metadata))
}

#' Convert a protein's Ensembl identifier to UniProt identifier
#' 
#' @param protein Character: Ensembl protein identifier
#' 
#' @return UniProt protein identifier
#' @export
#' @examples 
#' ensemblToUniprot("ENSP00000445929")
ensemblToUniprot <- function(protein) {
    external <- queryEnsembl(paste0("xrefs/id/", protein),
                             list("content-type"="application/json"),
                             grch37=TRUE)
    
    if (is.null(external)) return(NULL)
        
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
    tagList(uiOutput(ns("geneSelection")),
            uiOutput(ns("info")))
    
}

#' Interface when no information could be retrieved
#' @param output Shiny output
#' @param title Character: title of the message to show to the user
#' @param description Character: description of the message to show to the user
#' @importFrom shiny renderUI h3 br tags
#' @return NULL (this function is used to modify the Shiny session's state)
noinfo <- function(output, title=paste("No information available for the gene",
                                       "associated with this event."),
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
    if (ncol(feature) > 0) {
        for (col in 1:ncol(feature))
            feature[[col]] <- as.character(feature[[col]])
        
        feature$start <- as.numeric(feature$start)
        feature$stop <- as.numeric(feature$stop)
    }
    
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
#' @examples
#' \dontrun{
#' plotProtein("P38398")
#' }
plotProtein <- function(protein) {
    cat("Retrieving protein annotation from UniProt...", fill=TRUE)
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
                   <br>{point.variant}{point.description}") %>%
        export_highcharts()
    
    # The diverse types of features available
    types <- unique(feature$type)
    
    cat("Plotting protein domains...", fill=TRUE)
    featureList <- NULL
    
    if (nrow(feature) == 0)
        stop("This protein has no annotated domains in UniProt.")
    
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
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
#' @examples
#' event <- "SE_12_-_7985318_7984360_7984200_7982602_SLC2A14"
#' info  <- queryEnsemblByEvent(event, species="human", assembly="hg19")
#' pos   <- parseSplicingEvent(event)$pos[[1]]
#' \dontrun{
#' plotTranscripts(info, pos)
#' }
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
#' @param species Character: species name (NULL by default)
#' @param assembly Character: assembly version (NULL by default)
#' @param grch37 Boolean: use version GRCh37 of the genome? FALSE by default
#' 
#' @importFrom shiny renderUI h2 h4 plotOutput
#' @return HTML elements to render gene, protein and transcript annotation
renderGeneticInfo <- function(ns, info, species=NULL, assembly=NULL, 
                              grch37=FALSE) {
    start <- as.numeric(info$start)
    end   <- as.numeric(info$end)
    
    if (!is.null(species))
        ensembl <- tags$a("Ensembl", icon("external-link"), target="_blank",
                          href=paste0("http://", if(grch37) { "grch37." }, 
                                      "ensembl.org/", species, "/", 
                                      "Gene/Summary?g=", info$id))
    else
        ensembl <- NULL
    
    links <- tagList(
        ensembl,
        tags$a("UCSC", icon("external-link"), target="_blank",
               href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks",
                           if(grch37) { "?db=hg19" }, "&position=chr",
                           info$seq_region_name, ":", start, "-", end)),
        tags$a("GeneCards", icon("external-link"), target="_blank",
               href=paste0("http://www.genecards.org/cgi-bin/",
                           "carddisp.pl?gene=", info$id)),
        if (species == "human") { 
            tags$a("Human Protein Atlas (Cancer Atlas)", 
                   icon("external-link"), target="_blank",
                   href=paste0("http://www.proteinatlas.org/", info$id, 
                               "/cancer"))
        },
        tags$a("VAST-DB", icon("external-link"), target="_blank",
               href=paste0("http://vastdb.crg.eu/wiki/Gene:", info$id, 
                           "@Genome:", assembly))
    )
    
    dtWidth  <- "width: 80px;"
    ddMargin <- "margin-left: 100px;"
    
    if (!is.null(species)) {
        if (!is.null(assembly))
            speciesInfo <- sprintf("%s (%s assembly)", species, assembly)
        else
            speciesInfo <- species
    } else {
        speciesInfo <- "No species defined"
    }
    
    genetic <- tagList(
        h2(info$display_name, tags$small(info$id)),
        tags$dl(class="dl-horizontal",
                tags$dt(style=dtWidth, "Species"),
                tags$dd(style=ddMargin, speciesInfo),
                tags$dt(style=dtWidth, "Location"),
                tags$dd(style=ddMargin, 
                        sprintf("Chromosome %s: %s-%s (%s strand)",
                                info$seq_region_name,
                                format(start, big.mark=",", scientific=FALSE),
                                format(end, big.mark=",", scientific=FALSE),
                                ifelse(info$strand == -1,
                                       "reverse", "forward"))),
                tags$dt(style=dtWidth, "Description"),
                tags$dd(style=ddMargin,
                        sprintf("%s (%s)", info$description, info$biotype)),
                tags$dt(style=dtWidth, "Links"),
                tags$dd(style=ddMargin,
                        tags$ul(class="list-inline", lapply(links, tags$li)) ))
    )
    
    tagList(
        fluidRow(column(6, genetic),
                 column(6, uiOutput(ns("articles")))),
        h4("Transcripts"), 
        plotOutput(ns("plotTranscripts"), height="200px"),
        uiOutput(ns("selectizeProtein")),
        highchartOutput(ns("plotProtein"), height="200px"))
}

#' Query information from Ensembl by a given alternative splicing event
#' 
#' @param event Character: alternative splicing event identifier
#' @param ... Arguments to pass to \code{queryEnsemblByGene}
#' 
#' @return Information from Ensembl
#' @export
#' @examples 
#' event <- c("SE_17_-_41251792_41249306_41249261_41246877_BRCA1")
#' queryEnsemblByEvent(event, species="human", assembly="hg19")
queryEnsemblByEvent <- function(event, ...) {
    gene <- parseEvent(event)$gene
    if (gene == "Hypothetical")
        stop("This event has no associated gene")
    return(queryEnsemblByGene(gene, ...))
}

#' Query information from Ensembl by a given gene
#' 
#' @param gene Character: gene identifier
#' @param species Character: species (can be NULL when handling an ENSEMBL
#' identifier)
#' @param assembly Character: assembly version (can be NULL when handling an
#' ENSEMBL identifier)
#' 
#' @return Information from Ensembl
#' @export
#' @examples 
#' queryEnsemblByGene("BRCA1", "human", "hg19")
#' queryEnsemblByGene("ENSG00000139618")
queryEnsemblByGene <- function(gene, species=NULL, assembly=NULL) {
    if ( grepl("^ENSG", gene) ) {
        path   <- paste0("lookup/id/", gene)
        info <- queryEnsembl(path, list(expand=1))
    } else {
        if (is.null(species) || is.null(assembly))
            stop("Species and assembly need to be non-NULL")
        grch37 <- assembly == "hg19"
        path <- paste0("lookup/symbol/", species, "/", gene)
        info <- queryEnsembl(path, list(expand=1), grch37=grch37)
    }
    return(info)
}

#' Return the interface to display an article
#' 
#' @param article PubMed article
#' 
#' @importFrom shiny tags h5
#' 
#' @return HTML to render an article's interface
articleUI <- function(article) {
    authors <- article$authors$name
    if (length(authors) > 2) {
        authors <- paste(authors[1], "et al.")
    } else if (length(authors) == 2) {
        authors <- paste(authors, collapse=" and ")
    }
    year <- strsplit(article$pubdate, " ")[[1]][[1]]
    description <- sprintf("%s (%s). %s, %s(%s).", authors, year, 
                           article$source, article$volume, article$issue)
    
    pmid <- article$articleids$value[1]
    tags$a(href=paste0("http://pubmed.gov/", pmid), target="_blank",
           class="list-group-item", h5(class="list-group-item-heading", 
                                       article$title, tags$small(description)))
}

#' Return the interface of relevant PubMed articles for a given gene
#' 
#' @param gene Character: gene
#' @param ... Arguments to pass to \code{queryPubMed} function
#' 
#' @return HTML interface of relevant PubMed articles
pubmedUI <- function(gene, ...) {
    pubmed <- queryPubMed(gene, ...)
    
    articles <- pubmed[-1]
    articleList <- lapply(articles, articleUI)
    
    search <- pubmed$search$querytranslation
    search <- gsub("[Abstract]", "[Title/Abstract]", search, fixed = TRUE)
    search <- paste0("http://www.ncbi.nlm.nih.gov/pubmed/?term=", search)
    
    articlesUI <- div(class="panel panel-default", 
                      div(class="panel-heading",
                          tags$b("Relevant PubMed articles", tags$a(
                              href=search, target="_blank", 
                              class="pull-right", "Show more articles",
                              icon("external-link")))),
                      div(class="list-group", articleList))
    return(articlesUI)
}

#' Server logic
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom Sushi plotGenes zoomsregion labelgenome
#' @importFrom highcharter highchart %>%
#' @importFrom shiny fixedRow safeError
#' @importFrom methods is
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
infoServer <- function(input, output, session) {
    ns <- session$ns
    
    observe({
        event <- getEvent()
        if (is.null(getInclusionLevels()))
            return(noinfo(output, "Quantify alternative splicing",
                          paste("To perform this analysis, alternative splicing",
                                "must be quantified first.")))
        else if (is.null(event) || event == "") return(noinfo(output))
        
        # Select gene in case there is more than one available
        gene <- parseEvent(event)$gene[[1]]
        output$geneSelection <- renderUI({
            if (length(gene) > 1) {
                fixedRow(
                    column(3, h5(style=" margin-top: 0;",
                                 "Select one of the genes that may be",
                                 "associated with the event:")),
                    column(3, selectizeInput(ns("selectedGene"), NULL, 
                                             choices=gene)))
            }
        })
    })
    
    observe({
        event <- getEvent()
        if (is.null(getInclusionLevels()) || is.null(event) || event == "")
            return(NULL)
        
        parsed <- parseEvent(event)
        species  <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37   <- assembly == "hg19"
        
        gene <- parsed$gene[[1]]
        if (length(gene) > 1)
            gene <- input$selectedGene
        
        info <- tryCatch(queryEnsemblByGene(gene, species=species, 
                                            assembly=assembly), error=return)
        
        # Handle errors
        if (is(info, "error")) {
            if (grepl("Species and assembly", info$message)) {
                warning("No species or assembly information.")
                return(NULL)
            } else if (grepl("no associated gene", info$message)) {
                return(noinfo(output))
            }
        }
        
        if (is.null(info)) {
            output$info <- renderUI({ 
                title <- "Ensembl API appears to give no response."
                description <- paste("Please try selecting another event or",
                                     "try again later.")
                noinfo <- h3(title, br(), tags$small(description))
                fluidRow(column(6, noinfo), 
                         column(6, uiOutput(ns("articles"))))
            })
            return(NULL)
        }
        
        output$info <- renderUI(
            renderGeneticInfo(ns, info, species, assembly, grch37))
        
        output$plotTranscripts <- renderPlot(
            plotTranscripts(info, parsed$pos[[1]]))
        
        # Show NULL so it doesn't show previous results when loading
        output$selectizeProtein <- renderUI("Loading...")
        output$plotProtein <- renderHighchart(NULL)
        
        output$selectizeProtein <- renderUI({
            transcripts <- info$Transcript$Translation$id
            if (is.null(transcripts))
                return(h4("No proteins retrieved from UniProt"))
            names(transcripts) <- info$Transcript$id
            transcripts <- transcripts[!is.na(transcripts)]
            
            tagList(
                fixedRow(
                    column(3, selectizeInput(ns("selectedTranscript"),
                                             label="Select transcript", 
                                             choices=transcripts, width="auto")),
                    column(3, selectizeInput(ns("selectedProtein"),
                                             label="Select protein",
                                             choices=NULL, width="auto")),
                    column(3, uiOutput(ns("proteinLink"), 
                                       class="inline_selectize"))))
        })
    })
    
    # Check UniProt proteins matching a given Ensembl transcript
    observe({
        ensembl  <- input$selectedTranscript
        species  <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37   <- assembly == "hg19"
        
        if (is.null(ensembl) || is.null(species) || is.null(assembly)) {
            updateSelectizeInput(session, "selectedProtein", 
                                 choices=c("No UniProt match"=""))
        } else {
            # Look up matching identifiers from Uniprot
            external <- queryEnsembl(paste0("xrefs/id/", ensembl),
                                     list("content-type"="application/json"),
                                     grch37=grch37)
            if (is.null(external)) {
                updateSelectizeInput(
                    session, "selectedProtein", 
                    choices=c("Ensembl API appears to be offline"=""))
            } else {
                db <- external[grepl("Uniprot", external$dbname), ]
                if (nrow(db) == 0) {
                    updateSelectizeInput(session, "selectedProtein", 
                                         choices=c("No UniProt match"=""))
                }
                uniprot <- db$primary_id
                names(uniprot) <- sprintf("%s (%s)", db$display_id,
                                          db$db_display_name)
                updateSelectizeInput(session, "selectedProtein", 
                                     choices=uniprot)
            }
        }
    })
    
    # Update protein links depending on chosen transcript and protein
    observe({
        ensembl <- input$selectedTranscript
        uniprot <- input$selectedProtein
        
        species  <- tolower(getSpecies())
        assembly <- getAssemblyVersion()
        grch37   <- assembly == "hg19"
        if (is.null(species) || is.null(assembly)) return(NULL)
        
        href <- paste0("http://", if(grch37) { "grch37." }, "ensembl.org/", 
                       species, "/Search/Results?q=", ensembl)
        ensemblLink <- tags$a("Ensembl", icon("external-link"), target="_blank",
                              class="btn btn-link", href=href)
        
        href <- paste0("http://www.uniprot.org/uniprot/", uniprot)
        links <- tagList(ensemblLink, 
                         tags$a("Uniprot", icon("external-link"), href=href,
                                target="_blank", class="btn btn-link"))
        output$proteinLink <- renderUI(links)
    })
    
    # Plot UniProt proteins
    observe({
        uniprot <- input$selectedProtein
        if (is.null(uniprot) || uniprot == "") {
            output$plotProtein <- renderHighchart(NULL)
        } else {
            hc <- tryCatch(plotProtein(uniprot), error=return)
            output$plotProtein <- renderHighchart({
                if (is(hc, "error"))
                    stop(safeError(hc$message))
                
                return(hc)
            })
        }
    })
    
    # Render relevant articles according to available gene
    output$articles <- renderUI({
        event <- getEvent()
        gene <- parseEvent(event)$gene[[1]]
        if (length(gene) > 1)
            gene <- input$selectedGene
        
        if (is.null(gene))
            return(NULL)
        else {
            category <- unlist(strsplit(getCategory(), " "))
            articles <- pubmedUI(gene, "cancer", category, top=3)
            return(articles)
        }
    })
}

attr(infoUI, "loader") <- "analysis"
attr(infoUI, "name") <- "Gene, transcript and protein information"
attr(infoServer, "loader") <- "analysis"