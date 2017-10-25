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



#' Query information from Ensembl by a given alternative splicing event
#' 
#' @param event Character: alternative splicing event identifier
#' @inheritDotParams queryEnsemblByGene -gene
#' 
#' @return Information from Ensembl
#' @export
#' @examples 
#' event <- c("SE_17_-_41251792_41249306_41249261_41246877_BRCA1")
#' queryEnsemblByEvent(event, species="human", assembly="hg19")
queryEnsemblByEvent <- function(event, ...) {
    gene <- parseEvent(event)$gene[[1]]
    if (gene == "Hypothetical")
        stop("This event has no associated gene")
    return(queryEnsemblByGene(gene, ...))
}

#' Query information from Ensembl by a given gene
#' 
#' @param gene Character: gene identifier
#' @param species Character: species (can be NULL when handling an Ensembl
#' identifier)
#' @param assembly Character: assembly version (can be NULL when handling an
#' Ensembl identifier)
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

#' Query the UniProt REST API
#'
#' @param molecule Character: protein or transcript to query
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
#' 
#' transcript <- "ENST00000488540"
#' format <- "xml"
#' psichomics:::queryUniprot(transcript, format)
queryUniprot <- function(molecule, format="xml") {
    url <- "http://www.uniprot.org"
    path <- paste0("uniprot/?query=", molecule, "&format=", format)
    resp <- GET(url, path=path)
    warn_for_status(resp)
    r <- content(resp, "text", encoding = "UTF8")
    return(r)
}

#' Query the Proteins REST API
#'
#' @param acession Character: UniProt accessions (up to 100 values)
#' @param ensembl Character: Ensembl gene, transcript or translation identifier
#' (up to 20 values)
#' @param gene Character: UniProt gene name (up to 20 values)
#' @param protein Character: UniProt protein name
#' @param taxid Character: organism taxon identifier (up to 20 values; Human
#' by default)
#' @param chromosome Character: chromosome (up to 20 values)
#' @param location Integer: genome location range (start and end)
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#'
#' @return Parsed response in XML format
#'
#' @examples
#' psichomics:::queryProteinsAPI("Q7TQA5")
#' psichomics:::queryProteinsAPI(gene="BRCA1")
#' 
#' ensembl <- c("ENSG00000012048", "ENSG00000139618")
#' psichomics:::queryProteinsAPI(ensembl=ensembl)
queryProteinsAPI <- function(accession=NULL, ensembl=NULL, gene=NULL, 
                             protein=NULL, taxid="9606", chromosome=NULL,
                             location=NULL) {
    # Prepare multiple values as expected by API
    multipleValues <- function(var, collapse=",")
        if (!is.null(var)) paste(var, collapse=collapse)
    
    accession  <- multipleValues(accession)
    ensembl    <- multipleValues(ensembl)
    gene       <- multipleValues(gene)
    protein    <- multipleValues(protein)
    taxid      <- multipleValues(taxid)
    chromosome <- multipleValues(chromosome)
    location   <- multipleValues(location, "-")
    
    url <- "https://www.ebi.ac.uk"
    path <- "proteins/api/coordinates"
    resp <- GET(url, path=path, query=list(
        size=-1, accession=accession, ensembl=ensembl, gene=gene, 
        protein=protein, taxid=taxid, chromosome=chromosome, location=location))
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

#' Convert an Ensembl identifier to the respective UniProt identifier
#' 
#' @param molecule Character: Ensembl identifier
#' 
#' @return UniProt protein identifier
#' @export
#' @examples 
#' gene <- "ENSG00000173262"
#' ensemblToUniprot(gene)
#' 
#' protein <- "ENSP00000445929"
#' ensemblToUniprot(protein)
ensemblToUniprot <- function(molecule) {
    if(length(molecule) != 1) stop("Only pass one Ensembl identifier")
    
    external <- queryEnsembl(paste0("xrefs/id/", molecule),
                             list("content-type"="application/json"),
                             grch37=TRUE)
    
    if (is.null(external)) return(NULL)
        
    db <- external[grepl("Uniprot", external$dbname), ]
    uniprot <- db$primary_id
    names(uniprot) <- sprintf("%s (%s)", db$display_id, db$db_display_name)
    return(uniprot)
}

#' @rdname appUI
#' @importFrom shiny uiOutput
infoUI <- function(id) {
    ns <- NS(id)
    
    placeholder <- c(
        "Type a gene symbol to retrieve that gene's information..."="")
    tagList(
        fluidRow(
            column(6, selectizeInput(
                ns("selectedGene"), NULL, choices=placeholder, width="100%",  
                options=list(
                    create=TRUE, createOnBlur=TRUE,
                    onFocus=I(paste0('function() { $("#', ns("selectedGene"), 
                                     '")[0].selectize.clear(); }')),
                    onChange=I(paste0('function(value) { $("#',
                                      ns("selectedGene"), 
                                      '")[0].selectize.blur(); }')))),
                   uiOutput(ns("genetic"))),
            column(6, uiOutput(ns("articles")))),
        uiOutput(ns("info")))
    
}

#' Interface when no information could be retrieved
#' @param output Shiny output
#' @param description Character: description of the message to show to the user
#' @inheritDotParams inlineDialog -description
#' 
#' @importFrom shiny renderUI h3 br tags
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
noinfo <- function(output, description=paste(
    "No information available for this gene."), ...) {
    output$info <- renderUI(
        errorDialog(description, style="width: 400px;", ...))
}

#' Parse XML from Protein's RESTful service
#'
#' @param xml response from Proteins API
#'
#' @importFrom XML xmlTreeParse xmlRoot xmlAttrs xmlToList xmlName
#' xmlValue
#' @importFrom plyr ldply
#' 
#' @return List containing protein length and data frame of protein features
parseProteinsXML <- function(xml) {
    doc <- xmlTreeParse(xml)
    root <- suppressWarnings( xmlRoot(doc) )
    
    parseSingleProteinXml <- function(index, root) {
        root <- root[[index]]
        # Extract protein name, length, function and features
        gnCoordinate   <- root[names(root) == "gnCoordinate"][[1]]
        featureNodes   <- gnCoordinate[names(gnCoordinate) == "feature"]
        length         <- nchar(xmlValue(root[["sequence"]]))
        proteinName    <- xmlValue(root[["name"]])
        recommendeName <- xmlValue(
            root[["protein"]][["recommendedName"]][["fullName"]])
        strand <- ifelse(
            xmlAttrs(gnCoordinate[[1]])[["reverse_strand"]] == "true",
            "Negative", "Positive")
        chr <- xmlAttrs(gnCoordinate[[1]])[["chromosome"]]
        
        # Prepare table of exons
        exons <- gnCoordinate[names(gnCoordinate) == "genomicLocation"][[1]]
        exonsList <- xmlToList(exons)
        exonsList[".attrs"] <- NULL
        exonsDF <- data.frame(t(matrix(unlist(exonsList), nrow=5)), 
                              stringsAsFactors=FALSE)
        colnames(exonsDF) <- c("start", "stop", "genomicStart", "genomicStop", 
                               "id")
        for (i in 1:4) exonsDF[ , i] <- as.numeric(exonsDF[ , i])
        
        # Convert list of XMLNodes to list of characters
        l <- lapply(featureNodes, function(feat) {
            attrs <- xmlAttrs(feat)
            
            location  <- feat[[match("location", names(feat))]]
            start <- as.numeric(xmlAttrs(location[[1]]))
            # If there's no stop position, simply sum 1 to the start position
            stop  <- tryCatch(as.numeric(xmlAttrs(location[[2]])),
                              error=function(e) start+1)
            
            genomeLocation  <- feat[[match("genomeLocation", names(feat))]]
            genomeStart <- as.numeric(xmlAttrs(genomeLocation[[1]]))
            # If there's no stop position, simply sum 1 to the start position
            genomeStop  <- tryCatch(as.numeric(xmlAttrs(genomeLocation[[2]])),
                                    error=function(e) start+1)
            description <- xmlValue(feat[[match("description", names(feat))]])
            
            # Get original and variant aminoacid
            variant <- match("variation", names(feat))
            if (!is.na(variant) && !is.null(variant)) {
                original  <- xmlToList(feat[[match("original", names(feat))]])
                variation <- xmlToList(feat[[variant]])
                variant <- sprintf("%s>%s: ", original, variation)
            } else {
                variant <- NULL
            }
            
            evidence <- feat[["evidence"]]
            evidenceCode <- evidenceDb <- evidenceId <- NULL
            if (!is.null(evidence)) {
                evidenceCode <- xmlAttrs(evidence)[["code"]]
                if (length(evidence) > 0) {
                    evidenceDb   <- xmlAttrs(evidence[[1]])[["type"]]
                    evidenceId   <- xmlAttrs(evidence[[1]])[["id"]]
                }
            }
            return(c(attrs, start=start, stop=stop, genomeStart=genomeStart, 
                     genomeStop=genomeStop, description=description,
                     variant=variant, evidenceCode=evidenceCode, 
                     evidenceDb=evidenceDb, evidenceId=evidenceId))
        })
        
        # Convert list of characters to data frame of characters
        feature <- ldply(l, rbind)
        if (ncol(feature) > 0) {
            for (col in 1:ncol(feature))
                feature[[col]]  <- as.character(feature[[col]])
            
            feature$start       <- as.numeric(feature$start)
            feature$stop        <- as.numeric(feature$stop)
            feature$genomeStart <- as.numeric(feature$genomeStart)
            feature$genomeStop  <- as.numeric(feature$genomeStop)
            
            feature$description <- gsub(".", "", feature$description, 
                                        fixed=TRUE)
        }
        return(list(name=proteinName, length=length, chr=chr, strand=strand, 
                    exons=exonsDF, feature=feature))
    }
    
    lapply(seq(root), parseSingleProteinXml, root)
}

#' Get genomic coordinates for protein features based on exon coordinates
#'
#' @param index Integer: index
#' @param features List of features
#' @param exons List of exons
#'
#' @return Parsed data frame with features repeated for each exon that spans 
#' that same feature
getExonsSpanningFeature <- function(index, features, exons) {
    feature <- features[index, ]
    # Get exons spanning the feature
    firstExon <- feature[["genomeStart"]] <= exons$genomicStop & 
        feature[["genomeStart"]] >= exons$genomicStart
    firstExon <- which(firstExon)
    lastExon <- feature[["genomeStop"]] <= exons$genomicStop & 
        feature[["genomeStop"]] >= exons$genomicStart
    lastExon <- which(lastExon)
    if (length(firstExon) == 0 || length(lastExon) == 0)
        return(NULL)
    
    # Create table based on exons spanning the feature
    exons <- exons[seq(firstExon, lastExon), ]
    table <- feature[rep(1, nrow(exons)), ]
    table$genomeStart[-1] <- exons$genomicStart[-1]
    table$genomeStop[1:(nrow(exons) - 1)] <- 
        exons$genomicStop[1:(nrow(exons) - 1)]
    return(table)
}

#' Parse protein features from an XML response from the Proteins API
#' 
#' @param proteinList Parsed list containing protein features and respective
#' transcript's exons
#' 
#' @importFrom plyr rbind.fill
#' 
#' @return Data frame with parsed protein features
parseProteinFeatures <- function(proteinList) {
    res <- c()
    for (protein in proteinList) {
        features <- protein$feature
        if (!is.null(features) && nrow(features) > 0) {
            all <- lapply(seq(nrow(features)), getExonsSpanningFeature,
                          features, protein$exons)
            all <- rbind.fill(all)
            all$name   <- protein$name
            all$length <- protein$length
            all$chr    <- protein$chr
            all$strand <- protein$strand
            elem <- setNames(list(all), protein$name)
            res <- c(res, elem)
        }
    }
    return(rbind.fill(res))
}

#' Plot protein features
#'
#' @param features Parsed protein features
#' @param shiny Boolean: is the function running in a Shiny session? FALSE by
#' default
#'
#' @return Plot of protein features
plotProteinFeatures <- function(features, shiny=FALSE) {
    hc <- highchart() %>% 
        hc_chart(type="proteinPlot", zoomType="x") %>% 
        hc_tooltip(followPointer=TRUE) %>% 
        hc_title(text="Protein domains") %>%
        hc_yAxis(visible=FALSE)
    
    # Decompose data frame into list of lists
    df <- features
    colnames(df)[match(c("type", "genomeStart", "genomeStop"),
                       colnames(df))] <- c("name", "x", "x2")
    df$width <- 10
    df$y <- match(df$name, unique(df$name))
    ll <- lapply(split(df, seq_along(df[,1])), as.list)
    names(ll) <- NULL
    ll <- split(ll, df$name)
    
    for (name in names(ll)) {
        hc <- hc_add_series(hc, name=name, pointWidth=20, data=ll[[name]])
    }
    plottableXranges(hc, shiny=shiny)
}

#' Plot protein domains based on a given gene
#' 
#' @param gene Character: gene symbol
#' 
#' @return Plot of protein domains
plotProteinDomainsFrom <- function(gene) {
    xml <- parseProteinsXML(queryProteinsAPI(gene=gene))
    features <- parseProteinFeatures(xml)
    plotProteinFeatures(features)
}

#' Parse XML from UniProt's RESTful service
#'
#' @param xml response from UniProt
#'
#' @importFrom XML xmlTreeParse xmlRoot xmlAttrs xmlToList xmlName xmlChildren
#' @importFrom plyr ldply
#' @return List containing protein length and data frame of protein features
parseUniprotXML <- function(xml) {
    doc <- xmlTreeParse(xml)
    root <- xmlRoot(doc)[[1]]
    
    # Extract protein name, length, function and features
    names         <- vapply(xmlChildren(root), xmlName, character(1))
    comments      <- root[names == "comment"]
    role          <- lapply(comments, xmlAttrs) == "function"
    role          <- tryCatch(toString(comments[role][[1]][[1]][[1]]),
                              error=function(cond) NULL)
    featureNodes  <- root[names == "feature"]
    length        <- as.numeric(xmlAttrs(
        root[names == "sequence"][[1]])["length"])
    proteinName   <- tryCatch(
        toString(root[names == "protein"][[1]][[1]][[1]][[1]]),
        error=function(cond) NULL)
    
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
    
    return(list(name=proteinName, length=length, role=role, feature=feature))
}

#' Plot protein features
#'
#' @param molecule Character: UniProt protein or Ensembl transcript identifier
#'
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip
#' hc_add_series
#'
#' @return \code{highcharter} object
#' @export
#' @examples
#' \dontrun{
#' protein <- "P38398"
#' plotProtein(protein)
#' 
#' transcript <- "ENST00000488540"
#' plotProtein(transcript)
#' }
plotProtein <- function(molecule) {
    display("Retrieving protein annotation from UniProt...")
    xml     <- queryUniprot(molecule, "xml")
    if (xml == "") return(NULL)
    parsed  <- parseUniprotXML(xml)
    name    <- parsed$name
    length  <- parsed$length
    role    <- parsed$role
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
    
    display("Plotting protein domains...")
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
    
    attr(hc, "protein") <- parsed
    return(hc)
}

#' HTML code to plot a X-ranges series
#' 
#' @param hc \code{highcharter} object
#' @inheritParams plotTranscripts
#' 
#' @importFrom shiny tagList tags includeScript div
#' @importFrom htmltools browsable
#' @importFrom jsonlite toJSON
#' 
#' @return HTML elements
plottableXranges <- function(hc, shiny=FALSE) {
    hc <- toJSON(hc$x$hc_opts, auto_unbox=TRUE)
    hc <- gsub('"---|---"', "", hc)
    
    extended <- includeScript(insideFile("shiny", "www", "highcharts.ext.js"))
    
    if (shiny) {
        # No need to load Highcharts in Shiny
        container <- tagList(extended, div(id="container"))
    } else {
        container <- tagList(
            tags$script(src="https://code.highcharts.com/highcharts.js"),
            tags$script(src="https://code.highcharts.com/modules/exporting.js"),
            extended,
            div(id="container", style="height: 100vh;"))
    }
    
    browsable(tagList(
        container, 
        tags$script(sprintf("Highcharts.chart('container', %s)", hc))))
}

#' Plot transcripts
#' 
#' @param info Information retrieved from Ensembl
#' @param eventPosition Numeric: coordinates of the alternative splicing event;
#' NULL by default
#' @param shiny Boolean: is the function running in a Shiny session? FALSE by
#' default
#' 
#' @importFrom highcharter highchart hc_chart hc_title hc_legend hc_xAxis
#' hc_yAxis hc_plotOptions hc_tooltip hc_series %>%
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
#' @export
#' 
#' @examples
#' event <- "SE_12_-_7985318_7984360_7984200_7982602_SLC2A14"
#' info  <- queryEnsemblByEvent(event, species="human", assembly="hg19")
#' pos   <- parseSplicingEvent(event)$pos[[1]]
#' \dontrun{
#' plotTranscripts(info, pos)
#' }
plotTranscripts <- function(info, eventPosition=NULL, shiny=FALSE) {
    data <- list()
    for (i in 1:nrow(info$Transcript)) {
        transcript <- info$Transcript[i, ]
        name    <- transcript$id
        display <- transcript$display_name
        strand  <- ifelse(transcript$strand == 1, "forward", "reverse")
        chr     <- transcript$seq_region_name
        start   <- transcript$start
        end     <- transcript$end
        biotype <- gsub("_", " ", transcript$biotype)
        
        # Prepare exons
        elements <- list()
        exons <- transcript$Exon[[1]]
        for (k in 1:nrow(exons)) {
            exon  <- exons[k, ]
            start <- exon$start
            end   <- exon$end
            elements <- c(elements,
                          list(list(name="exon", x=start, x2=end, y=0,
                                    width=20)))
        }
        
        if (nrow(exons) > 1) {
            # Prepare introns
            introns <- NULL
            introns$start <- head(sort(exons$end), length(exons$end) - 1)
            introns$end   <- sort(exons$start)[-1]
            for (j in 1:length(introns$start)) {
                start <- introns$start[[j]]
                end   <- introns$end[[j]]
                elements <- c(elements,
                              list(list(name="intron", x=start, x2=end, y=0,
                                        width=5)))
            }
        }
        data <- c(data, list(list(name=name, borderRadius=0, pointWidth=10,
                                  display=display, strand=strand, chr=chr, 
                                  biotype=biotype, data=elements)))
    }
    
    # Plot transcripts
    hc <- highchart() %>%
        hc_chart(type="transcriptPlot", zoomType="x") %>%
        hc_title(text="") %>%
        hc_legend(enabled=FALSE) %>%
        hc_xAxis(title=list(text="Position (nucleotides)"), showFirstLabel=TRUE,
                 showLastLabel=TRUE) %>%
        hc_yAxis(title=list(text=""), visible=FALSE) %>%
        hc_plotOptions(series=list(borderWidth=0.5)) %>%
        hc_tooltip(followPointer=TRUE)
    hc <- do.call("hc_series", c(list(hc), data))
    
    if (!is.null(eventPosition)) {
        # Draw splicing event if event information is provided
        eventStart <- eventPosition[1]
        eventEnd   <- eventPosition[2]
        hc <- hc_xAxis(hc, plotBands=list(color="#7cb5ec50", from=eventStart, 
                                          to=eventEnd, label=list(
                                              text="Splicing Event", y=10,
                                              style=list(fontWeight="bold"))))
    }
    
    if (shiny)
        hc <- hc %>% hc_plotOptions(series=list(cursor="pointer", events=list(
            click="---function(e) { setTranscript(this.name); }---")))
    plottableXranges(hc, shiny)
}

#' Render genetic information
#' 
#' @param output Shiny output
#' @param ns Namespace function
#' @param info Information as retrieved from Ensembl
#' @param species Character: species name (NULL by default)
#' @param assembly Character: assembly version (NULL by default)
#' @param grch37 Boolean: use version GRCh37 of the genome? FALSE by default
#' 
#' @importFrom shiny renderUI h2 h3 plotOutput
#' @return HTML elements to render gene, protein and transcript annotation
renderGeneticInfo <- function(output, ns, info, species=NULL, assembly=NULL, 
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
            tags$a("Human Protein Atlas", icon("external-link"), 
                   target="_blank",
                   href=paste0("http://www.proteinatlas.org/", info$id, 
                               "/pathology"))
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
    
    output$genetic <- renderUI({ tagList(
        h2(style="margin-top: 0px;", info$display_name, tags$small(info$id)),
        tags$dl(
            class="dl-horizontal", style="margin-bottom: 0px;",
            tags$dt(style=dtWidth, "Species"),
            tags$dd(style=ddMargin, speciesInfo),
            tags$dt(style=dtWidth, "Location"),
            tags$dd(style=ddMargin,
                    sprintf("Chromosome %s: %s-%s (%s strand)",
                            info$seq_region_name,
                            format(start, big.mark=",", scientific=FALSE),
                            format(end, big.mark=",", scientific=FALSE),
                            ifelse(info$strand == -1, "reverse", "forward"))),
            tags$dt(style=dtWidth, "Description"),
            tags$dd(style=ddMargin,
                    sprintf("%s (%s)", info$description, info$biotype)),
            tags$dt(style=dtWidth, "Links"),
            tags$dd(style=ddMargin, 
                    tags$ul(class="list-inline", lapply(links, tags$li)) )))
    })
    
    tagList(
        h3("Transcripts"), 
        uiOutput(ns("plotTranscripts")),
        h3("Protein domains"),
        uiOutput(ns("selectProtein")),
        uiOutput(ns("proteinError")),
        highchartOutput(ns("plotProtein"), height="200px"))
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
    description <- sprintf("%s (%s)", authors, year)
    
    source <- article$source
    if (source != "") description <- sprintf("%s. %s", description, source)
    
    volume <- article$volume
    issue <- article$issue
    if (volume != "" && issue != "")
        description <- sprintf("%s, %s(%s)", description, volume, issue)
    else if (volume != "")
        description <- sprintf("%s, %s", description, volume)
    
    description <- sprintf("%s.", description)
    pmid <- article$articleids$value[1]
    tags$a(href=paste0("http://pubmed.gov/", pmid), target="_blank",
           class="list-group-item", h5(class="list-group-item-heading", 
                                       article$title, tags$small(description)))
}

#' Return the interface of relevant PubMed articles for a given gene
#' 
#' @param gene Character: gene
#' @inheritDotParams queryPubMed -primary
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

#' Render protein information
#' 
#' @param protein Character: protein identifier
#' @param transcript Character: Ensembl identifier of the protein's respective
#' transcript
#' @param species Character: species
#' @param assembly Character: assembly
#'
#' @return HTML elements
renderProteinInfo <- function(protein, transcript, species, assembly) {
    if (!is.null(protein)) {
        # Prepare protein name and length
        name <- sprintf("%s (%s aminoacids)", protein$name, protein$length)
        name <- column(2, tags$label("Protein name"), tags$ul(
            class="list-inline", tags$li(style="padding-top: 7px;", name)))
        
        # Prepare protein role
        if (is.null(protein$role) || protein$role == "")
            role <- helpText("No annotated function", style="margin: 0;")
        else
            role <- protein$role
        role <- column(5, tags$label("Protein function"),
                       tags$ul(class="list-inline",
                               tags$li(style="padding-top: 7px;", role)))
    }
    
    # Prepare external links
    grch37   <- assembly == "hg19"
    href <- paste0("http://", if(grch37) { "grch37." }, "ensembl.org/", 
                   species, "/Transcript/Summary?t=", transcript)
    ensemblLink <- tags$a("Ensembl", icon("external-link"), href=href,
                          target="_blank")
    
    href <- paste0("http://www.uniprot.org/uniprot/?query=", transcript)
    uniprotLink <- tags$a("UniProt", icon("external-link"), href=href,
                          target="_blank")
    links <- column(2, tags$label("External links"),
                    tags$ul(class="list-inline",
                            tags$li(style="padding-top: 7px;", ensemblLink),
                            tags$li(style="padding-top: 7px;", uniprotLink)))
    
    if (!is.null(protein))
        return(tagList(name, role, links))
    else
        return(links)
}

#' @rdname appServer
#' 
#' @importFrom highcharter highchart %>%
#' @importFrom shiny fixedRow safeError
#' @importFrom methods is
#' @importFrom shinyjs hide show
infoServer <- function(input, output, session) {
    ns <- session$ns
    
    observe({
        species  <- tolower(getSpecies())
        if (length(species) == 0) species <- "human"
        
        assembly <- getAssemblyVersion()
        if (is.null(assembly)) assembly <- "hg19"
        
        grch37   <- assembly == "hg19"
        gene <- input$selectedGene
            
        if (gene == "") return(NULL)
        info <- tryCatch(queryEnsemblByGene(gene, species=species, 
                                            assembly=assembly), error=return)
            
        if (is(info, "error") || is.null(info)) {
            output$info <- renderUI({
                title <- "No response from Ensembl"
                description <- "Please select another gene or try again later."
                fluidRow(column(
                    6, errorDialog(title, description, bigger=TRUE)))
            })
            
            hide("genetic")
            hide("articles")
            return(NULL)
        } else {
            show("genetic")
            show("articles")
        }
        
        output$info <- renderUI(
            renderGeneticInfo(output, ns, info, species, assembly, grch37))
        
        output$plotTranscripts <- renderUI({
            info <- queryEnsemblByGene(gene, species=species, assembly=assembly)
            plotTranscripts(info, shiny=TRUE)
        })
        
        # Render relevant articles according to available gene
        output$articles <- renderUI({
            number <- 3
            category <- getCategory()
            if (!is.null(category)) {
                category <- unlist(strsplit(getCategory(), " "))
                articles <- pubmedUI(gene, "cancer", category, top=number)
            } else {
                articles <- pubmedUI(gene, "cancer", top=number)
            }
            return(articles)
        })
        
        # Show NULL so it doesn't show previous results when loading
        output$selectProtein <- renderUI("Loading...")
        output$proteinError  <- renderUI(NULL)
        output$plotProtein   <- renderHighchart(NULL)
        
        output$selectProtein <- renderUI({
            transcripts <- info$Transcript$id
            tagList(
                fixedRow(
                    column(3, selectizeInput(ns("selectedTranscript"),
                                             label="Select a transcript", 
                                             choices=transcripts,
                                             width="auto")),
                    uiOutput(ns("proteinInfo"))))
        })
    })
    
    # Render UniProt protein domains and information
    observe({
        transcript <- input$selectedTranscript
        species    <- tolower(getSpecies())
        assembly   <- getAssemblyVersion()
        
        if (is.null(transcript) || transcript == "") {
            output$proteinInfo  <- renderUI(NULL)
            output$proteinError <- renderUI(NULL)
            output$plotProtein  <- renderHighchart(NULL)
        } else {
            hc <- tryCatch(plotProtein(transcript), error=return)
            output$proteinInfo  <- renderUI({
                if (is.null(species) || is.null(assembly)) return(NULL)
                
                if (!is(hc, "error") && !is.null(hc)) {
                    protein <- attr(hc, "protein")
                    renderProteinInfo(protein, transcript, species, assembly)
                } else {
                    renderProteinInfo(protein=NULL, transcript, species,
                                      assembly)
                }
            })
            output$proteinError <- renderUI({
                if (is(hc, "error"))
                    stop(safeError(hc$message))
                else if (is.null(hc))
                    helpText("Protein information from UniProt for this",
                             "transcript is not available.")
            })
            output$plotProtein <- renderHighchart(hc)
        }
    })
}

attr(infoUI, "loader") <- "analysis"
attr(infoUI, "name") <- "Gene, transcript and protein information"
attr(infoServer, "loader") <- "analysis"

#############################################################################
#' 
#' #' @importFrom httr GET content http_error config
#' #' @importFrom jsonlite fromJSON
#' #' 
#' #' @examples 
#' #' getLiftOverMap(from="hg38", to="hg19")
#' getLiftOverMap <- function(from="hg38", to="hg19") {
#'     fixCoordName <- function(var) switch(var, "hg38"="GRCh38", "hg19"="GRCh37")
#'     from <- fixCoordName(from)
#'     to   <- fixCoordName(to)
#'     
#'     map <- NULL
#'     for (chr in c(1:22, "X", "Y")) {
#'         base  <- "http://grch37.rest.ensembl.org"
#'         path  <- sprintf("/map/human/%s/%s:1..9999999999/%s", from, chr, to)
#'         query <- list("content-type"="application/json")
#'         resp <- tryCatch(GET(base, path=path, query=query, config=timeout(10)), 
#'                          error=return)
#'         if (is(resp, "error") || http_error(resp) || is.null(resp)) # Time out
#'             return(NULL) # for instance, time out
#'         r <- content(resp, "text", encoding = "UTF8")
#'         map$original <- rbind(map$original, fromJSON(r)[[1]]$original)
#'         map$mapped   <- rbind(map$mapped,   fromJSON(r)[[1]]$mapped)
#'     }
#'     return(map)
#' }
#' 
#' massLiftOver <- function(chr, start, end, map=NULL, from="hg38", to="hg19") {
#'     if (length(chr) != length(start) || length(chr) != length(end))
#'         stop("The arguments chr, start and end must have the same length")
#'     
#'     if (is.null(map)) map <- getLiftOverMap(from=from, to=to)
#'     
#'     # compWithOtherGroup <- function(x, FUN, y) Reduce("|", lapply(x, FUN, y))
#'     # 
#'     # # Get indexes based on original coordinates
#'     # original <- map$original
#'     # index <- original$seq_region_name %in% chr &
#'     #     compWithOtherGroup(start, "<=", original$end) &
#'     #     compWithOtherGroup(end, ">=", original$start)
#'     
#'     liftOver <- function(i, chr, start, end, map) {
#'         original <- map$original
#'         index <- which(original$seq_region_name %in% chr[i] &
#'                            start[i] <= original$end & end[i] >= original$start)
#'         mini <- min(index)
#'         maxi <- max(index)
#'         mapped <- map$mapped
#'         
#'         coord <- NULL
#'         coord$chr   <- mapped[mini, ]$seq_region_name
#'         coord$start <- mapped[mini, ]$start + abs(start[i] - 
#'                                                       original[mini, ]$start)
#'         coord$end   <-  mapped[maxi, ]$end - abs(end[i] - original[maxi, ]$end)
#'         return(coord)
#'     }
#'     lapply(seq(chr), liftOver, chr, start, end, map)
#' }
#' 
#' # library(httr)
#' # library(jsonlite)
#' # map <- getLiftOverMap(from="hg38", to="hg19")
#' # liftOver(chr, start, end, from="hg38", to="hg19", map)
#' 
#' library(httr)
#' library(jsonlite)
#' map <- getLiftOverMap(from="hg19", to="hg38")
#' 
#' chr <- c(1, 19, 10, 22, 12, 17, 22, 3, 2, 12)
#' start <- c(45270173, 35996888, 73048489, 31495882, 124857156, 35879174, 
#'            31495882, 56600768, 216236738, 123003598)
#' end   <- c(45270938, 35998360, 73050672, 31500302, 124862783, 35880641, 
#'            31500302, 56600978, 216238045, 123005932)
#' coords <- massLiftOver(chr, start, end, from="hg19", to="hg38", map)
#' 
#' chr   <- sapply(coords, "[[", "chr")
#' start <- sapply(coords, "[[", "start")
#' end   <- sapply(coords, "[[", "end")