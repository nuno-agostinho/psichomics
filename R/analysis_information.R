## TODO: check if API is alive before querying data

#' Query the Ensembl REST API
#'
#' @param path Character: API path
#' @param query Character: API query
#' @param grch37 Boolean: query the Ensembl GRCh37 API? if \code{FALSE}, query
#' the most recent API
#'
#' @importFrom httr GET timeout
#' @importFrom jsonlite fromJSON
#'
#' @return Parsed response or \code{NULL} if no response
#' @keywords internal
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

#' Query information from Ensembl
#'
#' @param gene Character: gene
#' @param species Character: species (may be \code{NULL} for an Ensembl
#' identifier)
#' @param assembly Character: assembly version (may be NULL for an Ensembl
#' identifier)
#'
#' @family functions to retrieve external information
#' @return Information from Ensembl
#' @export
#'
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

#' @rdname queryEnsemblByGene
#'
#' @param event Character: alternative splicing event
#' @inheritParams parseSplicingEvent
#'
#' @export
#' @examples
#' event <- "SE_17_-_41251792_41249306_41249261_41246877_BRCA1"
#' queryEnsemblByEvent(event, species="human", assembly="hg19")
queryEnsemblByEvent <- function(event, species=NULL, assembly=NULL, data=NULL) {
    gene <- parseSplicingEvent(event, data=data)$gene[[1]]
    stopifnot(
        "Could not parse AS event to identify cognate gene"=!is.null(gene),
        "AS event has no associated gene"=gene != "Hypothetical")
    return(queryEnsemblByGene(gene, species, assembly))
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
#' @keywords internal
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

#' Query the PubMed REST API
#'
#' @param primary Character: primary search term
#' @param ... Character: other relevant search terms
#' @param top Numeric: number of articles to retrieve
#' @param field Character: field of interest where to look for terms
#' (\code{abstract} by default)
#' @param sort Character: sort by a given parameter (\code{relevance} by
#' default)
#'
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#'
#' @return Parsed response
#' @keywords internal
#'
#' @examples
#' psichomics:::queryPubMed("BRCA1", "cancer", "adrenocortical carcinoma")
queryPubMed <- function(primary, ..., top=3, field="abstract",
                        sort="relevance") {
    args  <- unlist(list(...))
    if (!is.null(args)) {
        terms <- sprintf("%s AND (%s)", primary, paste(args, collapse=" OR "))
    } else {
        terms <- primary
    }
    url <- "https://eutils.ncbi.nlm.nih.gov"
    query <- list(db="pubmed", term=terms, retmax=top, tool="psichomics",
                  field=field, sort=sort,
                  email="nunodanielagostinho@gmail.com", retmode="json")
    resp <- GET(url, path="entrez/eutils/esearch.fcgi", query=query)
    warn_for_status(resp)
    search <- content(resp, "text", encoding = "UTF8")
    search <- fromJSON(search)[[2]]

    # Get summary information on the articles
    ids <- paste(search$idlist, collapse="+")
    query <- list(db="pubmed", tool="psichomics", id=ids,
                  email="nunodanielagostinho@gmail.com", retmode="json")
    resp <- GET(url, path="entrez/eutils/esummary.fcgi", query=query)
    warn_for_status(resp)
    metadata <- content(resp, "text", encoding = "UTF8")
    metadata <- fromJSON(metadata)[[2]][-1]
    return(c(search=list(search), metadata))
}

#' Convert from Ensembl to UniProt identifier
#'
#' @param protein Character: Ensembl identifier
#'
#' @family functions to retrieve external information
#' @return UniProt protein identifier
#' @export
#'
#' @examples
#' gene <- "ENSG00000173262"
#' ensemblToUniprot(gene)
#'
#' protein <- "ENSP00000445929"
#' ensemblToUniprot(protein)
ensemblToUniprot <- function(protein) {
    if(length(protein) != 1) stop("Only pass one Ensembl identifier")

    external <- queryEnsembl(paste0("xrefs/id/", protein),
                             list("content-type"="application/json"),
                             grch37=TRUE)

    if (is.null(external)) return(NULL)

    db <- external[grepl("Uniprot", external$dbname), ]
    uniprot <- db$primary_id
    names(uniprot) <- sprintf("%s (%s)", db$display_id, db$db_display_name)
    return(uniprot)
}

#' Return the interface of relevant PubMed articles for a given gene
#'
#' @param ns Namespace function
#' @param gene Character: gene
#' @inheritDotParams queryPubMed -primary
#'
#' @return HTML interface of relevant PubMed articles
#' @keywords internal
pubmedUI <- function(ns, gene, ...) {
    terms <- unlist(as.list(...))
    if (!is.null(gene)) gene <- setNames("gene", gene)
    terms <- c(gene, setNames(terms, terms))

    selectTerms <- selectizeInput(
        ns("articleTerms"), label=NULL, choices=terms, selected=terms,
        multiple=TRUE, width="auto", options=list(
            create=TRUE, createOnBlur=TRUE, persist=TRUE,
            plugins=list('remove_button'), placeholder="Add keywords..."))
    selectTerms[[2]]$style <- paste(selectTerms[[2]]$style, "margin-bottom: 0;")

    div(class="panel panel-default",
        div(class="panel-heading",
            tags$b("Relevant PubMed articles",
                   uiOutput(ns("articleSearch"), inline=TRUE))),
        div(class="list-group", tags$li(class="list-group-item", selectTerms),
            uiOutput(ns("articleList"))))
}

#' @rdname appUI
#' @importFrom shiny uiOutput
#' @importFrom R.utils capitalize
infoUI <- function(id) {
    ns <- NS(id)

    renderSelectize <- function(option_create) {
        render <- sprintf(
            "{ option_create: function (data, escape) {
                   return '<div class=\"create\">' + %s + '</div>'; }}",
            option_create)
        return(I(render))
    }

    species <- c(paste("human", c("hg19", "hg38")),
                 paste("mouse", c("mm9", "mm10")),
                 "rat rn6",
                 "zebrafish danRer11",
                 "Drosophila_melanogaster dm6",
                 "Caenorhabditis_elegans ce11",
                 "Saccharomyces_cerevisiae sacCer3")
    #setNames(paste(names(species), species),
    #         capitalize(sprintf("%s (%s assembly)", names(species), species)))

    onFocus <- I('function() { this.clear(); }')
    renderSpecies <- renderSelectize(
        option_create=paste(
            "'Search in <strong>' + escape(data.input).split(\" \")[0] +",
            "' (' + escape(data.input).split(\" \")[1] + ' assembly)' +",
            "'</strong>&hellip;'"))
    selectSpecies <- selectizeInput(
        ns("selectedSpecies"), "Species", width="100%",
        choices=species, selected=species[[1]],
        options=list(placeholder="Search for a species...", highlight=FALSE,
                     create=TRUE, onFocus=onFocus, render=renderSpecies))

    renderGene <- renderSelectize(
        option_create=paste(
            "'Search for <strong>' + escape(data.input) + '</strong>&hellip;'"))
    selectGene <- selectizeGeneInput(
        ns("selectedGene"), create=TRUE, createOnBlur=TRUE, render=renderGene,
        placeholder="Search for a gene symbol...")
    tagList(
        fluidRow(
            column(6, fluidRow(column(6, selectGene), column(6, selectSpecies)),
                   uiOutput(ns("genetic"))),
            column(6, pubmedUI(ns, gene=NULL, "cancer"))),
        uiOutput(ns("info")))
}

#' Parse XML from UniProt REST service
#'
#' @param xml response from UniProt
#'
#' @importFrom XML xmlTreeParse xmlRoot xmlAttrs xmlToList xmlName xmlChildren
#' @importFrom plyr ldply
#' @return List containing protein length and data frame of protein features
#' @keywords internal
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
        for (col in seq(ncol(feature))) {
            feature[[col]] <- as.character(feature[[col]])
        }
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
#' @family functions to retrieve external information
#' @return \code{highcharter} object
#' @export
#' @examples
#' protein <- "P38398"
#' plotProtein(protein)
#'
#' transcript <- "ENST00000488540"
#' plotProtein(transcript)
plotProtein <- function(molecule) {
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

    featureList <- NULL
    if (nrow(feature) == 0) {
        stop("No annotated domains were found for this protein in UniProt")
    }

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
#' @keywords internal
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
            extended, div(id="container", style="height: 100vh;"))
    }

    browsable(tagList(
        container,
        tags$script(sprintf("Highcharts.chart('container', %s)", hc))))
}

plotBand <- function(colour, from, to, gradient=NULL, text=NULL) {
    coords <- sort(c(from, to))
    from   <- coords[[1]]
    to     <- coords[[2]]

    if (!is.null(text))
        label <- list(text=text, y=10, style=list(fontWeight="bold"))
    else
        label <- list(y=10)

    if (is.null(gradient)) {
        res <- list(color=colour, from=from, to=to, label=label)
    } else {
        noColour <- "rgba(255, 255, 255, 0)"
        if (gradient == "invert") {
            firstColour <- colour
            lastColour  <- noColour
        } else {
            firstColour <- noColour
            lastColour  <- colour
        }

        res <- list(
            color=list(
                linearGradient=list(x1=1, x2=0, y1=1, y2=1),
                stops=list(c(0, firstColour), c(1, lastColour))),
            from=from, to=to, label=label)
    }
    return(res)
}

plotASeventRegion <- function(hc, event, data=NULL) {
    parsed <- parseSplicingEvent(event, coords=TRUE, data=data)[1, , drop=FALSE]
    if (is.null(parsed)) return(NULL)

    con1 <- sort(parsed$constitutive1[[1]])
    alt1 <- sort(parsed$alternative1[[1]])
    alt2 <- sort(parsed$alternative2[[1]])
    con2 <- sort(parsed$constitutive2[[1]])
    if (length(c(con1, alt1, alt2, con2)) == 0) return(NULL)
    type <- parsed$type

    if (type %in% c("AFE exon", "ALE exon")) type <- gsub(" exon", "", type)
    pretty <- prettifyEventType(type)

    if (type %in% c("MXE", "A3SS", "A5SS", "AFE", "ALE")) {
        text <- paste(pretty, "(alternative regions in orange and blue)")
    } else if (type %in% c("SE")) {
        text <- paste(pretty, "(alternative region in orange)")
    }

    orange <- "#ffb153cc"
    blue   <- "#7db6eccc"
    grey   <- "#d3d3d388"
    orangeBand <- blueBand <- greyBand <- NULL

    if (type == "SE") {
        orangeBand <- plotBand(orange, alt1[[1]], alt1[[2]])
        greyBand   <- plotBand(grey,   con1,      con2, text=text)
    } else if (type == "MXE") {
        orangeBand <- plotBand(orange, alt1[[1]], alt1[[2]])
        blueBand   <- plotBand(blue,   alt2[[1]], alt2[[2]])
        greyBand   <- plotBand(grey,   con1,      con2, text=text)
    } else if (type %in% c("A3SS", "A5SS", "AFE", "ALE")) {
        # Shift a given position
        shift <- function(pos, FUN, by=200) { FUN(as.numeric(pos), by) }

        if (type == "A3SS")      greyPos <- range(con1, alt1, alt2)
        else if (type == "A5SS") greyPos <- range(alt1, alt2, con2)
        else if (type == "AFE")  greyPos <- range(alt1, alt2, con2)
        else if (type == "ALE")  greyPos <- range(con1, alt1, alt2)

        plusStrand      <- parsed$strand == "+"
        downstreamMinus <- type %in% c("A3SS", "ALE") && !plusStrand
        downstreamPlus  <- type %in% c("A3SS", "ALE") && plusStrand
        upstreamPlus    <- type %in% c("A5SS", "AFE") && plusStrand
        upstreamMinus   <- type %in% c("A5SS", "AFE") && !plusStrand

        shiftBy <- min(abs(alt1 - alt2), 200)
        if (downstreamPlus || upstreamMinus) {
            gradient   <- "normal"
            orangePos <- c(alt1, shift(alt1, `+`, by=shiftBy))
            bluePos   <- c(alt2, shift(alt2, `+`, by=shiftBy))
        } else if (upstreamPlus || downstreamMinus) {
            gradient   <- "invert"
            orangePos <- c(shift(alt1, `-`, by=shiftBy), alt1)
            bluePos   <- c(shift(alt2, `-`, by=shiftBy), alt2)
        }
        greyBand   <- plotBand(grey, greyPos[[1]], greyPos[[2]], text=text)
        orangeBand <- plotBand(orange, orangePos[[1]], orangePos[[2]], gradient)
        blueBand   <- plotBand(blue, bluePos[[1]], bluePos[[2]], gradient)
    }
    # Ordered by stack order starting with the gray as background
    hc <- hc_xAxis(hc, plotBands=list(greyBand, blueBand, orangeBand))
    return(hc)
}

#' Plot transcripts
#'
#' @param info Information retrieved from Ensembl
#' @param eventPosition Numeric: coordinates of the alternative splicing event
#' (ignored if \code{event} is set)
#' @param event Character: identifier of the alternative splicing event to plot
#' @param eventData Object containing event information to be parsed
#' @param shiny Boolean: is the function running in a Shiny session?
#'
#' @importFrom highcharter highchart hc_chart hc_title hc_legend hc_xAxis
#' hc_yAxis hc_plotOptions hc_tooltip hc_series
#'
#' @family functions to retrieve external information
#' @inherit psichomics return
#' @export
#'
#' @examples
#' event <- "SE_12_-_7985318_7984360_7984200_7982602_SLC2A14"
#' info  <- queryEnsemblByEvent(event, species="human", assembly="hg19")
#' \dontrun{
#' plotTranscripts(info, event=event)
#' }
plotTranscripts <- function(info, eventPosition=NULL, event=NULL,
                            eventData=NULL, shiny=FALSE) {
    data <- list()
    for (i in seq(nrow(info$Transcript))) {
        transcript <- info$Transcript[i, ]
        name    <- transcript$id
        display <- transcript$display_name
        strand  <- ifelse(transcript$strand == 1, "+", "-")
        chr     <- transcript$seq_region_name
        start   <- transcript$start
        end     <- transcript$end
        biotype <- gsub("_", " ", transcript$biotype)

        # Prepare exons
        elements <- list()
        exons <- transcript$Exon[[1]]
        for (k in seq(nrow(exons))) {
            exon  <- exons[k, ]
            start <- exon$start
            end   <- exon$end
            len   <- abs(end - start)
            elements <- c(elements, list(list(
                name="exon", x=start, x2=end, y=0, width=20, length=len)))
        }

        if (nrow(exons) > 1) {
            # Prepare introns
            introns <- NULL
            introns$start <- head(sort(exons$end), length(exons$end) - 1)
            introns$end   <- sort(exons$start)[-1]
            for (j in seq(length(introns$start))) {
                start <- introns$start[[j]]
                end   <- introns$end[[j]]
                len   <- abs(end - start)
                elements <- c(elements, list(list(
                    name="intron", x=start, x2=end, y=0, width=5, length=len)))
            }
        }
        data <- c(data, list(list(name=name, borderRadius=0, pointWidth=10,
                                  display=display, strand=strand, chr=chr,
                                  biotype=biotype, data=elements)))
    }

    # Plot transcripts
    hc <- highchart() %>%
        hc_chart(type="xrange", zoomType="x") %>%
        hc_title(text="") %>%
        hc_legend(enabled=FALSE) %>%
        hc_xAxis(title=list(text="Position (nucleotides)"), showFirstLabel=TRUE,
                 showLastLabel=TRUE) %>%
        hc_yAxis(title=list(text=""), visible=FALSE) %>%
        hc_plotOptions(series=list(borderWidth=0.5)) %>%
        hc_tooltip(followPointer=TRUE)
    hc <- do.call("hc_series", c(list(hc), data))

    if (!is.null(event)) {
        plotRegion <- plotASeventRegion(hc, event, eventData)
        if (!is.null(plotRegion)) hc <- plotRegion
    } else if (!is.null(eventPosition)) {
        # Draw region only if splicing event position is provided
        eventStart <- eventPosition[1]
        eventEnd   <- eventPosition[2]
        hc <- hc_xAxis(hc, plotBands=list(
            color="#7cb5ec50", from=eventStart,
            to=eventEnd, label=list(text="Splicing Event", y=10,
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
#' @param species Character: species name
#' @param assembly Character: assembly version
#' @param grch37 Boolean: use version GRCh37 of the genome?
#' @param eventDiagram Diagram of selected alternative splicing event
#'
#' @importFrom shiny renderUI h2 h3 plotOutput
#'
#' @return HTML elements to render gene, protein and transcript annotation
#' @keywords internal
renderGeneticInfo <- function(output, info, species=NULL, assembly=NULL,
                              grch37=FALSE, eventDiagram=NULL, gene=NULL) {
    start <- as.numeric(info$start)
    end   <- as.numeric(info$end)

    item <- function(dt, dd, dtWidth="width: 80px;",
                     ddMargin="margin-left: 100px;",
                     ddStyle=NULL) {
        tagList(tags$dt(style=dtWidth, dt),
                tags$dd(style=ddMargin, style=ddStyle, dd))
    }

    unavailableItem <- function(dt, dd, ...) {
        info <- tagList(icon("times-circle"), dd)
        return(item(dt, info, ddStyle="color: gray;"))
    }

    if (is.null(info)) {
        msg <- "Information from Ensembl currently unavailable"
        genomicLocation <- unavailableItem("Location", msg)
        description     <- unavailableItem("Description", msg)
    } else {
        genomicLocation <- item(
            "Location", sprintf(
                "Chromosome %s: %s-%s (%s strand)",
                info$seq_region_name,
                format(start, big.mark=",", scientific=FALSE),
                format(end, big.mark=",", scientific=FALSE),
                ifelse(info$strand == -1, "reverse", "forward")))
        description <- item("Description",
                            sprintf("%s (%s)", info$description, info$biotype))
    }
    links <- prepareExternalLinks(info, species, assembly, grch37, gene)

    tagList(
        h2(style="margin-top: 0px;", info$display_name, tags$small(info$id)),
        tags$dl(class="dl-horizontal", style="margin-bottom: 0px;",
                genomicLocation,
                description,
                item("AS diagram", eventDiagram),
                item("Links",
                     tags$ul(class="list-inline", lapply(links, tags$li)))))
}

#' Return the interface to display an article
#'
#' @param article PubMed article
#'
#' @importFrom shiny tags h5
#'
#' @return HTML to render an article's interface
#' @keywords internal
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

    decodeHTMLentities <- function(char) {
        char <- gsub("&lt;", "<", char, fixed=TRUE)
        char <- gsub("&gt;", ">", char, fixed=TRUE)
        return(HTML(char))
    }

    tags$a(href=paste0("http://pubmed.gov/", pmid), target="_blank",
           class="list-group-item", h5(class="list-group-item-heading",
                                       decodeHTMLentities(article$title),
                                       tags$small(description)))
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
#' @keywords internal
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
    grch37      <- if (assembly == "hg19") "grch37." else ""
    href        <- sprintf("http://%sensembl.org/%s/Transcript/Summary?t=%s",
                           grch37, species, transcript)
    ensemblLink <- tags$a("Ensembl", icon("external-link"), href=href,
                          target="_blank")

    href        <- paste0("http://www.uniprot.org/uniprot/?query=", transcript)
    uniprotLink <- tags$a("UniProt", icon("external-link"), href=href,
                          target="_blank")

    links <- column(2, tags$label("External links"),
                    tags$ul(class="list-inline",
                            tags$li(style="padding-top: 7px;", ensemblLink),
                            tags$li(style="padding-top: 7px;", uniprotLink)))
    if (!is.null(protein)) {
        res <- tagList(name, role, links)
    } else {
        res <- links
    }
    return(res)
}

prepareExternalLinks <- function(info, species, assembly, grch37, gene) {
    linkTo <- function(title, href) {
        tags$a(title, icon("external-link"), target="_blank", href=href)
    }

    url <- list()
    if (!is.null(info)) {
        gene  <- info$id
        chr   <- info$seq_region_name
        start <- as.numeric(info$start)
        end   <- as.numeric(info$end)

        ucscPos <- sprintf("chr%s:%s-%s", chr, start, end)
        url$humanProteinAtlas <- sprintf(
            "http://www.proteinatlas.org/%s/pathology", gene)
        url$vastdb <- sprintf("http://vastdb.crg.eu/wiki/Gene:%s@Genome:%s",
                              gene, assembly)
    } else {
        ucscPos <- gene
        url$humanProteinAtlas <- sprintf(
            "http://www.proteinatlas.org/search/%s", gene)
        url$vastdb <- paste0("http://vastdb.crg.eu/wiki/search?query=", gene)
    }
    url$ensembl   <- sprintf("http://%sensembl.org/%s/Gene/Summary?g=%s",
                             ifelse(grch37, "grch37.", ""), species, gene)
    url$ucsc      <- sprintf(
        "https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s",
        assembly, ucscPos)
    url$geneCards <- paste0(
        "http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene)

    isHuman <- !is.null(species) && species == "human"
    links <- tagList(
        if (!is.null(species) && species != "")
            linkTo("Ensembl", url$ensembl),
        linkTo("UCSC", url$ucsc),
        if (isHuman) linkTo("GeneCards", url$geneCards),
        if (isHuman) linkTo("Human Protein Atlas", url$humanProteinAtlas),
        linkTo("VastDB", url$vastdb))
    return(links)
}

parseSpeciesAssembly <- function(info) {
    speciesRegex <- "(.*) (.*)"
    species      <- tolower(gsub(speciesRegex, "\\1", info))
    assembly     <- tolower(gsub(speciesRegex, "\\2", info))
    if (is.null(species) || length(species) == 0 || species == "") return(NULL)
    if (is.null(assembly)) assembly <- ""
    return(list(species=species, assembly=assembly))
}

plotSplicingEventHelper <- function(ASevent, data=NULL) {
    diagram <- suppressWarnings(plotSplicingEvent(ASevent, data=data)[[1]])
    if (diagram != "") {
        res <- HTML(diagram)
    } else {
        res <- NULL
    }
    return(res)
}

#' @rdname appServer
#'
#' @importFrom highcharter highchart %>%
#' @importFrom shiny fixedRow safeError
#' @importFrom methods is
#' @importFrom shinyjs hide show runjs
infoServer <- function(input, output, session) {
    ns <- session$ns

    # Update species/assembly according to the ones used for loading datasets
    observe({
        species  <- getSpecies()
        assembly <- getAssemblyVersion()
        selected <- paste(species, assembly)

        if (is.null(selected) || selected == "" || length(selected) == 0) {
            selected <- "human hg19"
        }
        updateSelectizeInput(session, "selectedSpecies", selected=selected)
    })

    # Update gene according to selected splicing event
    observe({
        ASevent <- getASevent()
        if (!is.null(ASevent)) {
            gene <- parseSplicingEvent(ASevent)$gene[[1]]
        } else {
            gene <- NULL
        }
        # Update gene choices according to available genes
        updateSelectizeInput(session, "selectedGene", choices=getGenes(),
                             selected=gene, server=TRUE)
    })

    output$eventDiagram <- renderUI({
        ASevent <- getASevent()
        if (!is.null(ASevent)) {
            res <- plotSplicingEventHelper(ASevent)
            if (is.null(res)) res <- "Diagram not supported"
        } else {
            res <- "No alternative splicing event selected or available"
        }
        return(res)
    })

    reactiveQueryPubMed <- reactive({
        terms <- input$articleTerms
        terms <- gsub("^gene_", "", terms)
        pubmed <- queryPubMed(terms[1], terms[-1])
        return(pubmed)
    })

    observeEvent(input$articleTerms, {
        pubmed <- tryCatch(reactiveQueryPubMed(), error=return)
        if (is(pubmed, "error")) return(NULL)

        articles <- pubmed[-1]
        if (length(articles) >= 1) {
            articleList <- lapply(articles, articleUI)
        } else {
            articleList <- helpText(class="list-group-item",
                                    "No articles match your search terms")
        }
        output$articleList <- renderUI(articleList)

        search <- pubmed$search$querytranslation
        search <- gsub("[Abstract]", "[Title/Abstract]", search, fixed=TRUE)
        search <- paste0("http://www.ncbi.nlm.nih.gov/pubmed/?term=", search)
        link <- tags$a(href=search, target="_blank", class="pull-right",
                       "Show more articles", icon("external-link"))
        output$articleSearch <- renderUI(link)
    })

    # Update gene used in search
    observe({
        gene <- input$selectedGene
        if (gene == "") return(NULL)
        terms <- isolate(input$articleTerms)

        genePrefix <- "gene_"
        terms <- terms[!startsWith(terms, genePrefix)]
        selected <- c(setNames(paste0(genePrefix, gene), gene), terms)
        updateSelectizeInput(session, "articleTerms", choices=selected,
                             selected=selected)
    })

    observe({
        parsed   <- parseSpeciesAssembly(input$selectedSpecies)
        if (is.null(parsed)) return(NULL)
        species  <- parsed$species
        assembly <- parsed$assembly
        grch37   <- assembly == "hg19"
        gene     <- input$selectedGene

        if (gene == "") return(NULL)
        info <- tryCatch(queryEnsemblByGene(gene, species=species,
                                            assembly=assembly), error=return)

        if (is(info, "error") || is.null(info)) {
            hide("info")
        } else {
            show("info")
        }

        output$genetic <- renderUI({
            if (is(info, "error")) info <- NULL
            geneticInfo <- renderGeneticInfo(
                output, info, species, assembly, grch37,
                eventDiagram=uiOutput(ns("eventDiagram")), gene=gene)

            if (is.null(info)) {
                error <- errorDialog(
                    "No response from Ensembl", bigger=TRUE,
                    "Please select another gene or try again later.")
                res <- tagList(geneticInfo, error)
            } else {
                res <- geneticInfo
            }
            return(res)
        })

        output$info <- renderUI({
            tagList(
                h3("Transcripts"),
                uiOutput(ns("plotTranscripts")),
                h3("Protein domains"),
                uiOutput(ns("selectProtein")),
                uiOutput(ns("proteinError")),
                highchartOutput(ns("plotProtein"), height="200px"))
        })

        output$plotTranscripts <- renderUI({
            info <- queryEnsemblByGene(gene, species=species, assembly=assembly)
            plotTranscripts(info, event=isolate(getASevent()), shiny=TRUE)
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

    plotProteinReactive <- reactive({ plotProtein(input$selectedTranscript) })

    # Render UniProt protein domains and information
    observe({
        transcript <- input$selectedTranscript
        parsed     <- parseSpeciesAssembly(input$selectedSpecies)
        if (is.null(parsed)) return(NULL)
        species    <- parsed$species
        assembly   <- parsed$assembly

        if (is.null(transcript) || transcript == "") {
            output$proteinInfo  <- renderUI(NULL)
            output$proteinError <- renderUI(NULL)
            output$plotProtein  <- renderHighchart(NULL)
        } else {
            hc <- tryCatch(plotProteinReactive(), error=return)
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
                             "transcript not available")
            })
            output$plotProtein <- renderHighchart(hc)
        }
    })
}

attr(infoUI, "loader") <- "analysis"
attr(infoUI, "name") <- "Gene, transcript and protein information"
attr(infoServer, "loader") <- "analysis"
