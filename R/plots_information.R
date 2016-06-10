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
    stop_for_status(resp)
    r <- content(resp, "text", encoding = "UTF8")
    return(fromJSON(r))
}

infoUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("info"))
}

#' @importFrom highcharter highchart %>%
infoServer <- function(input, output, session) {
    gene <- "BRCA2"
    query <- paste0("lookup/symbol/human/", gene)
    info <- queryEnsembl(path, list(expand=1))
    
    start <- as.numeric(info$start)
    end   <- as.numeric(info$end)
    
    output$info <- renderUI({
        tagList(
            h2(info$display_name, tags$small(info$id)),
            sprintf("Chromosome %s: %s-%s (%s strand)",
                    info$seq_region_name, 
                    format(start, big.mark=",", scientific=FALSE),
                    format(end, big.mark=",", scientific=FALSE),
                    ifelse(info$strand == -1,"reverse", "forward")),
            br(), sprintf("%s (%s)", info$description,
                          info$biotype), br(),
            tags$a("Ensembl", icon("external-link"), target="_blank",
                   class="btn btn-link",
                   href=paste0("http://grch37.ensembl.org/human/",
                               "Gene/Summary?g=", info$id)),
            tags$a("UCSC", icon("external-link"), target="_blank",
                   class="btn btn-link",
                   href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks",
                               "?db=hg19&position=chr", 
                               info$seq_region_name, ":", start, "-", end))
        )
    })
}

attr(infoUI, "loader") <- "plots"
attr(infoUI, "name") <- "Information"
attr(infoServer, "loader") <- "plots"