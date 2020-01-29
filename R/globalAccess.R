## Functions to get and set globally accessible variables

#' @include utils.R 
NULL

# Global variable with all the data of a session
sharedData <- reactiveValues()

# Global variable to keep progress for CLI version
.hidden <- new.env()

#' Get or set hidden globally accessible elements
#' 
#' @return Getters return hidden globally accessible data, whereas setters 
#' return \code{NULL} as they are only used to modify the state of hidden
#' elements
#' @keywords internal
getHidden <- function() .hidden$elem

#' @rdname getHidden
#' @param val Value to attribute
setHidden <- function(val) .hidden$elem <- val

#' Get or set globally accessible elements
#' 
#' @param category Character: data category
#' @param ... Arguments to identify a variable
#' @param sep Character to separate identifiers
#' 
#' @note Needs to be called inside a reactive function
#' @family functions to get and set global variables
#' 
#' @return Getters return globally accessible data, whereas setters return 
#' \code{NULL} as they are only used to modify the Shiny session's state
#' @keywords internal
getGlobal <- function(category=getCategory(), ..., sep="_") {
    sharedData[[paste(category, ..., sep=sep)]]
}

#' @rdname getGlobal
#' @param value Value to attribute to an element
setGlobal <- function(category=getCategory(), ..., value, sep="_") {
    sharedData[[paste(category, ..., sep=sep)]] <- value
}

#' Get global data
#' 
#' @return Variable containing all data of interest
#' @keywords internal
getData <- reactive(sharedData$data)

#' @rdname getGlobal
#' @param data List of data frame or matrix to set as data
setData <- function(data) setGlobal("data", value=data)

#' @rdname getGlobal
#' @param name Character: data table name
setDataTable <- function(name, value, category=getCategory())
    sharedData$data[[category]][[name]] <- value

#' @rdname getGlobal
getAutoNavigation <- reactive(sharedData$autoNavigation)

#' @rdname getGlobal
#' @param auto Boolean: enable automatic navigation of browser history?
setAutoNavigation <- function(auto) setGlobal("autoNavigation", value=auto)

#' @rdname getGlobal
getCores <- reactive(sharedData$cores)

#' @rdname getGlobal
#' @param integer Integer: value of the setting
setCores <- function(integer) setGlobal("cores", value=integer)

#' @rdname getGlobal
#' @importFrom shiny isRunning
getSignificant <- function() {
    if (isRunning()) reactive(sharedData$significant)()
    else return(3)
}

#' @rdname getGlobal
setSignificant <- function(integer) setGlobal("significant", value=integer)

#' @rdname getGlobal
#' @importFrom shiny isRunning
getPrecision <- function() {
    if (isRunning()) reactive(sharedData$precision)()
    else return(3)
}

#' @rdname getGlobal
setPrecision <- function(integer) setGlobal("precision", value=integer)

#' @rdname getGlobal
#' @keywords internal
getASevents <- function() {
    psi <- getInclusionLevels()
    if (!is.null(psi)) {
        choices <- rownames(psi)
        names(choices) <- parseSplicingEvent(choices, char=TRUE)
        return( sort(choices) )
    }
}

#' @rdname getGlobal
#' @keywords internal
getASevent <- reactive(sharedData$event)

#' @rdname getGlobal
#' @param event Character: alternative splicing event
setASevent <- function(event) setGlobal("event", value=event)

#' @rdname getGlobal
getEvent <- getASevent

#' @rdname getGlobal
setEvent <- setASevent

#' @rdname getGlobal
getGenes <- function() {
    genes <- NULL
    
    # Retrieve genes from gene expression
    geneExpr <- getGeneExpression()
    if (!is.null(geneExpr)) {
        original <- sort(unique(unlist(lapply(geneExpr, rownames))))
        genes    <- gsub("\\|.*", "", original) # Process TCGA gene symbols
        unknown  <- genes == "?"
        genes[unknown] <- original[unknown]
    }
    
    # Retrieve genes based on AS events
    ASevents <- getASevents()
    if (!is.null(ASevents))
        genes <- c(unlist(parseSplicingEvent(ASevents)$gene), genes)
    
    if (!is.null(genes)) {
        genes   <- sort(unique(genes))
        # Show unknown genes last
        unknown <- gsub("\\|.*", "", genes) == "?"
        genes   <- c(genes[!unknown], genes[unknown])
    }
    return(genes)
}

#' Get curated, literature-based gene lists
#'
#' Available gene lists:
#' \itemize{
#'   \item{\strong{Sebestyen et al., 2016}: 1350 genes encoding RNA-binding 
#'   proteins, 167 of which are splicing factors}
#' }
#'
#' @family functions for data grouping
#' @return List of genes
#' @export
#'
#' @examples
#' getGeneList()
getGeneList <- function() {
    prepareCitation <- function(attr) {
        if (length(attr$Author) == 1)
            authors <- attr$Author[[1]]
        else
            authors <- paste(attr$Author[[1]], "et al.")
        
        title   <- attr$`Article Title`
        journal <- attr$Journal
        year    <- attr$Date
        volume  <- attr$Volume
        issue   <- attr$Issue
        pages   <- attr$Pages
        
        sprintf("%s (%s). %s. %s, %s(%s), %s", 
                authors, year, title, journal, volume, issue, pages)
    }

    # Sebestyen et al. 2016
    rbps     <- readFile("Sebestyen_et_al_2016.RDS")
    rbpSF    <- rbps$`RNA-binding proteins that are splicing factors`
    rbpNonSF <- rbps$`RNA-binding proteins that are not splicing factors`    
    sebestyen2016 <- list(
        "RNA-binding protein splicing factors"=rbpSF,
        "RNA-binding proteins"=sort(c(rbpSF, rbpNonSF)))
    attr(sebestyen2016, "citation") <- prepareCitation(attributes(rbps))
    
    res <- list("Sebestyen et al. 2016"=sebestyen2016)
    class(res) <- c("geneList", class(res))
    return(res)
}

#' Print gene list
#' 
#' @param object \code{geneList}
#' 
#' @return Print available gene lists
#' @keywords internal
print.geneList <- function(object) {
    for (set in names(object)) {
        cat(sprintf(set), fill=TRUE)
        for (item in names(object[[set]])) {
            ll <- object[[set]][[item]]
            sample <- 4
            genes <- paste(head(ll, n=sample), collapse=", ")
            if (length(ll) > sample) genes <- paste0(genes, ", ...")
            cat(sprintf("  -> %s [%s genes]: %s", 
                        item, length(ll), genes), fill=TRUE)
        }
        cat(fill=TRUE)
        cat("Source:", attr(object[[set]], "citation"), fill=TRUE)
        
        consoleWidth <- options("width")
        cat(paste(rep("=", consoleWidth), collapse=""), fill=TRUE)
    }
}

#' @rdname getGlobal
getCategories <- reactive(names(getData()))

#' @rdname getGlobal
getCategory <- reactive(sharedData$category)

#' @rdname getGlobal
setCategory <- function(category) setGlobal("category", value=category)

#' @rdname getGlobal
getCategoryData <- reactive(
    if(!is.null(getCategory())) getData()[[getCategory()]])

#' @rdname getGlobal
getActiveDataset <- reactive(sharedData$activeDataset)

#' @rdname getGlobal
#' @param dataset Character: dataset name
setActiveDataset <- function(dataset) setGlobal("activeDataset", value=dataset)

#' @rdname getGlobal
#' @param attrs Character: name of attributes to retrieve (if \code{NULL}, the
#' whole dataset is returned)
getClinicalData <- function(attrs=NULL) {
    clinical <- getCategoryData()[["Clinical data"]]
    attrs <- attrs[attrs != ""]
    if (!is.null(attrs)) {
        cols <- lapply(attrs, grep, colnames(clinical), fixed=TRUE)
        cols <- unique(unlist(cols))
        if (length(cols) > 0) {
            clinical <- clinical[ , cols, drop=FALSE]
        } else {
            clinical <- NULL
        }
    }
    return(clinical)
}

#' @rdname getGlobal
getSubjectId <- function() {
    clinical <- getClinicalData()
    if ( !is.null(clinical) ) {
        return( rownames(clinical) )
    } else {
        return(NULL)
    }
}

#' @rdname getGlobal
getSubjectAttributes <- function() {
    clinical <- getClinicalData()
    if ( !is.null(clinical) ) {
        subjectAttrs <- colnames(clinical)
        attr(subjectAttrs, "default") <- attr(clinical, "show")
        return(subjectAttrs)
    } else {
        return(NULL)
    }
}

#' @rdname getGlobal
getSampleInfo <- reactive(getCategoryData()[["Sample metadata"]])

#' @rdname getGlobal
setSampleInfo <- function(value, category = getCategory())
    setDataTable("Sample metadata", value, category)

#' @rdname getGlobal
getSampleId <- function() {
    sampleInfo <- getSampleInfo()
    if ( !is.null(sampleInfo) ) {
        return( rownames(sampleInfo) )
    } else {
        return(NULL)
    }
}

#' @rdname getGlobal
getSampleAttributes <- function() {
    sampleInfo <- getSampleInfo()
    if ( !is.null(sampleInfo) ) {
        sampleAttrs <- colnames(sampleInfo)
        attr(sampleAttrs, "default") <- attr(sampleInfo, "show")
        return(sampleAttrs)
    } else {
        return(NULL)
    }
}

#' @rdname getGlobal
getJunctionQuantification <- function(category=getCategory()) {
    if (!is.null(category)) {
        data <- getData()[[category]]
        match <- sapply(data, attr, "dataType") == "Junction quantification"
        if (any(match)) return(data[match])
    }
}

#' @rdname getGlobal
#' @param item Character: name of specific item to retrieve (if \code{NULL}, the
#' whole list is returned)
#' @param EList Boolean: return gene expression datasets as \code{EList} if
#' possible or as data frames?
getGeneExpression <- function(item=NULL, category=getCategory(), EList=FALSE) {
    if (!is.null(category)) {
        data <- getData()[[category]]
        match <- sapply(data, attr, "dataType") == "Gene expression"
        if (any(match)) {
            df <- data[match]
            if (!is.null(item)) {
                res <- df[[item]]
                if (!EList && is(res, "EList")) {
                    # Convert EList object to data frame
                    res <- inheritAttrs(data.frame(res), res)
                }
            } else if (!EList) {
                # Convert EList objects to data frames
                res <- lapply(data[match], function(i) {
                    if (is(i, "EList")) data.frame(i) else i
                })
            }
            return(res)
        }
    }
}

#' @rdname getGlobal
#' @param geneExpr Data frame or matrix: normalised gene expression
setNormalisedGeneExpression <- function(geneExpr, category=getCategory()) {
    ns  <- names(getData()[[category]])
    num <- gsub("Gene expression \\(normalised.*?([0-9]*).*\\)", "\\1", ns)
    num <- suppressWarnings(as.integer(num))
    
    if (any(!is.na(num))) {
        num <- max(num, na.rm=TRUE)
        num <- paste0(" ", num + 1)
    } else if ("Gene expression (normalised)" %in% ns) {
        num <- " 1"
    } else {
        num <- ""
    }
    ns <- sprintf("Gene expression (normalised%s)", num)
    setDataTable(ns, geneExpr, category)
}

#' @rdname getGlobal
getInclusionLevels <- reactive(getCategoryData()[["Inclusion levels"]])

#' @rdname getGlobal
#' @param incLevels Data frame or matrix: inclusion levels
setInclusionLevels <- function(incLevels, category=getCategory())
    setDataTable("Inclusion levels", incLevels, category)

#' @rdname getGlobal
getPCA <- function(category=getCategory())
    getGlobal(category, "PCA")

#' @rdname getGlobal
#' @param pca \code{prcomp} object (principal component analysis)
setPCA <- function(pca, category=getCategory())
    setGlobal(category, "PCA", value=pca)

#' @rdname getGlobal
getICA <- function(category=getCategory())
    getGlobal(category, "ICA")

#' @rdname getGlobal
#' @param ica Object containing independent component analysis
setICA <- function(ica, category=getCategory())
    setGlobal(category, "ICA", value=ica)

#' @rdname getGlobal
getCorrelation <- function(category=getCategory())
    getGlobal(category, "correlation")

#' @rdname getGlobal
#' @param correlation \code{prcomp} object (correlation analyses)
setCorrelation <- function(correlation, category=getCategory())
    setGlobal(category, "correlation", value=correlation)

#' @rdname getGlobal
getGroupIndependenceTesting <- function(category=getCategory())
    getGlobal(category, "groupIndependenceTesting")

#' @rdname getGlobal
#' @param groupIndependenceTesting Object containing group independence testing
#' results
setGroupIndependenceTesting <- function(groupIndependenceTesting, 
                                        category=getCategory()) {
    setGlobal(category, "groupIndependenceTesting", 
              value=groupIndependenceTesting)
}

#' @rdname getGlobal
getSpecies <- function(category=getCategory())
    getGlobal(category, "species")

#' @rdname getGlobal
#' @param species Character: species
setSpecies <- function(species, category=getCategory())
    setGlobal(category, "species", value=species)

#' @rdname getGlobal
getAssemblyVersion <- function(category=getCategory())
    getGlobal(category, "assemblyVersion")

#' @rdname getGlobal
#' @param assembly Character: assembly version
setAssemblyVersion <- function(assembly, category=getCategory())
    setGlobal(category, "assemblyVersion", value=assembly)

#' @rdname getGlobal
getAnnotationName <- function(category=getCategory())
    getGlobal(category, "annotName")

#' @rdname getGlobal
#' @param annotName Character: annotation name
setAnnotationName <- function(annotName, category=getCategory())
    setGlobal(category, "annotName", value=annotName)

#' @rdname getGlobal
getURLtoDownload <- function() getGlobal("URLtoDownload")

#' @rdname getGlobal
#' @param url Character: URL links to download
setURLtoDownload <- function(url) setGlobal("URLtoDownload", value=url)

#' Get or set clinical matches from a given data type
#' @inherit getGlobal
#' @param dataset Character: data set name
#' 
#' @family functions to get and set global variables
#' @keywords internal
getClinicalMatchFrom <- function(dataset, category=getCategory())
    getGlobal(category, dataset, "clinicalMatch")

#' @rdname getClinicalMatchFrom
#' @param matches Vector of integers: clinical matches of dataset
setClinicalMatchFrom <- function(dataset, matches, category=getCategory())
    setGlobal(category, dataset, "clinicalMatch", value=matches)

#' Get or set groups
#' 
#' @inherit getGlobal
#' 
#' @param type Character: type of groups (either \code{Patients}, 
#' \code{Samples}, \code{ASevents} or \code{Genes})
#' @param complete Boolean: return all the information on groups (\code{TRUE})
#' or just the group names and respective indexes (\code{FALSE})?
#' 
#' @family functions to get and set global variables
#' @keywords internal
getGroups <- function(type=c("Patients", "Samples", "ASevents", "Genes"), 
                      complete=FALSE, category=getCategory()) {
    type <- match.arg(type)
    if (type %in% c("Patients", "Samples") )
        groups <- getGlobal(category, "sampleGroups")
    else if (type %in% c("ASevents", "Genes"))
        groups <- getGlobal(category, "ASeventGroups")
    
    # Return all data if requested
    if (complete) return(groups)
    
    # Check if data of interest is available
    if (!type %in% colnames(groups)) return(NULL)
    
    # If available, return data of interest
    g <- groups[ , type, drop=TRUE]
    if (length(g) == 1) names(g) <- rownames(groups)
    
    # Return colour lookup table for groups
    if ("Colour" %in% colnames(groups)) {
        colour <- groups[ , "Colour", drop=TRUE]
        colour <- setNames(unlist(colour), names(colour))
        attr(g, "Colour") <- colour
    }
    return(g)
}

#' @rdname getGroups
#' @param groups Matrix: groups of dataset
setGroups <- function(type=c("Patients", "Samples", "ASevents", "Genes"), 
                      groups, category=getCategory()) {
    type <- match.arg(type)
    if (type %in% c("Patients", "Samples") )
        type <- "sampleGroups"
    else if (type %in% c("ASevents", "Genes"))
        type <- "ASeventGroups"
    
    setGlobal(category, type, value=groups)
}


# Plot points or regions --------------------------------------------------

#' Get or set points or regions for plots
#' 
#' @inherit getGlobal
#' 
#' @family functions to get and set global variables
#' @keywords internal
getHighlightedPoints <- function(id, category=getCategory())
    getGlobal(category, id, "highlighted")

#' @rdname getHighlightedPoints
#' @param events Integer: index of events
setHighlightedPoints <- function(id, events, category=getCategory())
    setGlobal(category, id, "highlighted", value=events)

#' @rdname getHighlightedPoints
#' @param id Character: identifier
getZoom <- function(id, category=getCategory())
    getGlobal(category, id, "zoom")

#' @rdname getHighlightedPoints
#' @param zoom Integer: range of X and Y coordinates for zooming
setZoom <- function(id, zoom, category=getCategory())
    setGlobal(category, id, "zoom", value=zoom)

#' @rdname getHighlightedPoints
getSelectedPoints <- function(id, category=getCategory())
    getGlobal(category, id, "selected")

#' @rdname getHighlightedPoints
setSelectedPoints <- function(id, events, category=getCategory())
    setGlobal(category, id, "selected", value=events)

#' @rdname getHighlightedPoints
getLabelledPoints <- function(id, category=getCategory())
    getGlobal(category, id, "labelled")

#' @rdname getHighlightedPoints
setLabelledPoints <- function(id, events, category=getCategory())
    setGlobal(category, id, "labelled", value=events)


# Differential expression --------------------------------------------------

#' Get or set differential expression' elements for a data category
#' 
#' @inherit getGlobal
#' 
#' @family functions to get and set global variables
#' @keywords internal
getDifferentialExpression <- function(category=getCategory())
    getGlobal(category, "differentialExpression")

#' @rdname getDifferentialExpression
#' @param differential Data frame or matrix: differential analyses table
setDifferentialExpression <- function(differential, category=getCategory())
    setGlobal(category, "differentialExpression", value=differential)

#' @rdname getDifferentialExpression
getDifferentialExpressionFiltered <- function(category=getCategory())
    getGlobal(category, "differentialExpressionFiltered")

#' @rdname getDifferentialExpression
setDifferentialExpressionFiltered <- function(differential, 
                                            category=getCategory())
    setGlobal(category, "differentialExpressionFiltered", value=differential)

#' @rdname getDifferentialExpression
getDifferentialExpressionSurvival <- function(category=getCategory())
    getGlobal(category, "differentialExpressionSurvival")

#' @rdname getDifferentialExpression
#' @param survival Data frame or matrix: differential analyses' survival data
setDifferentialExpressionSurvival <- function(survival, category=getCategory())
    setGlobal(category, "differentialExpressionSurvival", value=survival)

#' @rdname getDifferentialExpression
getDifferentialExpressionResetPaging <- function(category=getCategory())
    getGlobal(category, "diffExpressionResetPaging")

#' @rdname getDifferentialExpression
#' @param reset Character: reset paging of differential analyses table?
setDifferentialExpressionResetPaging <- function(reset, category=getCategory())
    setGlobal(category, "diffExpressionResetPaging", value=reset)

#' @rdname getDifferentialExpression
getDifferentialExpressionColumns <- function(category=getCategory())
    getGlobal(category, "diffExpressionCols")

#' @rdname getDifferentialExpression
#' @param columns Character: differential analyses' column names
setDifferentialExpressionColumns <- function(columns, category=getCategory())
    setGlobal(category, "diffExpressionCols", value=columns)


# Differential splicing ---------------------------------------------------

#' Get or set differential splicing' elements for a data category
#' 
#' @inherit getGlobal
#' 
#' @family functions to get and set global variables
#' @keywords internal
getDifferentialSplicing <- function(category=getCategory())
    getGlobal(category, "differentialSplicing")

#' @rdname getDifferentialSplicing
#' @param differential Data frame or matrix: differential analyses table
setDifferentialSplicing <- function(differential, category=getCategory())
    setGlobal(category, "differentialSplicing", value=differential)

#' @rdname getDifferentialSplicing
getDifferentialSplicingFiltered <- function(category=getCategory())
    getGlobal(category, "differentialSplicingFiltered")

#' @rdname getDifferentialSplicing
setDifferentialSplicingFiltered <- function(differential, 
                                            category=getCategory())
    setGlobal(category, "differentialSplicingFiltered", value=differential)

#' @rdname getDifferentialSplicing
getDifferentialSplicingSurvival <- function(category=getCategory())
    getGlobal(category, "diffSplicingSurv")

#' @rdname getDifferentialSplicing
#' @param survival Data frame or matrix: differential analyses' survival data
setDifferentialSplicingSurvival <- function(survival, category=getCategory())
    setGlobal(category, "diffSplicingSurv", value=survival)

#' @rdname getDifferentialSplicing
getDifferentialSplicingResetPaging <- function(category=getCategory())
    getGlobal(category, "diffSplicingResetPaging")

#' @rdname getDifferentialSplicing
#' @param reset Character: reset paging of differential analyses table?
setDifferentialSplicingResetPaging <- function(reset, category=getCategory())
    setGlobal(category, "diffSplicingResetPaging", value=reset)

#' @rdname getDifferentialSplicing
getDifferentialSplicingColumns <- function(category=getCategory())
    getGlobal(category, "diffSplicingCols")

#' @rdname getDifferentialSplicing
#' @param columns Character: differential analyses' column names
setDifferentialSplicingColumns <- function(columns, category=getCategory())
    setGlobal(category, "diffSplicingCols", value=columns)