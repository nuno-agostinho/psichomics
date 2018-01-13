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
#' return NULL as they are only used to modify the state of hidden elements
getHidden <- function() .hidden$elem

#' @rdname getHidden
#' @param val Value to attribute
setHidden <- function(val) .hidden$elem <- val

#' Get or set globally accessible elements
#' 
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' @param ... Arguments to identify a variable
#' @param sep Character to separate identifiers
#' 
#' @note Needs to be called inside a reactive function
#' 
#' @seealso \code{\link{getEvent}}, \code{\link{getClinicalMatchFrom}},
#' \code{\link{getGroups}} and \code{\link{getDifferentialAnalyses}}
#' 
#' @return Getters return globally accessible data, whereas setters return NULL 
#' as they are only used to modify the Shiny session's state
getGlobal <- function(category=getCategory(), ..., sep="_") {
    sharedData[[paste(category, ..., sep=sep)]]
}

#' @rdname getGlobal
#' @param value Value to attribute to an element
setGlobal <- function(category=getCategory(), ..., value, sep="_") {
    sharedData[[paste(category, ..., sep=sep)]] <- value
}

#' Get global data
#' @return Variable containing all data of interest
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

#' @inherit getGlobal
getASevents <- function() {
    psi <- getInclusionLevels()
    if (!is.null(psi)) {
        choices <- rownames(psi)
        names(choices) <- parseSplicingEvent(choices, char=TRUE)
        return( sort(choices) )
    }
}

#' @inherit getGlobal
getASevent <- reactive(sharedData$event)

#' @rdname getEvent
#' @param event Character: alternative splicing event
setASevent <- function(event) setGlobal("event", value=event)

#' @inherit getGlobal
getEvent <- getASevent

#' @rdname getEvent
setEvent <- setASevent

#' @inherit getGlobal
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

#' @rdname getEvent
getCategories <- reactive(names(getData()))

#' @rdname getEvent
getCategory <- reactive(sharedData$category)

#' @rdname getEvent
setCategory <- function(category) setGlobal("category", value=category)

#' @rdname getEvent
getCategoryData <- reactive(
    if(!is.null(getCategory())) getData()[[getCategory()]])

#' @rdname getEvent
getActiveDataset <- reactive(sharedData$activeDataset)

#' @rdname getEvent
#' @param dataset Character: dataset name
setActiveDataset <- function(dataset) setGlobal("activeDataset", value=dataset)

#' @rdname getEvent
#' @param attrs Character: name of attributes to retrieve (if NULL, the whole 
#' dataset is returned)
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

#' @rdname getEvent
getPatientId <- function() {
    clinical <- getClinicalData()
    if ( !is.null(clinical) ) {
        return( rownames(clinical) )
    } else {
        return(NULL)
    }
}

#' @rdname getEvent
getPatientAttributes <- function() {
    clinical <- getClinicalData()
    if ( !is.null(clinical) ) {
        patientAttrs <- colnames(clinical)
        attr(patientAttrs, "default") <- attr(clinical, "show")
        return(patientAttrs)
    } else {
        return(NULL)
    }
}

#' @rdname getEvent
getSampleInfo <- reactive(getCategoryData()[["Sample metadata"]])

#' @rdname getEvent
setSampleInfo <- function(value, category = getCategory())
    setDataTable("Sample metadata", value, category)

#' @rdname getEvent
getSampleId <- function() {
    sampleInfo <- getSampleInfo()
    if ( !is.null(sampleInfo) ) {
        return( rownames(sampleInfo) )
    } else {
        return(NULL)
    }
}

#' @rdname getEvent
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

#' @rdname getEvent
getJunctionQuantification <- function(category=getCategory()) {
    if (!is.null(category)) {
        data <- getData()[[category]]
        match <- sapply(data, attr, "dataType") == "Junction quantification"
        if (any(match)) return(data[match])
    }
}

#' @rdname getEvent
getGeneExpression <- function(category=getCategory()) {
    if (!is.null(category)) {
        data <- getData()[[category]]
        match <- sapply(data, attr, "dataType") == "Gene expression"
        if (any(match)) return(data[match])
    }
}

#' @rdname getEvent
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

#' @rdname getEvent
getInclusionLevels <- reactive(getCategoryData()[["Inclusion levels"]])

#' @rdname getEvent
#' @param incLevels Data frame or matrix: inclusion levels
setInclusionLevels <- function(incLevels, category=getCategory())
    setDataTable("Inclusion levels", incLevels, category)

#' @rdname getEvent
getPCA <- function(category=getCategory())
    getGlobal(category, "PCA")

#' @rdname getEvent
#' @param pca \code{prcomp} object (principal component analysis)
setPCA <- function(pca, category=getCategory())
    setGlobal(category, "PCA", value=pca)

#' @rdname getEvent
getICA <- function(category=getCategory())
    getGlobal(category, "ICA")

#' @rdname getEvent
#' @param ica Object containing independent component analysis
setICA <- function(ica, category=getCategory())
    setGlobal(category, "ICA", value=ica)

#' @rdname getEvent
getGroupIndependenceTesting <- function(category=getCategory())
    getGlobal(category, "groupIndependenceTesting")

#' @rdname getEvent
#' @param groupIndependenceTesting Object containing group independence testing
#' results
setGroupIndependenceTesting <- function(groupIndependenceTesting, 
                                        category=getCategory()) {
    setGlobal(category, "groupIndependenceTesting", 
              value=groupIndependenceTesting)
}

#' @rdname getEvent
getSpecies <- function(category=getCategory())
    getGlobal(category, "species")

#' @rdname getEvent
#' @param species Character: species
setSpecies <- function(species, category=getCategory())
    setGlobal(category, "species", value=species)

#' @rdname getEvent
getAssemblyVersion <- function(category=getCategory())
    getGlobal(category, "assemblyVersion")

#' @rdname getEvent
#' @param assembly Character: assembly version
setAssemblyVersion <- function(assembly, category=getCategory())
    setGlobal(category, "assemblyVersion", value=assembly)

#' @rdname getEvent
getAnnotationName <- function(category=getCategory())
    getGlobal(category, "annotName")

#' @rdname getEvent
#' @param annotName Character: annotation name
setAnnotationName <- function(annotName, category=getCategory())
    setGlobal(category, "annotName", value=annotName)

#' @rdname getEvent
getURLtoDownload <- function() getGlobal("URLtoDownload")

#' @rdname getEvent
#' @param url Character: URL links to download
setURLtoDownload <- function(url) setGlobal("URLtoDownload", value=url)

#' Get or set clinical matches from a given data type
#' @inherit getGlobal
#' @param dataset Character: data set name (e.g. "Junction quantification")
getClinicalMatchFrom <- function(dataset, category=getCategory())
    getGlobal(category, dataset, "clinicalMatch")

#' @rdname getClinicalMatchFrom
#' @param matches Vector of integers: clinical matches of dataset
setClinicalMatchFrom <- function(dataset, matches, category=getCategory())
    setGlobal(category, dataset, "clinicalMatch", value=matches)

#' Get or set groups
#' @inherit getGlobal
#' 
#' @param type Character: type of groups (either "Patients", "Samples", 
#' "ASevents" or "Genes")
#' @param complete Boolean: return all the information on groups (TRUE) or just 
#' the group names and respective indexes (FALSE)? FALSE by default
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
#' @inherit getGlobal
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
#' @inherit getGlobal
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
#' @inherit getGlobal
getDifferentialAnalyses <- function(category=getCategory())
    getGlobal(category, "differentialAnalyses")

#' @rdname getDifferentialAnalyses
#' @param differential Data frame or matrix: differential analyses table
setDifferentialAnalyses <- function(differential, category=getCategory())
    setGlobal(category, "differentialAnalyses", value=differential)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesFiltered <- function(category=getCategory())
    getGlobal(category, "differentialAnalysesFiltered")

#' @rdname getDifferentialAnalyses
setDifferentialAnalysesFiltered <- function(differential, 
                                            category=getCategory())
    setGlobal(category, "differentialAnalysesFiltered", value=differential)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesSurvival <- function(category=getCategory())
    getGlobal(category, "diffAnalysesSurv")

#' @rdname getDifferentialAnalyses
#' @param survival Data frame or matrix: differential analyses' survival data
setDifferentialAnalysesSurvival <- function(survival, category=getCategory())
    setGlobal(category, "diffAnalysesSurv", value=survival)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesResetPaging <- function(category=getCategory())
    getGlobal(category, "diffAnalysesResetPaging")

#' @rdname getDifferentialAnalyses
#' @param reset Character: reset paging of differential analyses table?
setDifferentialAnalysesResetPaging <- function(reset, category=getCategory())
    setGlobal(category, "diffAnalysesResetPaging", value=reset)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesColumns <- function(category=getCategory())
    getGlobal(category, "diffAnalysesCols")

#' @rdname getDifferentialAnalyses
#' @param columns Character: differential analyses' column names
setDifferentialAnalysesColumns <- function(columns, category=getCategory())
    setGlobal(category, "diffAnalysesCols", value=columns)