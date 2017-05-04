## Functions to get and set globally accessible variables

#' @include utils.R 
NULL

# Global variable with all the data of a session
sharedData <- reactiveValues()

#' Get or set globally accessible elements
#' 
#' @param ... Arguments to identify a variable
#' @param sep Character to separate identifiers
#' @param category Character: data category (e.g. "Carcinoma 2016"); by default,
#' it uses the selected data category
#' 
#' @note Needs to be called inside a reactive function
#' 
#' @seealso \code{\link{psichomics:::getEvent}},
#' \code{\link{psichomics:::getClinicalMatchFrom}},
#' \code{\link{psichomics:::getGroupsFrom}} and
#' \code{\link{psichomics:::getDifferentialAnalyses}}
#' 
#' @return Getters return globally accessible data, whereas setters return NULL 
#' as they are only used to modify the Shiny session's state
getGlobal <- function(..., sep="_") sharedData[[paste(..., sep=sep)]]

#' @rdname getGlobal
#' @param value Value to attribute to an element
setGlobal <- function(..., value, sep="_") {
    sharedData[[paste(..., sep=sep)]] <- value
}

#' Get global data
#' @return Variable containing all data of interest
getData <- reactive(sharedData$data)

#' @rdname getGlobal
#' @param data Data frame or matrix to set as data
setData <- function(data) setGlobal("data", value=data)

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
getSignificant <- reactive(sharedData$significant)

#' @rdname getGlobal
setSignificant <- function(integer) setGlobal("significant", value=integer)

#' @rdname getGlobal
getPrecision <- reactive(sharedData$precision)

#' @rdname getGlobal
setPrecision <- function(integer) setGlobal("precision", value=integer)

#' Set or get central elements
#' @inherit getGlobal
getEvent <- reactive(sharedData$event)

#' @rdname getEvent
#' @param event Character: alternative splicing event
setEvent <- function(event) setGlobal("event", value=event)

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
getClinicalData <- reactive(getCategoryData()[["Clinical data"]])

#' @rdname getEvent
getJunctionQuantification <- function(category=getCategory()) {
    if (!is.null(category)) {
        data <- getData()[[category]]
        match <- sapply(data, attr, "dataType") == "Junction quantification"
        if (any(match)) return(data[match])
    }
}

#' @rdname getEvent
getInclusionLevels <- reactive(getCategoryData()[["Inclusion levels"]])

#' @rdname getEvent
#' @param incLevels Data frame or matrix: inclusion levels
setInclusionLevels <- function(incLevels, category=getCategory())
    sharedData$data[[category]][["Inclusion levels"]] <- incLevels

#' @rdname getEvent
getInclusionLevelsPCA <- function(category=getCategory())
    getGlobal(category, "inclusionLevelsPCA")

#' @rdname getEvent
#' @param pca \code{prcomp} object (PCA) of inclusion levels
setInclusionLevelsPCA <- function(pca, category=getCategory())
    setGlobal(category, "inclusionLevelsPCA", value=pca)

#' @rdname getEvent
getSampleInfo <- reactive(getCategoryData()[["Sample metadata"]])

#' @rdname getEvent
#' @param info Data frame or matrix: sample information
setSampleInfo <- function(info, category=getCategory())
    sharedData$data[[category]][["Sample metadata"]] <- value

#' @rdname getEvent
getPatientId <- function(category=getCategory())
    getGlobal(category, "patients")

#' @rdname getEvent
#' @param patients Character: identifier of patients
setPatientId <- function(patients, category=getCategory())
    setGlobal(category, "patients", value=patients)

#' @rdname getEvent
getSampleId <- function(category=getCategory())
    getGlobal(category, "samples")

#' @rdname getEvent
#' @param samples Character: identifier of samples
setSampleId <- function(samples, category=getCategory())
    setGlobal(category, "samples", value=samples)

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

#' Get or set groups from a given data type
#' @inherit getGlobal
#' 
#' @param dataset Character: data set name (e.g. "Clinical data")
#' @param complete Boolean: return all the information on groups (TRUE) or just 
#' the group names and respective indexes (FALSE)? FALSE by default
#' @param samples Boolean: show groups by samples (TRUE) or patients (FALSE)?
#' FALSE by default
getGroupsFrom <- function(dataset, category=getCategory(), complete=FALSE,
                          samples=FALSE) {
    groups <- getGlobal(category, dataset, "groups")
    
    # Return all data if requested
    if (complete) return(groups)
    
    if (samples)
        col <- "Samples"
    else
        col <- "Patients"
    
    # Check if data of interest is available
    if (!col %in% colnames(groups)) return(NULL)
    
    # If available, return data of interest
    g <- groups[ , col, drop=TRUE]
    if (length(g) == 1) names(g) <- rownames(groups)
    return(g)
}

#' @rdname getGroupsFrom
#' @param groups Matrix: groups of dataset
setGroupsFrom <- function(dataset, groups, category=getCategory())
    setGlobal(category, dataset, "groups", value=groups)

#' Get or set differential analyses' elements for a data category
#' @inherit getGlobal
#' @param differential Data frame or matrix: differential analyses table
getDifferentialAnalyses <- function(category=getCategory())
    getGlobal(category, "differentialAnalyses")

#' @rdname getDifferentialAnalyses
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
getDifferentialAnalysesHighlightedEvents <- function(category=getCategory())
    getGlobal(category, "differentialAnalysesHighlighted")

#' @rdname getDifferentialAnalyses
#' @param events Integer: index of events
setDifferentialAnalysesHighlightedEvents <- function(events, 
                                                     category=getCategory())
    setGlobal(category, "differentialAnalysesHighlighted", value=events)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesZoom <- function(category=getCategory())
    getGlobal(category, "differentialAnalysesZoom")

#' @rdname getDifferentialAnalyses
#' @param zoom Integer: range of X and Y coordinates for zooming
setDifferentialAnalysesZoom <- function(zoom, category=getCategory())
    setGlobal(category, "differentialAnalysesZoom", value=zoom)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesSelected <- function(category=getCategory())
    getGlobal(category, "differentialAnalysesSelected")

#' @rdname getDifferentialAnalyses
setDifferentialAnalysesSelected <- function(events, category=getCategory())
    setGlobal(category, "differentialAnalysesSelected", value=events)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesSurvival <- function(category=getCategory())
    getGlobal(category, "diffAnalysesSurv")

#' @rdname getDifferentialAnalyses
#' @param survival Data frame or matrix: differential analyses' survival data
setDifferentialAnalysesSurvival <- function(survival, category=getCategory())
    setGlobal(category, "diffAnalysesSurv", value=survival)

#' @rdname getDifferentialAnalyses
getDifferentialAnalysesColumns <- function(category=getCategory())
    getGlobal(category, "diffAnalysesCols")

#' @rdname getDifferentialAnalyses
#' @param columns Character: differential analyses' column names
setDifferentialAnalysesColumns <- function(columns, category=getCategory())
    setGlobal(category, "diffAnalysesCols", value=columns)