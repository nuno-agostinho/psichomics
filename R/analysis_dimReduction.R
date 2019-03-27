#' Reduce dimensionality after processing missing values from data frame
#' 
#' @param data Data frame: data
#' @param type Character: dimensionality reduction technique (\code{pca} or 
#' \code{ica})
#' @param naTolerance Integer: percentage of tolerated missing values per column
#' (deprecated)
#' @param missingValues Integer: number of tolerated missing values per column
#' to be replaced with the mean of the values of that same column (5% of total
#' rows by default)
#' @param scale. Boolean: scale variables?
#' @param ... Extra parameters passed to FUN
#' @inheritParams base::scale
#' 
#' @importFrom stats prcomp
#' @importFrom fastICA fastICA
#' @importFrom miscTools rowMedians colMedians
#' 
#' @return PCA result in a \code{prcomp} object or ICA result
#' object
#' @keywords internal
reduceDimensionality <- function(data, type=c("pca", "ica"), center=TRUE, 
                                 scale.=FALSE, naTolerance=NULL, 
                                 missingValues=round(0.05 * ncol(data)), ...) {
    # # Get individuals (rows) with less than a given percentage of NAs
    # nas <- rowSums(is.na(data))
    # # hist(nas/ncol(data)*100)
    # data <- data[nas/ncol(data)*100 <= naTolerance, , drop=FALSE]
    # if (nrow(data) == 0) return(NULL)
    
    # # Replace NAs with the medians for each individual (row)
    # medians <- rowMedians(data, na.rm=TRUE)
    # data[is.na(data)] <- rep(medians, sum(is.na(data)))
    
    # Get loadings (columns) with less than a given percentage of NAs
    nas <- colSums(is.na(data))
    
    if (!is.null(naTolerance)) {
        warning("The argument 'naTolerance' is deprecated:",
                "use 'missingValues' instead.")
        data <- data[ , nas/nrow(data) * 100 <= naTolerance, drop=FALSE]
    } else {
        data <- data[ , nas <= missingValues, drop=FALSE]
    }
    
    if (ncol(data) == 0) {
        warning("Empty data input. ",
                "Try increasing the tolerance for missing values.")
        return(NULL)
    }
    
    # Replace NAs with the medians for each loading (column)
    medians <- colMedians(data, na.rm=TRUE)
    nas <- colSums(is.na(data))
    data[is.na(data)] <- rep(medians, nas)
    
    if (type == "pca") {
        # Perform principal component analysis
        res <- tryCatch(prcomp(data, center=center, scale.=scale., ...), 
                        error=return)
    } else if (type == "ica") {
        # Perform independent component analysis
        data <- scale(data, scale=scale., center=center)
        res <- tryCatch(fastICA(data, ...), error=return)
        
        # Rename colnames
        if (!is(res, "error") && !is.null(res$S)) {
            colnames(res$S) <- paste0("IC", seq(ncol(res$S)))
        }
        
        # # Fix rownames for C implementation of fastICA
        # res <- tryCatch(fastICA(data, method="C", ...), error=return)
        # if (!is(res, "error")) {
        #     if (!is.null(res$X)) rownames(res$X) <- rownames(data)
        #     if (!is.null(res$S)) rownames(res$S) <- rownames(data)
        # }
    }
    
    # Result is useless if it only has one point
    if ("x" %in% names(res) && nrow(res$x) == 1) res <- NULL
    return(res)
}

#' Add clusters to \code{highchart} object
#' 
#' Clusters are added as coloured polygons.
#' 
#' @param hc \code{highchart} object
#' @param data Data frame
#' @param clustering Character: group of each sample
#' 
#' @importFrom grDevices chull
#' 
#' @return \code{highcharter} object
#' @keywords internal
plotClusters <- function(hc, data, clustering) {
    for ( each in sort(unique(clustering)) ) {
        df <- data[clustering == each, , drop=FALSE]
        df <- df[chull(df), , drop=FALSE] # cluster points' convex hull
        
        colour <- JS(paste0(
            "Highcharts.Color(Highcharts.getOptions().",
            "colors[", each, "]).setOpacity(0.3).get()"))
        
        if (nrow(df) <= 2) {
            hc <- hc %>% hc_add_series(
                df, zIndex=-1, color=colour,
                name=paste("Cluster", each), lineWidth=8,
                marker=list(radius=8, symbol="circle"))
        } else {
            hc <- hc %>% hc_add_series(
                df, type="polygon", zIndex=-1, color=colour,
                name=paste("Cluster", each))
        }
    }
    return(hc)
}

#' @rdname appUI
#' @importFrom shiny NS
dimReductionUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "dimReduction", 
                             priority=c("pcaUI", "icaUI"))
    return(uiList)
}

#' @rdname appServer
#' 
#' @importFrom shiny observe observeEvent renderPlot
#' @importFrom shinyjs hide show
dimReductionServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("dimReduction",
                                 priority=c("pcaServer", "icaServer"))
}

attr(dimReductionUI, "loader") <- "analysis"
attr(dimReductionUI, "name") <- "Dimensionality reduction techniques"
attr(dimReductionServer, "loader") <- "analysis"