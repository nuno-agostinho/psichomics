# Test if event diagrams colour the same reference exon
plotEventTranscripts <- function(e, species="human", assembly="hg19") {
    info <- queryEnsemblByEvent(e, species=species, assembly=assembly)
    plotTranscripts(info, event=e)
}

getReferenceExonFromPlotEventTranscripts <- function(events, ...,
                                                     colour="#ffb153") {
    getReferenceExonForOneEvent <- function(e, ...) {
        hc    <- plotASeventRegion(highchart(), e)
        bands <- hc$x$hc_opts$xAxis$plotBands

        # Find coordinates of orange exon
        index <- grep(colour, bands, fixed=TRUE)
        found <- bands[[index]]
        return(c(found$from, found$to))
    }
    pos <- pbapply::pblapply(events, getReferenceExonForOneEvent, ...)
    names(pos) <- events
    return(pos)
}

getReferenceExonFromPlotSplicingEvent <- function(..., colour="#ffb153") {
    p <- plotSplicingEvent(...)
    pattern <- sprintf(".*<rect.*%s.*>([0-9]*?)</text>.*", colour)
    pos <- pbapply::pblapply(p, function(i) {
        as.numeric(gsub(i[[1]], pattern=pattern, replacement="\\1"))
    })
    return(pos)
}

getReferenceExonFromEventID <- function(...) {
    parsed <- parseSplicingEvent(..., coords=TRUE)
    setNames(parsed$alternative1, rownames(parsed))
}

compareReferenceExons <- function(events, eventData=NULL) {
    message("Extracting coordinates of exonic references from...")
    message("  --> event IDs")
    ref         <- getReferenceExonFromEventID(events, data=eventData)
    message("  --> splicing diagrams")
    diagrams    <- getReferenceExonFromPlotSplicingEvent(events, data=eventData)
    message("  --> transcript plots...")
    transcripts <- getReferenceExonFromPlotEventTranscripts(events)

    comp <- function(event, ref, a, b) {
        any(a[[event]] %in% ref[[event]]) && any(b[[event]] %in% ref[[event]])
    }
    res <- pbapply::pbsapply(names(ref), comp, ref, transcripts, diagrams)
    names(res) <- names(ref)
    return(res)
}

compareReferenceExonsInData <- function(data, assembly="hg19") {
    junctionQuant <- data[[1]]$`Junction quantification`
    annot <- loadAnnotation(listSplicingAnnotations(assembly=assembly)[[1]])
    psi <- quantifySplicing(annot, junctionQuant)
    events <- rownames(psi)
    cmp <- compareReferenceExons(events, attr(psi, "rowData"))
    return(cmp)
}

res  <- list()
res$TCGA$data <- loadTCGAdata(cohort="BRCA")
res$TCGA$cmp  <- compareReferenceExonsInData(res$TCGA$data, "hg19")
print(length(res$TCGA$cmp))
print(all(res$TCGA$cmp))

res$SRA$data <- loadSRAproject("SRP053101")
res$SRA$cmp  <- compareReferenceExonsInData(res$SRA$data, "hg38")
print(length(res$SRA$cmp))
print(all(res$SRA$cmp))

options(timeout=1000)
res$GTEx$data <- loadGtexData(tissue="Adipose Tissue")
res$GTEx$cmp  <- compareReferenceExonsInData(res$GTEx$data, "hg38")
print(length(res$GTEx$cmp))
print(all(res$GTEx$cmp))
