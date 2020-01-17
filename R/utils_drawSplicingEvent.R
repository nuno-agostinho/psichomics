drawRect <- function(x1, x2, y=10, height=20, fill="#faa000", 
                     stroke="#ad7001") {
    style <- paste("stroke-width: 1.5px;", "fill: %s", "stroke: %s", 
                   sep="; ")
    tag("rect", c(class="diagram", x=x1, y=y, width=x2 - x1, height=height,
                  style=sprintf(style, fill, stroke)))
}

drawPath <- function(x1, x2, y=20, type=c("line", "curve"), stroke="black") {
    type  <- match.arg(type)
    style <- paste("fill: none", paste("stroke:", stroke),
                   "stroke-width: 1.5px", sep="; ")
    if (type == "line") {
        tag("path", c(class="diagram", style=style,
                      d=sprintf("M %s %s L %s %s", x1, y, x2, y)))
    } else if (type == "curve") {
        curve <- 45
        curveAnchor <- round((x2 - x1) / 3)
        tag("path", c(class="diagram", style=style, d=sprintf(
            "M %s %s C %s %s %s %s %s %s",
            x1, y, x1 + curveAnchor, curve, x2 - curveAnchor, curve, x2, y)))
    }
}

drawText <- function(text, x, y=8, anchor="start", baseline="auto", 
                     withinRect=FALSE) {
    if (withinRect)
        class <- "diagram"
    else
        class <- "diagram outside"
    
    style <- paste("font-size: 10px", "font-family: Helvetica, sans-serif",
                   "text-anchor: %s", "dominant-baseline: %s", sep="; ")
    tag("text", c(class=class, x=x, y=y, style=sprintf(style, anchor, baseline),
                  text))
}

#' Prepare SVG diagram of alternative splicing events
#' 
#' @param parsed Alternative splicing event
#' @param type Character: alternative splicing event type
#' @param class Character: class of SVG parent tag
#' @param style Character: style of SVG parent tag
#' @param showText Boolean: display coordinates and exon length (if available)
#' @param showPath Boolean: display alternative splicing junctions
#' @param showAlternative1 Boolean: show alternative exon 1 and respective
#'   splicing junctions and text?
#' @param showAlternative2 Boolean: show alternative exon 2 and respective 
#'   splicing junctions and text? (only related with mutually exclusive exons)
#' @param constitutiveLength Numeric: length of constitutive exon(s)'s 
#'   representation
#' @param alternativeLength Numeric: length of alternative exon(s)'s
#'   representation
#' @param intronLength Numeric: length of intron's representation
#' @param constitutiveFill Character: fill colour of constitutive exons
#' @param constitutiveStroke Character: stroke colour of constitutive exons
#' @param alternative1Fill Character: fill colour of alternative exon 1
#' @param alternative1Stroke Character: stroke colour of alternative exon 1
#' @param alternative2Fill Character: fill colour of alternative exon 2 (only
#'   required for mutually exclusive exons)
#' @param alternative2Stroke Character: stroke colour of alternative exon 2
#'   (only required for mutually exclusive exons)
#'
#' @importFrom shiny tag tagAppendChildren tagList
#' 
#' @return Diagrams per alternative splicing event in SVG
#' @keywords internal
diagramSplicingEvent <- function(
    parsed, type, class="pull-right", style=NULL, showText=TRUE, showPath=TRUE,
    showAlternative1=TRUE, showAlternative2=TRUE,
    constitutiveLength=60, alternativeLength=NULL, intronLength=15, 
    constitutiveFill="lightgray", constitutiveStroke="darkgray", 
    alternative1Fill="#ffb153", alternative1Stroke="#faa000", 
    alternative2Fill="#caa06c", alternative2Stroke="#9d7039") {
    
    safeAreaLength <- 1
    
    isSE   <- type == "SE"
    isMXE  <- type == "MXE"
    isAFE  <- type == "AFE"
    isALE  <- type == "ALE"
    isA3SS <- type == "A3SS"
    isA5SS <- type == "A5SS"
    if (!isSE && !isMXE && !isAFE && !isALE && !isA3SS && !isA5SS)
        stop("Unsupported alternative splicing event type")
    
    if (is.null(alternativeLength)) {
        alternativeLength <- 60
        if (isSE || isMXE) alternativeLength <- 110
    }
    
    # Exon coordinates in the diagram
    C1s <- safeAreaLength
    C1e <- C1s + constitutiveLength
    if (isA5SS) {
        A1s <- C1e
    } else {
        A1s <- C1e + intronLength
    }
    A1e <- A1s + alternativeLength
    
    if (!showAlternative1) A1e <- A1s
    
    if (isA3SS) {
        C2s <- A1e
    } else if (isMXE) {
        A2s <- A1e + intronLength
        A2e <- A2s + alternativeLength
        if (!showAlternative2) A2e <- A2s
        C2s <- A2e + intronLength
    } else {
        C2s <- A1e + intronLength
    }
    C2e <- C2s + constitutiveLength
    diagramWidth <- C2e + safeAreaLength
    
    # Prepare exons
    exon <- tagList()
    exon$C1 <- drawRect(C1s, C1e, fill=constitutiveFill,
                        stroke=constitutiveStroke)
    if (showAlternative1) {
        exon$A1 <- drawRect(A1s, A1e, fill=alternative1Fill, 
                            stroke=alternative1Stroke)
    }
    if (isMXE && showAlternative2) {
        exon$A2 <- drawRect(A2s, A2e, fill=alternative2Fill,
                            stroke=alternative2Stroke)
    }
    exon$C2 <- drawRect(C2s, C2e, fill=constitutiveFill, 
                        stroke=constitutiveStroke)
    
    # Prepare connecting lines between exons
    drawCurve <- function(x1, x2, y=30, ...)
        drawPath(x1, x2, y=y, type="curve", ...)
    
    if (showPath) {
        path <- tagList()
        if (!isMXE) {
            path$C1e_C2s <- drawCurve(C1e, C2s, stroke=constitutiveStroke)
        }
        if ((isMXE || isAFE) && showAlternative1) {
            path$A1e_C2s <- drawCurve(A1e, C2s, stroke=alternative1Stroke)
        }
        if ((isMXE || isALE) && showAlternative1) {
            path$C1e_A1s <- drawCurve(C1e, A1s, stroke=alternative1Stroke)
        }
        if ((isSE || isA3SS) && showAlternative1) {
            path$C1e_A1s <- drawPath(C1e, A1s, stroke=alternative1Stroke)
        }
        if ((isSE || isA5SS) && showAlternative1) {
            path$A1e_C2s <- drawPath(A1e, C2s, stroke=alternative1Stroke)
        }
        if (isMXE && showAlternative2) {
            path$C1e_A2s <- drawCurve(C1e, A2s, stroke=alternative2Stroke)
            path$A2e_C2s <- drawCurve(A2e, C2s, stroke=alternative2Stroke)
        }
    } else {
        path <- NULL
    }
    
    # Prepare text for genome coordinates and exon length
    if (showText) {
        text     <- tagList()
        text$C1e <- drawText("%s", C1e, anchor="end")
        
        if (showAlternative1) {
            if (isSE || isMXE || isALE || isA3SS) {
                text$A1s <- drawText("%s", A1s)
            }
            if (isSE || isMXE || isAFE || isA5SS) {
                text$A1e <- drawText("%s", A1e, anchor="end")
            }
        }
        
        writeLengthWithinRect <- function(x1, x2, y=20) {
            drawText("%s nts", x=mean(c(x1, x2)), withinRect=TRUE,
                     y=y, anchor="middle", baseline="middle")
        }
        if ((isSE || isMXE || isA5SS || isA3SS) && showAlternative1) {
            text$A1len <- writeLengthWithinRect(A1s, A1e)
        }
        if (isMXE && showAlternative2) {
            text$A2s <- drawText("%s", A2s)
            text$A2e <- drawText("%s", A2e, anchor="end")
            text$A2len <- writeLengthWithinRect(A2s, A2e)
        }
        text$C2s <- drawText("%s", C2s)
    } else {
        text <- NULL
    }
    
    # Finalise SVG
    svg <- tag("svg", 
               c(height="50px", class=class, style=style, width=diagramWidth))
    svg <- tagAppendChildren(svg, exon, path, text)
    
    if (isSE || isA5SS || isA3SS) {
        A1len <- sapply(parsed, function(i)
            as.numeric(i[[5]]) - as.numeric(i[[6]]))
        A1len <- abs(A1len)
        
        if (isSE) {
            svgFinal <- sprintf(
                as.character(svg),
                sapply(parsed, "[[", 4), 
                sapply(parsed, "[[", 5), sapply(parsed, "[[", 6), A1len,
                sapply(parsed, "[[", 7))
        } else if (isA5SS || isA3SS) {
            svgFinal <- sprintf(
                as.character(svg), sapply(parsed, "[[", 4),
                sapply(parsed, "[[", 5), A1len, sapply(parsed, "[[", 6))
        }
    } else if (isMXE) {
        A1len <- sapply(parsed, function(i)
            as.numeric(i[[5]]) - as.numeric(i[[6]]))
        A1len <- abs(A1len)
        A2len <- sapply(parsed, function(i)
            as.numeric(i[[7]]) - as.numeric(i[[8]]))
        A2len <- abs(A2len)
        
        svgFinal <- sprintf(
            as.character(svg), sapply(parsed, "[[", 4),
            sapply(parsed, "[[", 5), sapply(parsed, "[[", 6), A1len,
            sapply(parsed, "[[", 7), sapply(parsed, "[[", 8), A2len,
            sapply(parsed, "[[", 9))
    } else {
        svgFinal <- sprintf(as.character(svg), sapply(parsed, "[[", 4),
                            sapply(parsed, "[[", 5), sapply(parsed, "[[", 6))
    }
    names(svgFinal) <- names(parsed)
    return(svgFinal)
}

#' Plot diagram of alternative splicing events
#' 
#' @param ASevent Character: alternative splicing event identifiers
#' @inheritParams diagramSplicingEvent
#' @param raw Boolean: if \code{FALSE}, plot the events; if \code{TRUE}, return
#'   the respective SVG raw code
#' 
#' @importFrom shiny tag tagList
#' 
#' @return List of SVG (one for each alternative splicing event)
#' @export
#' 
#' @examples
#' ASevent <- c(
#'     "SE_9_+_6486925_6492303_6492401_6493826_UHRF2",
#'     "SE_11_+_86925_92303_92401_93826_TESTING/MHG2",
#'     "A5SS_15_+_63353472_63353987_63354414_TPM1",
#'     "A3SS_3_-_145796903_145794682_145795711_PLOD2",
#'     "AFE_17_-_15165746_15168471_15164078_PMP22",
#'     "ALE_18_-_5395066_5394792_5393477_EPB41L3",
#'     "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1")
#' diagram <- plotSplicingEvent(ASevent)
#' 
#' diagram
#' diagram[[6]]
#' diagram[["A3SS_3_-_145796903_145794682_145795711_PLOD2"]]
plotSplicingEvent <- function(
    ASevent, raw=FALSE, class=NULL, style=NULL, 
    showText=TRUE, showPath=TRUE, showAlternative1=TRUE, showAlternative2=TRUE,
    constitutiveLength=60, alternativeLength=NULL, intronLength=15,
    constitutiveFill="lightgray", constitutiveStroke="darkgray",
    alternative1Fill="#ffb153", alternative1Stroke="#faa000",
    alternative2Fill="#caa06c", alternative2Stroke="#9d7039") {
    
    # Custom, faster parser of alternative splicing event type
    parsed <- strsplit(ASevent, "_")
    names(parsed) <- ASevent
    type <- sapply(parsed, "[[", 1)
    
    parsed <- split(parsed, type)
    svg <- NULL
    for (each in names(parsed)) {
        svg <- c(svg, diagramSplicingEvent(
            parsed[[each]], each, class=class, style=style,
            showPath=showPath, showText=showText,
            showAlternative1=showAlternative1,
            showAlternative2=showAlternative2,
            constitutiveLength=constitutiveLength, 
            alternativeLength=alternativeLength, intronLength=intronLength,
            constitutiveFill=constitutiveFill, 
            constitutiveStroke=constitutiveStroke,
            alternative1Fill=alternative1Fill,
            alternative1Stroke=alternative1Stroke, 
            alternative2Fill=alternative2Fill,
            alternative2Stroke=alternative2Stroke))
    }
    
    # Order results based on original input
    svg <- svg[ASevent]
    
    # Plottable SVG images
    if (!raw) svg <- lapply(lapply(svg, HTML), browsable)
    return(svg)
}