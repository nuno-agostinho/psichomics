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

writeLengthWithinRect <- function(textLen, x1, x2, y=20) {
    drawText(textLen, x=mean(c(x1, x2)), withinRect=TRUE, y=y, anchor="middle",
             baseline="middle")
}

drawRect <- function(x1, x2, text1=NULL, text2=NULL, textLen=NULL, y=10,
                     height=20, fill="#faa000", stroke="#ad7001",
                     showText=TRUE) {
    style <- paste("stroke-width: 1.5px;", "fill: %s", "stroke: %s",
                   sep="; ")
    if (!showText) text1 <- text2 <- textLen <- NULL
    text1   <- if (!is.null(text1))   drawText(text1, x1)
    text2   <- if (!is.null(text2))   drawText(text2, x2, anchor="end")
    textLen <- if (!is.null(textLen)) writeLengthWithinRect(textLen, x1, x2)
    tagList(
        tag("rect", c(class="diagram", x=x1, y=y, width=x2 - x1, height=height,
                      style=sprintf(style, fill, stroke))),
        text1, text2, textLen)
}

# Prepare connecting lines between exons
drawCurve <- function(x1, x2, y=30, ...) {
    drawPath(x1, x2, y=y, type="curve", ...)
}

prepareSVG <- function(exon, path, class=NULL, style=NULL, width=NULL) {
    svg <- tag("svg",
               c(height="50px", class=class, style=style, width=width))
    svg <- tagAppendChildren(svg, exon, path)
    svg <- as.character(svg)
    return(svg)
}

#' Prepare SVG diagram of alternative splicing events
#'
#' @param parsed Alternative splicing event
#' @param type Character: alternative splicing event type
#' @param class Character: class of SVG parent tag
#' @param style Character: style of SVG parent tag
#' @param showText Boolean: display coordinates and length (if available)
#' @param showPath Boolean: display alternative splicing junctions
#' @param showAlternative1 Boolean: show alternative exon 1 and respective
#'   splicing junctions and text?
#' @param showAlternative2 Boolean: show alternative exon 2 and respective
#'   splicing junctions and text? (only related with mutually exclusive exons)
#' @param constitutiveWidth Numeric: width of constitutive exon(s)
#' @param alternativeWidth Numeric: width of alternative exon(s)
#' @param intronWidth Numeric: width of intron's representation
#' @param constitutiveFill Character: fill colour of constitutive exons
#' @param constitutiveStroke Character: stroke colour of constitutive exons
#' @param alternative1Fill Character: fill colour of alternative exon 1
#' @param alternative1Stroke Character: stroke colour of alternative exon 1
#' @param alternative2Fill Character: fill colour of alternative exon 2
#' @param alternative2Stroke Character: stroke colour of alternative exon 2
#'
#' @importFrom shiny tag tagAppendChildren tagList
#'
#' @return Diagrams per alternative splicing event in SVG
#' @keywords internal
diagramSplicingEvent <- function(
    parsed, type, class="pull-right", style=NULL,
    showText=TRUE, showPath=TRUE, showAlternative1=TRUE, showAlternative2=TRUE,
    constitutiveWidth=NULL, alternativeWidth=NULL, intronWidth=NULL,
    constitutiveFill="lightgray", constitutiveStroke="darkgray",
    alternative1Fill="#ffb153", alternative1Stroke="#faa000",
    alternative2Fill="#caa06c", alternative2Stroke="#9d7039") {

    isVASTTOOLS <- isTRUE(any(parsed$source == "vast-tools"))

    isSE   <- type == "SE"
    isMXE  <- type == "MXE"
    isAFE  <- type == "AFE"
    isALE  <- type == "ALE"
    isA3SS <- type == "A3SS"
    isA5SS <- type == "A5SS"
    isRI   <- type == "RI"
    if (!isSE && !isMXE && !isAFE && !isALE && !isA3SS && !isA5SS && !isRI) {
        stop("Unsupported alternative splicing event type")
    } else if (isVASTTOOLS && (isA3SS || isA5SS)) {
        msg <- paste("Alt3 and Alt5 alternative splicing events from",
                     "VAST-TOOLS are currentyly not supported")
        stop(msg)
    }

    # Prepare element width
    if (is.null(alternativeWidth)) {
        alternativeWidth <- 60
        if (isSE || isMXE) alternativeWidth <- 110
    }
    if (is.null(intronWidth)) {
        intronWidth <- 15
        if (isRI && showAlternative1) intronWidth <- 110
    }
    if (is.null(constitutiveWidth)) {
        constitutiveWidth <- 60
        if (isRI) constitutiveWidth <- 110
    }

    pos_C1s <- pos_C1e <- pos_A1s <- pos_A1e <- pos_A2s <- pos_A2e <-
        pos_C2s <- pos_C2e <- A1len <- A2len <- NULL

    safeAreaWidth <- 1
    exon <- tagList()
    path <- tagList()
    if (isSE) {
        pos_C1e  <- as.numeric(parsed$constitutive1)
        pos_A1s  <- as.numeric(sapply(parsed$alternative1, "[[", 1))
        pos_A1e  <- as.numeric(sapply(parsed$alternative1, "[[", 2))
        pos_C2s  <- as.numeric(parsed$constitutive2)
        A1len    <- paste(abs(pos_A1e - pos_A1s), "nts")

        C1s <- safeAreaWidth
        C1e <- C1s + constitutiveWidth
        A1s <- C1e + if (showAlternative1) intronWidth else 0
        A1e <- A1s + if (showAlternative1) alternativeWidth else 0
        C2s <- A1e + intronWidth
        C2e <- C2s + constitutiveWidth
        diagramWidth <- C2e + safeAreaWidth

        exon$C1 <- drawRect(C1s, C1e, text2="%1$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        exon$A1 <- drawRect(A1s, A1e, text1="%2$s", text2="%3$s",
                            textLen="%4$s", showText=showText,
                            fill=alternative1Fill, stroke=alternative1Stroke)
        exon$C2 <- drawRect(C2s, C2e, text1="%5$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        path$C1e_C2s <- drawCurve(C1e, C2s, stroke=constitutiveStroke)
        path$C1e_A1s <- drawCurve(C1e, A1s, stroke=alternative1Stroke)
        path$A1e_C2s <- drawCurve(A1e, C2s, stroke=alternative1Stroke)
        if (!showPath) path <- NULL

        if (!showAlternative1) {
            exon$A1 <- NULL
            path[grep("A1", names(path), fixed=TRUE)] <- NULL
        }

        # Replace with genomic positions per event
        svg      <- prepareSVG(exon, path, class, style, diagramWidth)
        svgFinal <- sprintf(svg, pos_C1e, pos_A1s, pos_A1e, A1len, pos_C2s)
    } else if (isMXE) {
        pos_C1e  <- as.numeric(parsed$constitutive1)
        pos_A1s  <- as.numeric(sapply(parsed$alternative1,  "[[", 1))
        pos_A1e  <- as.numeric(sapply(parsed$alternative1,  "[[", 2))
        pos_A2s  <- as.numeric(sapply(parsed$alternative2,  "[[", 1))
        pos_A2e  <- as.numeric(sapply(parsed$alternative2,  "[[", 2))
        pos_C2s  <- as.numeric(parsed$constitutive2)

        A1len    <- paste(abs(pos_A1e - pos_A1s), "nts")
        A2len    <- paste(abs(pos_A2e - pos_A2s), "nts")

        isA1ref  <- !xor(parsed$strand == "+", pos_A1e <= pos_A2s)
        hideA1   <- (!showAlternative1 && isA1ref) ||
            (!showAlternative2 && !isA1ref)
        hideA2   <- (!showAlternative2 && isA1ref) ||
            (!showAlternative1 && !isA1ref)

        C1s <- safeAreaWidth
        C1e <- C1s + constitutiveWidth
        A1s <- C1e + if (!hideA1) intronWidth else 0
        A1e <- A1s + if (!hideA1) alternativeWidth else 0
        A2s <- A1e + if (!hideA2) intronWidth else 0
        A2e <- A2s + if (!hideA2) alternativeWidth else 0
        C2s <- A2e + intronWidth
        C2e <- C2s + constitutiveWidth
        diagramWidth <- C2e + safeAreaWidth

        exon$C1 <- drawRect(C1s, C1e, text2="%1$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        exon$A1 <- drawRect(A1s, A1e, text1="%2$s", text2="%3$s",
                            textLen="%4$s", showText=showText,
                            fill="%9$s", stroke="%10$s")
        exon$A2 <- drawRect(A2s, A2e, text1="%5$s", text2="%6$s",
                            textLen="%7$s", showText=showText,
                            fill="%11$s", stroke="%12$s")
        exon$C2 <- drawRect(C2s, C2e, text1="%8$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        path$C1e_A1s <- drawCurve(C1e, A1s, stroke="%10$s")
        path$A1e_C2s <- drawCurve(A1e, C2s, stroke="%10$s")
        path$C1e_A2s <- drawCurve(C1e, A2s, stroke="%12$s")
        path$A2e_C2s <- drawCurve(A2e, C2s, stroke="%12$s")
        if (!showPath) path <- NULL

        if (hideA1) {
            exon$A1 <- NULL
            path[grep("A1", names(path), fixed=TRUE)] <- NULL
        }
        if (hideA2) {
            exon$A2 <- NULL
            path[grep("A2", names(path), fixed=TRUE)] <- NULL
        }

        # Replace with genomic positions per event and correct reference exon
        svg      <- prepareSVG(exon, path, class, style, diagramWidth)
        svgFinal <- sprintf(
            svg, pos_C1e,
            ifelse(isA1ref, pos_A1s, pos_A2s),
            ifelse(isA1ref, pos_A1e, pos_A2e),
            ifelse(isA1ref, A1len, A2len),
            ifelse(isA1ref, pos_A2s, pos_A1s),
            ifelse(isA1ref, pos_A2e, pos_A1e),
            ifelse(isA1ref, A2len, A1len),
            pos_C2s,
            ifelse(isA1ref, alternative1Fill, alternative2Fill),
            ifelse(isA1ref, alternative1Stroke, alternative2Stroke),
            ifelse(isA1ref, alternative2Fill, alternative1Fill),
            ifelse(isA1ref, alternative2Stroke, alternative1Stroke))
    } else if (isAFE || isA5SS) {
        pos_A1e  <- as.numeric(parsed$alternative1)
        pos_A2e  <- as.numeric(parsed$alternative2)
        pos_C2s  <- as.numeric(parsed$constitutive2)

        isA2ref  <- !xor(parsed$strand == "+", pos_A1e >= pos_A2e)
        hideA1   <- (!showAlternative1 && isA2ref) ||
            (!showAlternative2 && !isA2ref)
        hideA2   <- (!showAlternative2 && isA2ref) ||
            (!showAlternative1 && !isA2ref)

        A1s <- safeAreaWidth
        A1e <- A1s + if (!hideA1) alternativeWidth else 0
        A2s <- A1e + if (!hideA1 && isAFE) intronWidth else 0
        A2e <- A2s + if (!hideA2) alternativeWidth else 0
        C2s <- A2e + if (!hideA1 || !hideA2) intronWidth else 0
        C2e <- C2s + constitutiveWidth
        diagramWidth <- C2e + safeAreaWidth

        exon$A1 <- drawRect(A1s, A1e, text2="%1$s", showText=showText,
                            fill="%4$s", stroke="%5$s")
        exon$A2 <- drawRect(A2s, A2e, text2="%2$s", showText=showText,
                            fill="%6$s", stroke="%7$s")
        exon$C2 <- drawRect(C2s, C2e, text1="%3$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        path$A1e_C2s <- drawCurve(A1e, C2s, stroke="%5$s")
        path$A2e_C2s <- drawCurve(A2e, C2s, stroke="%7$s")
        if (!showPath) path <- NULL

        if (hideA1) {
            exon$A1 <- NULL
            path[grep("A1", names(path), fixed=TRUE)] <- NULL
        }
        if (hideA2) {
            exon$A2 <- NULL
            path[grep("A2", names(path), fixed=TRUE)] <- NULL
        }

        # Replace with genomic positions per event and correct reference exon
        svg      <- prepareSVG(exon, path, class, style, diagramWidth)
        svgFinal <- sprintf(
            svg,
            ifelse(isA2ref, pos_A2e, pos_A1e),
            ifelse(isA2ref, pos_A1e, pos_A2e),
            pos_C2s,
            ifelse(isA2ref, alternative1Fill, constitutiveFill),
            ifelse(isA2ref, alternative1Stroke, constitutiveStroke),
            ifelse(isA2ref, constitutiveFill, alternative1Fill),
            ifelse(isA2ref, constitutiveStroke, alternative1Stroke))
    } else if (isALE || isA3SS) {
        pos_C1s  <- as.numeric(parsed$constitutive1)
        pos_A1e  <- as.numeric(parsed$alternative1)
        pos_A2e  <- as.numeric(parsed$alternative2)

        isA1ref  <- !xor(parsed$strand == "+", pos_A1e <= pos_A2e)

        hideA1 <- (!showAlternative1 && isA1ref) ||
            (!showAlternative2 && !isA1ref)
        hideA2 <- (!showAlternative2 && isA1ref) ||
            (!showAlternative1 && !isA1ref)

        C1s <- safeAreaWidth
        C1e <- C1s + constitutiveWidth
        A1s <- C1e + if (!hideA1) intronWidth else 0
        A1e <- A1s + if (!hideA1) alternativeWidth else 0
        A2s <- A1e + if (!hideA2 && isALE) intronWidth else 0
        A2e <- A2s + if (!hideA2) alternativeWidth else 0
        diagramWidth <- A2e + safeAreaWidth

        exon$C1 <- drawRect(C1s, C1e, text2="%1$s", showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        exon$A1 <- drawRect(A1s, A1e, text1="%2$s", showText=showText,
                            fill="%4$s", stroke="%5$s")
        exon$A2 <- drawRect(A2s, A2e, text1="%3$s", showText=showText,
                            fill="%6$s", stroke="%7$s")
        path$C1e_A1s <- drawCurve(C1e, A1s, stroke="%5$s")
        path$C1e_A2s <- drawCurve(C1e, A2s, stroke="%7$s")
        if (!showPath) path <- NULL

        if (hideA1) {
            exon$A1 <- NULL
            path[grep("A1", names(path), fixed=TRUE)] <- NULL
        }
        if (hideA2) {
            exon$A2 <- NULL
            path[grep("A2", names(path), fixed=TRUE)] <- NULL
        }

        # Replace with genomic positions per event and correct reference exon
        svg      <- prepareSVG(exon, path, class, style, diagramWidth)
        svgFinal <- sprintf(
            svg,
            pos_C1s,
            ifelse(isA1ref, pos_A1e, pos_A2e),
            ifelse(isA1ref, pos_A2e, pos_A1e),
            ifelse(isA1ref, alternative1Fill, constitutiveFill),
            ifelse(isA1ref, alternative1Stroke, constitutiveStroke),
            ifelse(isA1ref, constitutiveFill, alternative1Fill),
            ifelse(isA1ref, constitutiveStroke, alternative1Stroke))
    } else if (isRI) {
        pos_C1s  <- as.numeric(sapply(parsed$constitutive1, "[[", 1))
        pos_C1e  <- as.numeric(sapply(parsed$constitutive1, "[[", 2))
        pos_C2s  <- as.numeric(sapply(parsed$constitutive2, "[[", 1))
        pos_C2e  <- as.numeric(sapply(parsed$constitutive2, "[[", 2))
        A1len    <- paste(abs(pos_C1e - pos_C2s), "nts")

        C1s <- safeAreaWidth
        C1e <- C1s + constitutiveWidth
        A1s <- C1e
        A1e <- A1s + intronWidth
        C2s <- A1e
        C2e <- C2s + constitutiveWidth
        diagramWidth <- C2e + safeAreaWidth

        exon$C1 <- drawRect(C1s, C1e, text1="%1$s", text2="%2$s",
                            showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)
        exon$A1 <- drawRect(A1s, A1e, textLen="%3$s", showText=showText,
                            fill=alternative1Fill, stroke=alternative1Stroke,
                            height=12, y=14)
        exon$C2 <- drawRect(C2s, C2e, text1="%4$s", text2="%5$s",
                            showText=showText,
                            fill=constitutiveFill, stroke=constitutiveStroke)

        path$C1e_C2s <- drawCurve(C1e, C2s, stroke=constitutiveStroke)
        if (!showPath) path <- NULL

        if (!showAlternative1) exon$A1 <- NULL

        # Replace with genomic positions per event
        svg      <- prepareSVG(exon, path, class, style, diagramWidth)
        svgFinal <- sprintf(svg, pos_C1s, pos_C1e, A1len, pos_C2s, pos_C2e)
    }
    names(svgFinal) <- rownames(parsed)

    addClass <- function(char) {
        # char        <- HTML(char)
        class(char) <- c("splicingEventPlot", class(char))
        return(char)
    }
    svgFinal <- lapply(svgFinal, addClass)

    attr(svgFinal, "positions") <- paste(pos_C1s, pos_C1e, pos_A1s, pos_A1e,
                                         pos_A2s, pos_A2e, pos_C2s, pos_C2e,
                                         A1len, A2len)
    return(svgFinal)
}

#' Plot diagram of alternative splicing events
#'
#' @param ASevent Character: alternative splicing event identifiers
#' @inheritParams parseSplicingEvent
#' @inheritParams diagramSplicingEvent
#'
#' @importFrom shiny tag tagList
#'
#' @return List of SVG (one for each alternative splicing event)
#' @export
#'
#' @examples
#' events <- c(
#'   "A3SS_15_+_63353138_63353912_63353397_TPM1",
#'   "A3SS_11_-_61118463_61117115_61117894_CYB561A3",
#'   "A5SS_21_+_48055675_48056459_48056808_PRMT2",
#'   "A5SS_1_-_1274742_1274667_1274033_DVL1",
#'   "AFE_9_+_131902430_131901928_131904724_PPP2R4",
#'   "AFE_5_-_134686513_134688636_134681747_H2AFY",
#'   "ALE_12_+_56554104_56554410_56555171_MYL6",
#'   "ALE_8_-_38314874_38287466_38285953_FGFR1",
#'   "SE_9_+_6486925_6492303_6492401_6493826_UHRF2",
#'   "SE_19_-_5218431_5216778_5216731_5215606_PTPRS",
#'   "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1",
#'   "MXE_17_-_74090495_74087316_74087224_74086478_74086410_74085401_EXOC7")
#' diagram <- plotSplicingEvent(events)
#'
#' \dontrun{
#' diagram[["A3SS_3_-_145796903_145794682_145795711_PLOD2"]]
#' diagram[[6]]
#' diagram
#' }
plotSplicingEvent <- function(
    ASevent, data=NULL, showText=TRUE, showPath=TRUE,
    showAlternative1=TRUE, showAlternative2=TRUE,
    constitutiveWidth=NULL, alternativeWidth=NULL, intronWidth=NULL,
    constitutiveFill="lightgray", constitutiveStroke="darkgray",
    alternative1Fill="#ffb153", alternative1Stroke="#faa000",
    alternative2Fill="#caa06c", alternative2Stroke="#9d7039",
    class=NULL, style=NULL) {

    if (is.data.frame(ASevent)) {
        if (is.null(data)) data <- getSplicingEventData(ASevent)
        ASevent <- rownames(ASevent)
    }

    parsed <- parseSplicingEvent(ASevent, coords=TRUE, data=data)
    if (!is.null(parsed$type)) {
        parsed <- split(parsed, parsed$type)
    } else {
        type   <- rep("NULL", length(ASevent))
        parsed <- split(ASevent, type)
    }

    svg <- NULL
    positions <- NULL
    showWarning <- FALSE
    for (type in names(parsed)) {
        parsedType <- parsed[[type]]
        diagram <- tryCatch(
            diagramSplicingEvent(
                parsedType, type, class=class, style=style,
                showText=showText, showPath=showPath,
                showAlternative1=showAlternative1,
                showAlternative2=showAlternative2,
                constitutiveWidth=constitutiveWidth,
                alternativeWidth=alternativeWidth, intronWidth=intronWidth,
                constitutiveFill=constitutiveFill,
                constitutiveStroke=constitutiveStroke,
                alternative1Fill=alternative1Fill,
                alternative1Stroke=alternative1Stroke,
                alternative2Fill=alternative2Fill,
                alternative2Stroke=alternative2Stroke),
            error=return)

        if (is(diagram, "error")) {
            showWarning                <- TRUE

            len <- nrow(parsedType)
            if (is.null(len)) len <- length(parsedType)

            ns  <- rownames(parsedType)
            if (is.null(ns)) ns <- parsedType

            diagram                    <- rep("", len)
            attr(diagram, "positions") <- diagram
            names(diagram)             <- ns
        }
        svg <- c(svg, diagram)
        positions <- c(positions, attr(diagram, "positions"))
    }
    # Order results based on original input
    ord <- match(names(svg), ASevent)
    svg <- svg[ord]
    class(svg) <- c("splicingEventPlotList", class(svg))
    attr(svg, "positions") <- positions[ord]

    if (showWarning) {
        msg <- paste("Diagrams could not be rendered for one or more",
                     "alternative splicing event types.")
        warning(msg)
    }
    return(svg)
}

#' @importFrom htmltools knit_print.html
#' @export
print.splicingEventPlot <- function(x, ..., browse=interactive()) {
    x <- HTML(x)
    if (!browse) {
        return(knit_print.html(x))
    } else {
        print(x, ..., browse=browse)
        return(invisible(x))
    }
}

#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom shiny fluidPage
#' @export
print.splicingEventPlotList <- function(x, ..., browse=interactive()) {
    server <- function(input, output) {
        prepareData <- reactive(
            data.frame(cbind(names(x), x), stringsAsFactors=FALSE))
        output$eventsTable <- renderDataTable(
            prepareData(), rownames=FALSE, class="compact display",
            colnames=c("Alternative splicing event", "Diagram"),
            style="bootstrap",
            caption="Diagrams of alternative splicing events")
    }

    if (length(x) == 1) {
        print(x[[1]], ..., browse=browse)
    } else {
        runApp(list(ui=fluidPage(dataTableOutput("eventsTable")),
                    server=server))
    }
}
