% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_drawSplicingEvent.R
\name{diagramSplicingEvent}
\alias{diagramSplicingEvent}
\title{Prepare SVG diagram of alternative splicing events}
\usage{
diagramSplicingEvent(
  parsed,
  type,
  class = "pull-right",
  style = NULL,
  showText = TRUE,
  showPath = TRUE,
  showAlternative1 = TRUE,
  showAlternative2 = TRUE,
  constitutiveWidth = NULL,
  alternativeWidth = NULL,
  intronWidth = NULL,
  constitutiveFill = "lightgray",
  constitutiveStroke = "darkgray",
  alternative1Fill = "#ffb153",
  alternative1Stroke = "#faa000",
  alternative2Fill = "#caa06c",
  alternative2Stroke = "#9d7039"
)
}
\arguments{
\item{parsed}{Alternative splicing event}

\item{type}{Character: alternative splicing event type}

\item{class}{Character: class of SVG parent tag}

\item{style}{Character: style of SVG parent tag}

\item{showText}{Boolean: display coordinates and length (if available)}

\item{showPath}{Boolean: display alternative splicing junctions}

\item{showAlternative1}{Boolean: show alternative exon 1 and respective
splicing junctions and text?}

\item{showAlternative2}{Boolean: show alternative exon 2 and respective
splicing junctions and text? (only related with mutually exclusive exons)}

\item{constitutiveWidth}{Numeric: width of constitutive exon(s)}

\item{alternativeWidth}{Numeric: width of alternative exon(s)}

\item{intronWidth}{Numeric: width of intron's representation}

\item{constitutiveFill}{Character: fill colour of constitutive exons}

\item{constitutiveStroke}{Character: stroke colour of constitutive exons}

\item{alternative1Fill}{Character: fill colour of alternative exon 1}

\item{alternative1Stroke}{Character: stroke colour of alternative exon 1}

\item{alternative2Fill}{Character: fill colour of alternative exon 2}

\item{alternative2Stroke}{Character: stroke colour of alternative exon 2}
}
\value{
Diagrams per alternative splicing event in SVG
}
\description{
Prepare SVG diagram of alternative splicing events
}
\keyword{internal}
