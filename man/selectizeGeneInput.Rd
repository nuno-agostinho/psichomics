% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{selectizeGeneInput}
\alias{selectizeGeneInput}
\title{Create input to select a gene}
\usage{
selectizeGeneInput(
  id,
  label = "Gene",
  choices = NULL,
  multiple = FALSE,
  ...,
  placeholder = "Type to search for a gene..."
)
}
\arguments{
\item{id}{Character: identifier}

\item{label}{Display label for the control, or \code{NULL} for no label.}

\item{choices}{List of values to select from. If elements of the list are
named, then that name --- rather than the value --- is displayed to the
user. It's also possible to group related inputs by providing a named list
whose elements are (either named or unnamed) lists, vectors, or factors. In
this case, the outermost names will be used as the group labels (leveraging
the \verb{<optgroup>} HTML tag) for the elements in the respective sublist. See
the example section for a small demo of this feature.}

\item{multiple}{Is selection of multiple items allowed?}

\item{...}{Arguments passed to the \code{options} list of
\code{selectizeInput()}}

\item{placeholder}{Character: placeholder}
}
\value{
HTML elements
}
\description{
Create input to select a gene
}
\keyword{internal}
