% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_dimReduction_pca.R
\name{performPCA}
\alias{performPCA}
\title{Perform principal component analysis after processing missing values}
\usage{
performPCA(
  data,
  center = TRUE,
  scale. = FALSE,
  missingValues = round(0.05 * nrow(data)),
  ...
)
}
\arguments{
\item{data}{an optional data frame (or similar: see
    \code{\link[stats]{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}

\item{center}{a logical value indicating whether the variables
    should be shifted to be zero centered. Alternately, a vector of
    length equal the number of columns of \code{x} can be supplied.
    The value is passed to \code{scale}.}

\item{scale.}{a logical value indicating whether the variables should
    be scaled to have unit variance before the analysis takes
    place.  The default is \code{FALSE} for consistency with S, but
    in general scaling is advisable.  Alternatively, a vector of length
    equal the number of columns of \code{x} can be supplied.  The
    value is passed to \code{\link{scale}}.}

\item{missingValues}{Integer: number of tolerated missing values per column
to be replaced with the mean of the values of that same column}

\item{...}{Arguments passed on to \code{stats::prcomp}}
}
\value{
PCA result in a \code{prcomp} object
}
\description{
Perform principal component analysis after processing missing values
}
\examples{
performPCA(USArrests)
}
\seealso{
Other functions to analyse principal components: 
\code{\link{calculateLoadingsContribution}()},
\code{\link{plotPCA}()},
\code{\link{plotPCAvariance}()}
}
\concept{functions to analyse principal components}
