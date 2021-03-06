% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_information.R
\name{ensemblToUniprot}
\alias{ensemblToUniprot}
\title{Convert from Ensembl to UniProt identifier}
\usage{
ensemblToUniprot(protein)
}
\arguments{
\item{protein}{Character: Ensembl identifier}
}
\value{
UniProt protein identifier
}
\description{
Convert from Ensembl to UniProt identifier
}
\examples{
gene <- "ENSG00000173262"
ensemblToUniprot(gene)

protein <- "ENSP00000445929"
ensemblToUniprot(protein)
}
\seealso{
Other functions to retrieve external information: 
\code{\link{plotProtein}()},
\code{\link{plotTranscripts}()},
\code{\link{queryEnsemblByGene}()}
}
\concept{functions to retrieve external information}
