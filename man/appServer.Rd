% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app.R, R/analysis.R, R/analysis_correlation.R,
%   R/analysis_diffExpression.R, R/analysis_diffExpression_event.R,
%   R/analysis_diffExpression_table.R, R/analysis_diffSplicing.R,
%   R/analysis_diffSplicing_event.R, R/analysis_diffSplicing_table.R,
%   R/analysis_dimReduction.R, R/analysis_dimReduction_ica.R,
%   R/analysis_dimReduction_pca.R, R/analysis_information.R,
%   R/analysis_survival.R, R/analysis_template.R, R/data.R, R/data_firebrowse.R,
%   R/data_geNormalisationFiltering.R, R/data_gtex.R, R/data_inclusionLevels.R,
%   R/data_inclusionLevelsFilter.R, R/data_local.R, R/data_recount.R,
%   R/groups.R, R/help.R
\name{appServer}
\alias{appServer}
\alias{analysesServer}
\alias{diffEventServer}
\alias{correlationServer}
\alias{diffExpressionServer}
\alias{diffExpressionEventServer}
\alias{diffExpressionTableServer}
\alias{diffSplicingServer}
\alias{diffSplicingEventServer}
\alias{diffSplicingTableServer}
\alias{dimReductionServer}
\alias{icaServer}
\alias{pcaServer}
\alias{infoServer}
\alias{survivalServer}
\alias{templateServer}
\alias{dataServer}
\alias{firebrowseServer}
\alias{geNormalisationFilteringServer}
\alias{gtexDataServer}
\alias{inclusionLevelsServer}
\alias{inclusionLevelsFilterServer}
\alias{localDataServer}
\alias{recountDataServer}
\alias{groupsServer}
\alias{helpServer}
\title{Server logic}
\usage{
appServer(input, output, session)

analysesServer(input, output, session)

diffEventServer(ns, input, output, session, psi)

correlationServer(input, output, session)

diffExpressionServer(input, output, session)

diffExpressionEventServer(input, output, session)

diffExpressionTableServer(input, output, session)

diffSplicingServer(input, output, session)

diffSplicingEventServer(input, output, session)

diffSplicingTableServer(input, output, session)

dimReductionServer(input, output, session)

icaServer(input, output, session)

pcaServer(input, output, session)

infoServer(input, output, session)

survivalServer(input, output, session)

templateServer(input, output, session)

dataServer(input, output, session)

firebrowseServer(input, output, session)

geNormalisationFilteringServer(input, output, session)

gtexDataServer(input, output, session)

inclusionLevelsServer(input, output, session)

inclusionLevelsFilterServer(input, output, session)

localDataServer(input, output, session)

recountDataServer(input, output, session)

groupsServer(input, output, session)

helpServer(input, output, session)
}
\arguments{
\item{input}{Shiny input}

\item{output}{Shiny output}

\item{session}{Shiny session}
}
\value{
\code{NULL} (function is only used to modify the Shiny session's
state or internal variables)
}
\description{
Instructions to build the Shiny app
}
\keyword{internal}
