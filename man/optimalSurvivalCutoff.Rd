% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{optimalSurvivalCutoff}
\alias{optimalSurvivalCutoff}
\title{Calculate optimal data cutoff that best separates survival curves}
\usage{
optimalSurvivalCutoff(
  clinical,
  data,
  censoring,
  event,
  timeStart,
  timeStop = NULL,
  followup = "days_to_last_followup",
  session = NULL,
  filter = TRUE,
  survTime = NULL,
  lower = NULL,
  upper = NULL
)
}
\arguments{
\item{clinical}{Data frame: clinical data}

\item{data}{Numeric: data values}

\item{censoring}{Character: censor using \code{left}, \code{right},
\code{interval} or \code{interval2}}

\item{event}{Character: name of column containing time of the event of
interest}

\item{timeStart}{Character: name of column containing starting time of the
interval or follow up time}

\item{timeStop}{Character: name of column containing ending time of the
interval (only relevant for interval censoring)}

\item{followup}{Character: name of column containing follow up time}

\item{session}{Shiny session (only used for the visual interface)}

\item{filter}{Boolean or numeric: elements to use (all are used by default)}

\item{survTime}{\code{survTime} object: times to follow up, time start, time
stop and event (optional)}

\item{lower, upper}{Bounds in which to search (if \code{NULL}, bounds are set
to \code{lower = 0} and \code{upper = 1} if all data values are within that
interval; otherwise, \code{lower = min(data, na.rm = TRUE)} and
\code{upper = max(data, na.rm = TRUE)})}
}
\value{
List containing the optimal cutoff (\code{par}) and the corresponding
p-value (\code{value})
}
\description{
Uses \code{stats::optim} with the Brent method to test multiple cutoffs and
to find the minimum log-rank p-value.
}
\examples{
clinical <- read.table(text = "2549   NA ii  female
                                840   NA i   female
                                 NA 1204 iv    male
                                 NA  383 iv  female
                               1293   NA iii   male
                                 NA 1355 ii    male")
names(clinical) <- c("patient.days_to_last_followup",
                     "patient.days_to_death",
                     "patient.stage_event.pathologic_stage",
                     "patient.gender")
timeStart  <- "days_to_death"
event      <- "days_to_death"

psi <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
opt <- optimalSurvivalCutoff(clinical, psi, "right", event, timeStart)
}
\seealso{
Other functions to analyse survival: 
\code{\link{assignValuePerSubject}()},
\code{\link{getAttributesTime}()},
\code{\link{labelBasedOnCutoff}()},
\code{\link{plotSurvivalCurves}()},
\code{\link{plotSurvivalPvaluesByCutoff}()},
\code{\link{processSurvTerms}()},
\code{\link{survdiffTerms}()},
\code{\link{survfit.survTerms}()},
\code{\link{testSurvival}()}
}
\concept{functions to analyse survival}
