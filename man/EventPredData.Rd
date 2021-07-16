\name{EventPredData}
\alias{EventPredData}
\docType{data}
\title{Example data set for EventPred}
\description{
Example data set based on a Phase III oncology trial with an overall survival endpoint to illustrate event prediction in event-driven trials.
}
\arguments{No arguments}
\value{No return value}
\usage{data(EventPredData)}
\format{
  A data set with 150 records (one record per patient) and 4 variables:
  \describe{
    \item{\code{enrollment}}{Time of the patient's enrollment (months).}
    \item{\code{time}}{Time of the patient's death or last contact if the event is censored relative to the patient's enrollment (months).}
    \item{\code{event}}{Event indicator (1 if the patient died and 0 otherwise).}
    \item{\code{dropout}}{Patient dropout indicator (1 if the patient dropped out of the trial/was lost to follow up and 0 otherwise).}
    }
}
\keyword{datasets}
