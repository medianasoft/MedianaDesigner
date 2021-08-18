\name{EventPred}
\alias{EventPred}
\alias{EventPredR}
\alias{EventPredEventCount}
\title{
Simulation-based event prediction in trials with an event-driven design
}
\description{
This function performs event forecasting in trials with an event-driven design. Blinded event data at an interim analysis are used to forecast the number of events at pre-defined time points in the future. For an example of the function call, see \link{EventPredExample}.
}
\usage{
EventPred(parameters)
}
\arguments{
  \item{parameters}{List of the trial design and other parameters. The required elements are defined below: 
  \itemize{
    \item \code{data_set}{: Data frame that contains blinded patient enrollment, event and dropout data at the interim analysis with a single record per patient. The following :  
     \itemize{
      \item \code{enrollment}{: Time of the patient's enrollment.}
      \item \code{time}{: Time of the patient's death or last contact if the event is censored relative to the patient's enrollment.}
      \item \code{event}{: Event indicator (1 if the patient died and 0 otherwise).}
      \item \code{dropout}{: Patient dropout indicator (1 if the patient dropped out of the trial/was lost to follow up and 0 otherwise).}
    }
    An example is provided in the \code{\link{EventPredData}} data set.
    }
    \item \code{time_points}{: Numeric vector defining the future time points for computing event predictions. Each elements must be greater than the latest time point in the \code{data_set}.}
    \item \code{event_prior_distribution}{: Numeric vector with two elements defining the shape and rate parameters of the prior gamma distribution for the event hazard rate. Each element must be positive. The parameters could be computed using the \code{\link{EventPredPriorDistribution}} function.}
    \item \code{dropout_prior_distribution}{: Numeric vector with two elements defining the shape and rate parameters of the prior gamma distribution for the patient dropout hazard rate. Each element must be positive. The parameters could be computed using the \code{\link{EventPredPriorDistribution}} function.}
    \item \code{enrollment_prior_distribution}{: Numeric vector with two elements defining the shape and rate parameters of the prior gamma distribution for the patient enrollment rate. Each element must be positive. The parameters could be computed using the \code{\link{EventPredPriorDistribution}} function.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    }
  }
}
\value{
The function returns an object of class \samp{EventPredResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{interim_analysis}{Data frame containing the event information from the interim analysis data set.}
  \item{predictions}{Data frame containing the event prediction information at the future time points.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{EventPredApp}}
}
