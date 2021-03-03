\name{FutRule}
\alias{FutRule}
\alias{FutRuleC}
\title{
Simulation-based selection of an optimal futility stopping rule at an interim analysis
}
\description{
This function evaluates operating characteristics of a multi-arm trial design with a single interim analysis. A futility stopping rule will be applied at this interim look and the trial will be stopped early for futility if the predicted probability of success (conditional power) is less than a pre-defined futility threshold in all treatment arms. An optimal value of the futility threshold is derived by maximizing the sensitivity and specificity rates. For examples of the function call, see \link{FutRuleExample1}, \link{FutRuleExample2} or \link{FutRuleExample3}.
}
\usage{
FutRule(parameters)

}
\arguments{
  \item{parameters}{List of the trial design and other parameters. The required elements are defined below:
  \itemize{
    \item \code{endpoint_type}{: Character value defining the primary endpoint's type. Possible values:  
     \itemize{
      \item \code{"Normal"}: Normally distributed endpoint.
      \item \code{"Binary"}: Binary endpoint.
      \item \code{"Time-to-event"}: Time-to-event endpoint.
    }
    }
    \item \code{sample_size}{: Integer vector defining the number of enrolled patients in the trial arms (control and multiple experimental treatments). Each element must be positive.}
    \item \code{control_mean}{: Numeric value defining the mean of the primary endpoint in the control arm. This parameter is required only with normally distributed endpoints (\code{endpoint_type="Normal"}).}
    \item \code{control_sd}{: Numeric value defining the standard deviation of the primary endpoint in the control arm. This value must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_mean}{: Numeric vector defining the means of the primary endpoint in the experimental treatment arms. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_sd}{: Numeric vector defining the standard deviations of the primary endpoint in the experimental treatment arms. Each element must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{control_rate}{: Numeric value defining the proportion or response rate for the primary endpoint in the control arm. This value must be between 0 and 1. This parameter is required only with binary endpoints \cr (\code{endpoint_type="Binary"}).}
    \item \code{treatment_rate}{: Numeric vector defining the proportions or response rates for the primary endpoint in the experimental treatment arms. Each element must be between 0 and 1. This parameter is required only with binary endpoints.}
    \item \code{control_time}{: Numeric value defining the median time, e.g., the median survival time, for the exponentially distributed primary endpoint in the control arm. This value must be positive. This parameter is required only with time-to-event endpoints (\code{endpoint_type="Time-to-event"}).}
    \item \code{treatment_time}{: Numeric vector defining the median times, e.g., the median survival times, for the exponentially distributed primary endpoint in the experimental treatment arms. Each element must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{info_frac}{: Numeric value defining the information fraction at the interim analysis. This value must be between 0 and 1.}
    \item \code{event_count}{: Numeric value defining the total number of events at the final analysis. This value must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{dropout_rate}{: Numeric value defining the patient dropout rate. With normally distributed endpoints and binary endpoints, a uniform patient dropout process is assumed and thus this parameter defines the fraction of patients that will be excluded from the interim and final analyses. With time-to-event endpoints, the time to patient dropout is assumed to follow an exponential distribution and this parameter defines the annual dropout rate. This value must be between 0 and 1.}
    \item \code{enrollment_period}{: Numeric value defining the length of the patient enrollment period in months. This value must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{enrollment_parameter}{: Numeric value defining the median enrollment time in months. The patient enrollment process is assumed to be governed by a truncated exponential distribution and this parameter defines the time point by which 50\% of the patients are enrolled into the trial. This value must be between 0 and the length of the patient enrollment period. This parameter is required only with time-to-event endpoints.}.
    \item \code{alpha}{: Numeric value defining the overall one-sided Type I error rate. The default value is 0.025.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    }
  }
}
\value{
The function returns an object of class \samp{FutRuleResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{sim_summary}{List containing the operating characteristics of the futility stopping rule.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{FutRuleApp}}
}
