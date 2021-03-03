\name{ADPopSel}
\alias{ADPopSel}
\alias{ADPopSelC}
\title{
Simulation-based design of adaptive trials with data-driven population selection
}
\description{
This function performs a simulation-based evaluation of operating characteristics for adaptive trials with data-driven population selection. A two-arm Phase III clinical trial with two pre-defined patient populations (overall population and subpopulation of biomarker-positive patients) and two interim analyses is assumed. The first interim analysis supports early stopping for futility in the overall population and the second interim analysis identifies the most promising patient population or populations. For examples of the function call, see \link{ADPopSelExample1}, \link{ADPopSelExample2} or \link{ADPopSelExample3}.
}
\usage{
ADPopSel(parameters)

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
    \item \code{sample_size}{: Integer vector with two elements defining the number of enrolled patients in the two trial arms (control and experimental treatment). Each element must be positive.}
    \item \code{control_mean}{: Numeric value defining the mean of the primary endpoint in the control arm. This parameter is required only with normally distributed endpoints (\code{endpoint_type="Normal"}).}
    \item \code{control_sd}{: Numeric value defining the standard deviation of the primary endpoint in the control arm. This value must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_mean}{: Numeric value defining the mean of the primary endpoint in the experimental treatment arm. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_sd}{: Numeric value defining the standard deviation of the primary endpoint in the experimental treatment arm. This value must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{control_rate}{: Numeric value defining the proportion or response rate for the primary endpoint in the control arm. This value must be between 0 and 1. This parameter is required only with binary endpoints \cr (\code{endpoint_type="Binary"}).}
    \item \code{treatment_rate}{: Numeric value defining the proportion or response rate for the primary endpoint in the experimental treatment arm. This value must be between 0 and 1. This parameter is required only with binary endpoints.}
    \item \code{control_time}{: Numeric value defining the median time, e.g., the median survival time, for the exponentially distributed primary endpoint in the control arm. This value must be positive. This parameter is required only with time-to-event endpoints (\code{endpoint_type="Time-to-event"}).}
    \item \code{treatment_time}{: Numeric value defining the median time, e.g., the median survival time, for the exponentially distributed primary endpoint in the experimental treatment arm. This value must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{info_frac}{: Numeric vector with three elements defining the information fractions at the first interim analysis, at the second interim analysis and at the final analysis. The first and second elements must be between 0 and 1, and the third element must be 1.}
    \item \code{event_count}{: Numeric value defining the total number of events at the final analysis. This value must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{futility_threshold}{: Numeric value defining the futility threshold for conditional power at the first interim analysis. This value must be positive.}
    \item \code{influence}{: Numeric value defining the influence threshold for selecting the most promising population or populations at the second interim analysis. This value must be positive.}
    \item \code{interaction}{: Numeric value defining the interaction threshold for selecting the most promising population or populations at the second interim analysis. This value must be greater than 1.}
    \item \code{dropout_rate}{: Numeric value defining the patient dropout rate. With normally distributed endpoints and binary endpoints, a uniform patient dropout process is assumed and thus this parameter defines the fraction of patients that will be excluded from the interim and final analyses. With time-to-event endpoints, the time to patient dropout is assumed to follow an exponential distribution and this parameter defines the annual dropout rate. This value must be between 0 and 1.}
    \item \code{enrollment_period}{: Numeric value defining the length of the patient enrollment period in months. This value must be positive. This parameter is required only with time-to-event endpoints.}
    \item \code{enrollment_parameter}{: Numeric value defining the median enrollment time in months. The patient enrollment process is assumed to be governed by a truncated exponential distribution and this parameter defines the time point by which 50\% of the patients are enrolled into the trial. This value must be between 0 and the length of the patient enrollment period. This parameter is required only with time-to-event endpoints.}.
    \item \code{alpha}{: Numeric value defining the overall one-sided Type I error rate. The default value is 0.025.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    }
  }
}
\value{
The function returns an object of class \samp{ADPopSelResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{sim_results}{Data frame containing key descriptive statistics for each simulation run.}
  \item{sim_summary}{List containing the key operating characteristics of the adaptive design.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{ADPopSelApp}}
}
