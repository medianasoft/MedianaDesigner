\name{ADRand}
\alias{ADRand}
\alias{ADRandC}
\alias{ComputeDRFunctionParameters}
\alias{DRFunction}
\title{
Simulation-based design of adaptive trials with response-adaptive randomization
}
\description{
This function performs a simulation-based evaluation of operating characteristics for multi-stage adaptive Phase II trials that employ response-adaptive randomization. Patients are assigned to receive placebo or one of several doses of the experimental treatment. The efficacy data at the end of each stage are used to update the randomization scheme for the next stage with the ultimate goal of assigning most patients to the most effective doses. For an example of the function call, see \link{ADRandExample}.
}
\usage{
ADRand(parameters)
}
\arguments{
  \item{parameters}{List of the trial design and other parameters. The required elements are defined below: 
  \itemize{
    \item \code{endpoint_type}{: Character value defining the primary endpoint's type. Possible values:  
     \itemize{
      \item \code{"Normal"}: Normally distributed endpoint.
    }
    }
    \item \code{direction}{: Character value defining the direction of favorable outcome. Possible values: \code{"Higher"} (a higher value of the endpoint indicates a more favorable outcome) and \code{"Lower"} (a lower value of the endpoint indicates a more favorable outcome).}  
    \item \code{dose_levels}{: Integer vector defining the dose levels in the trial (0 corresponds to the placebo group). Each element must be non-negative.}
    \item \code{stage_sample_size}{: Integer vector defining the total number of patients enrolled in each stage. Each element must be positive.}
    \item \code{control_mean}{: Numeric value defining the mean of the primary endpoint in the placebo arm. }
    \item \code{control_sd}{: Numeric value defining the standard deviation of the primary endpoint in the placebo arm. This value must be positive. }
    \item \code{treatment_mean}{: Numeric vector defining the mean of the primary endpoint in each dosing arm.}
    \item \code{treatment_sd}{: Numeric vector defining the standard deviation of the primary endpoint in each dosing arm.}
    \item \code{treatment_period}{: Numeric value defining the length of the treatment period in the trial.}
    \item \code{ratio_placebo}{: Numeric value defining the fixed randomization ratio in the placebo arm. This value must be between 0 and 1.}
    \item \code{balance}{: Numeric value defining the degree of balance for adaptive randomization. This value must be between 0 and 3.}    
    \item \code{delta}{: Numeric value defining the threshold for a clinically meaningful improvement over placebo.}
    \item \code{exponential_model_parameter}{: Numeric value defining the non-linear parameter (delta) for the exponential dose-response model used in the MCPMod method. This value must be positive.}
    \item \code{emax_model_parameter}{: Numeric value defining the non-linear parameter (ED50) for the Emax dose-response model used in the MCPMod method. This value must be positive.}
    \item \code{logistic_model_parameters}{: Numeric vector with two elements defining the non-linear parameters (ED50 and delta) for the logistic dose-response model used in the MCPMod method. The values must be positive.}
    \item \code{enrollment_period}{: Numeric value defining the length of the patient enrollment period. This value must be positive. }
    \item \code{enrollment_parameter}{: Numeric value defining the median enrollment time. The patient enrollment process is assumed to be governed by a truncated exponential distribution and this parameter defines the time point by which 50\% of the patients are enrolled into the trial. This value must be between 0 and the length of the patient enrollment period.}
    \item \code{dropout_rate}{: Numeric value defining the patient dropout rate. A uniform patient dropout process is assumed and thus this parameter defines the fraction of patients that will be excluded from the interim and final analyses. This value must be between 0 and 1.}
    \item{alpha}{: Numeric value defining the overall one-sided Type I error rate. The default value is 0.025.}
    \item \code{random_seed}{: Integer value defining the random number generator seed. The default value is 49283.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    }
  }
}
\value{
The function returns an object of class \samp{ADRandResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{sim_results}{List with multiple data frames containing key descriptive statistics for each simulation run.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{ADRandApp}}
}
