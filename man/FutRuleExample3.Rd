\name{FutRuleExample3}
\alias{FutRuleExample3}
\title{
Simulation-based selection of an optimal futility stopping rule (time-to-event endpoint)
}
\description{
Consider a Phase III trial with a time-to-event primary efficacy endpoint. A single regimen of an experimental treatment will be evaluated versus a control. A futility assessment will be conducted at an interim look using conditional power and the trial will be stopped for futility if conditional power in the single treatment arm is below a pre-defined futility threshold. An optimal value of the futility threshold can be found by computing the sensitivity and specificity rates associated with the futility stopping rule and then identifying the threshold that simultaneously maximizes both rates.

The following design parameters will be assumed in the trial:
\itemize{
    \item A balanced design with 125 enrolled patients per trial arm will be used.
    \item The target number of events at the final analysis is 175.
    \item The annual patient dropout rate is equal to 5\% (the time to patient dropout is assumed to follow an exponential distribution).
    \item The length of the patient enrollment period is 12 months and the median enrollment time is expected to be 9 months.
    \item A late interim analysis with the information fraction of 60\% will be conducted. 
}
The calculations will be performed assuming the median event times of 7.5 months and 10.5 months in the control and treatment arms, respectively. 

Operating characteristics of the futility stopping rule (sensitivity and specificity rates) and an optimal futility threshold will be computed using the \code{\link{FutRule}} function based on 1,000 simulation runs. A list of all trial design parameters (\code{parameters}) needs to be set up as shown below and passed to this function. A detailed simulation report can be generated using the \code{\link{GenerateReport}} function and a graphical user interface can be launched by calling the \code{\link{FutRuleApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{FutRule}}, \code{\link{FutRuleExample1}}, \code{\link{FutRuleExample2}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Time-to-event"

# Direction of favorable outcome
parameters$direction = "Higher"

# Number of enrolled patients (control, treatment)
parameters$sample_size = c(125, 125)

# Target event count at the final analysis
parameters$event_count = 175

# Annual patient dropout rate
parameters$dropout_rate = 0.05

# Median time in the control arm 
parameters$control_time = 7.5

# Median time in the treatment arm
parameters$treatment_time = 10.5

# Information fraction
parameters$info_frac = 0.6

# Enrollment period
parameters$enrollment_period = 12

# Median enrollment time
parameters$enrollment_parameter = 9

# One-sided alpha level
parameters$alpha = 0.025

# Number of simulations
parameters$nsims = 1000

# Run simulations to compute characteristics of the futility stopping rule
results = FutRule(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("FutRule Time-to-event endpoint.docx", fileext=".docx"))
}