\name{FutRuleExample1}
\alias{FutRuleExample1}
\title{
Simulation-based selection of an optimal futility stopping rule (normally distributed endpoint)
}
\description{
Consider a multi-arm Phase II or Phase III trial with a single interim analysis and assume that the primary efficacy endpoint is normally distributed. A futility assessment will be performed for each treatment arm at this interim look using the predicted probability of success (conditional power). A treatment arm will be dropped at this decision point if conditional power is below a pre-defined threshold. An optimal value of this futility threshold can be found by computing the sensitivity and specificity rates associated with the futility stopping rule and then identifying the threshold that simultaneously maximizes both rates.

The following design parameters will be assumed in the trial:
\itemize{
    \item Two treatments (two doses of an experimental treatment) will be compared to a control in this trial and a balanced design with 100 enrolled patients per trial arm will be utilized.
    \item The patient dropout rate at the end of the treatment period is assumed to be 5\%.
    \item The information fraction at the interim analysis is 50\%. 
}
The calculations will be performed under the following set of treatment effect assumptions: 
\itemize{
    \item The mean and standard deviation of the primary efficacy endpoint in the control arm are set to 0 and 1, respectively. 
    \item The mean values in the two treatment arms are assumed to be equal to 0.25 and 0.3 with a common standard deviation of 1.
 }

Operating characteristics of the futility stopping rule such as the sensitivity and specificity rates as well as an optimal futility threshold will be computed using the \code{\link{FutRule}} function based on 1,000 simulation runs. A list of all trial design parameters (\code{parameters}) needs to be set up as shown below and passed to this function. A detailed simulation report can be generated using the \code{\link{GenerateReport}} function and a graphical user interface can be launched by calling the \code{\link{FutRuleApp}} function.
}
\seealso{
  \code{\link{FutRule}}, \code{\link{FutRuleExample2}}, \code{\link{FutRuleExample3}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Normal"

# Number of enrolled patients (control, two treatments)
parameters$sample_size = c(100, 100, 100)

# Patient dropout rate
parameters$dropout_rate = 0.05

# Mean and SD in the control arm 
parameters$control_mean = 0
parameters$control_sd = 1

# Means and SDs in the treatment arms 
parameters$treatment_mean = c(0.25, 0.30)
parameters$treatment_sd = c(1, 1)

# Information fraction
parameters$info_frac = 0.5

# One-sided alpha level
parameters$alpha = 0.025

# Number of simulations
parameters$nsims = 1000

# Run simulations to compute characteristics of the futility stopping rule
results = FutRule(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("FutRule Normally distributed endpoint.docx", fileext=".docx"))
}