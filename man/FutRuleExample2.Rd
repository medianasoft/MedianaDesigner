\name{FutRuleExample2}
\alias{FutRuleExample2}
\title{
Simulation-based selection of an optimal futility stopping rule (binary endpoint)
}
\description{
Consider a Phase II or Phase III trial with three treatment arms that correspond to two doses of an experimental treatment versus a control and assume that the primary efficacy endpoint is a binary endpoint (response rate, i.e., a higher value indicates a more favorable outcome). A futility assessment will be performed for each treatment arm at an interim look using conditional power and a treatment arm will be dropped if conditional power is below a pre-defined futility threshold. An optimal value of the futility threshold can be found by computing the sensitivity and specificity rates associated with the futility stopping rule and then identifying the threshold that simultaneously maximizes both rates.

The following design parameters will be assumed in the trial:
\itemize{
    \item A balanced design with 75 enrolled patients per trial arm will be utilized.
    \item The patient dropout rate at the end of the treatment period is assumed to be 5\%.
    \item An early interim analysis with the information fraction of 30\% will be conducted. 
}
The calculations will be performed under the following set of treatment effect assumptions: 
\itemize{
    \item The control response rate is set to 35\%. 
    \item The response rates in the three treatment arms are set to 45\%, 50\% and 55\%.
 }

Operating characteristics of the futility stopping rule (sensitivity and specificity rates) and an optimal futility threshold will be computed using the \code{\link{FutRule}} function based on 1,000 simulation runs. A list of all trial design parameters (\code{parameters}) needs to be set up as shown below and passed to this function. A detailed simulation report can be generated using the \code{\link{GenerateReport}} function and a graphical user interface can be launched by calling the \code{\link{FutRuleApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{FutRule}}, \code{\link{FutRuleExample1}}, \code{\link{FutRuleExample3}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Binary"

# Direction of favorable outcome
parameters$direction = "Higher"

# Number of enrolled patients (control, three treatments)
parameters$sample_size = c(75, 75, 75, 75)

# Dropout rate
parameters$dropout_rate = 0.05

# Response rate in the control arm 
parameters$control_rate = 0.35

# Response rates in the treatment arms 
parameters$treatment_rate = c(0.45, 0.5, 0.55)

# Information fraction
parameters$info_frac = 0.3

# One-sided alpha level
parameters$alpha = 0.025

# Number of simulations
parameters$nsims = 1000

# Run simulations to compute characteristics of the futility stopping rule
results = FutRule(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("FutRule Binary endpoint.docx", fileext=".docx"))
}