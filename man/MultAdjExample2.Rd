\name{MultAdjExample2}
\alias{MultAdjExample2}
\title{
Simulation-based power calculations in Phase III trials with multiple endpoints
}
\description{
Consider a Phase III trial for the treatment of agitation in patients with dementia. A single dose of an experimental treatment will be compared to placebo in this trial. The efficacy profile of the experimental treatment will be assessed using three primary endpoints, namely, three components of the Cohen-Mansfield Agitation Inventory, Community version (CMAI-C) rating scale: 
 \itemize{
    \item Aggressive behavior.
    \item Non-aggressive agitated behavior. 
    \item Verbally agitated behavior.
}
Changes in the individual endpoints will be assessed at Day 28. Each endpoint is normally distributed and, since higher scores on this rating scale indicate more severe symptoms, a lower value of each endpoint corresponds to a beneficial effect. The endpoints are assumed to be equally correlated with a common correlation coefficient of 0.4.

An important feature of this trial is that the trial will meet its primary objective if at least one endpoint test is significant. As a result, the overall Type I error rate will be inflated if the primary analysis is performed without a multiplicity adjustment. To preserve the overall Type I error rate at a one-sided 2.5\% level, the Holm procedure will be applied. The three endpoints are considered to be equally important and, by default, equal weights will be assigned to the endpoint-specific null hypotheses of no treatment effect.

The following design parameters will be used in this example:
\itemize{
    \item A total of 380 patients will be enrolled in the trial to be assigned to either placebo or experimental treatment using a 1:1 ratio. 
    \item The patient dropout rate at the end of the treatment period is expected be 5\%.
}

The \code{\link{MultAdj}} function will be invoked to run power calculations in this trial. A list of all design and multiplicity adjustment parameters (\code{parameters}) will be created and passed to this function. A simulation report will be created by calling the \code{\link{GenerateReport}} function. Power calculations can also be performed in this trial using a Shiny-based application that can be launched by calling the \code{\link{MultAdjApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{MultAdjApp}}, \code{\link{MultAdj}}, \code{\link{MultAdjExample1}}, \code{\link{MultAdjExample3}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Normal"

# Direction of beneficial effect
parameters$direction = "Lower"

# Number of dose-placebo comparisons
parameters$n_comparisons = 1

# Number of endpoints
parameters$n_endpoints = 3

# Number of enrolled patients (placebo, experimental treatment) 
parameters$sample_size = c(190, 190)

# Patient dropout rate
parameters$dropout_rate = 0.05

parameters$control_mean = c(0, 0, 0)
parameters$treatment_mean = c(-0.3, -0.3, -0.3)
parameters$control_sd = c(1, 1, 1)
parameters$treatment_sd = c(1, 1, 1)

# Endpoint correlation matrix
parameters$endpoint_correlation = matrix(c(1, 0.4, 0.4, 
                                           0.4, 1, 0.4,
                                           0.4, 0.4, 1), 3, 3)

# Multiple testing procedure
parameters$mult_test = "Holm"

# Overall one-sided Type I error rate
parameters$alpha = 0.025

# Number of simulations
parameters$nsims = 100

# Number of cores for parallel calculations
parameters$ncores = 1

# Run simulations to perform power calculations
results = MultAdj(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("MultAdj Multiple endpoints.docx", 
              fileext=".docx"))
}