\name{MultAdjExample3}
\alias{MultAdjExample3}
\title{
Simulation-based power calculations in Phase III trials with multiple endpoints and multiple dose-placebo comparisons
}
\description{
The advanced setting with several sources of multiplicity, i.e., multiplicity induced by the analysis of multiple endpoints as well as multiple dose-placebo comparisons, will be illustrated using a Phase III trial for the treatment of acute schizophrenia. The efficacy profiles of two doses of an experimental treatment (Dose L and Dose H) will be compared to that of placebo using the following two endpoints:  
\itemize{
	\item Primary endpoint: Change in the Positive and Negative Syndrome Scale (PANSS) total score.
	\item Key secondary endpoint: Change in the Clinical Global Impressions severity (CGI-S) score.
}
The efficacy evaluations will be performed at the end of the 4-week treatment period. Both endpoints are normally distributed and a lower value of each endpoint indicates a more favorable outcome. The correlation between the endpoint-specific test statistics is assumed to be 0.5.

The resulting four null hypotheses of no effect will be grouped into two endpoint-specific families:
\itemize{
	\item Family 1: H1 (No difference between Dose L and placebo on the primary endpoint) and H2 (No difference between Dose H and placebo on the primary endpoint).
	\item Family 2: H3 (No difference between Dose L and placebo on the key secondary endpoint) and H4 (No difference between Dose H and placebo on the key secondary endpoint).
}
Clinically relevant logical restrictions will be imposed on the testing strategy to ensure that H3 can be tested only if H1 is rejected and H4 can be tested only if H2 is rejected.

The Hochberg-based gatekeeping procedure will be utilized in this trial to control the overall Type I error rate at a one-sided 2.5\% level across the four hypotheses of interest. This gatekeeping procedure will be set up using Hochberg-based component procedures, i.e., a truncated version of the Hochberg procedure with the truncation parameter of 0.8 will be applied in Family 1 and the regular Hochberg procedure in Family 2. The modified mixture method will be used in this gatekeeping procedure.

The following design parameters will be used in this example:
\itemize{
    \item A total of 360 patients will be enrolled in the trial. The patients will be allocated to the three trial arms using a 1:1:1 randomization scheme. 
    \item The patient dropout rate at the end of the 4-week treatment period is expected be 10\%.
}

To compute power for the individual hypotheses as well as each of the two families, the design parameters as well as the parameters of the selected gatekeeping procedure (\code{parameters}) will passed to the \code{\link{MultAdj}} function. A simulation report will be created by calling the \code{\link{GenerateReport}} function. Power calculations can also be performed in this trial using a Shiny-based application, which can be launched by calling the \code{\link{MultAdjApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{MultAdjApp}}, \code{\link{MultAdj}}, \code{\link{MultAdjExample2}}, \code{\link{MultAdjExample3}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Normal"

# Direction of beneficial effect
parameters$direction = "Lower"

# Number of dose-placebo comparisons
parameters$n_comparisons = 2

# Number of endpoints
parameters$n_endpoints = 2

# Number of enrolled patients (placebo, Dose L and Dose H) 
parameters$sample_size = c(120, 120, 120)

# Patient dropout rate
parameters$dropout_rate = 0.10

# Endpoint information (rows corresponds to endpoints and 
# columns corresponds to dose-placebo comparisons)
parameters$control_mean = c(0, 0)
parameters$treatment_mean = matrix(c(-0.30, -0.35, 
                                     -0.35, -0.40), 2, 2, byrow = TRUE)

parameters$control_sd = c(1, 1)
parameters$treatment_sd = matrix(c(1, 1, 
                                   1, 1), 2, 2, byrow = TRUE)

# Endpoint correlation matrix
parameters$endpoint_correlation = matrix(c(1, 0.5, 
                                           0.5, 1), 2, 2)

# Component procedure to be used in the gatekeeping procedure
parameters$mult_test = "Hochberg"

# Mixture method used in the gatekeeping procedure
parameters$mult_method = "Modified"

# Truncation parameters in the gatekeeping procedure
parameters$mult_test_gamma = c(0.8, 1)

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
              tempfile("MultAdj Multiple endpoints and dose-placebo comparisons.docx", 
              fileext=".docx"))
}