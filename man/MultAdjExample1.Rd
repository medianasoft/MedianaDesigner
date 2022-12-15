\name{MultAdjExample1}
\alias{MultAdjExample1}
\title{
Simulation-based power calculations in Phase III trials with multiple dose-placebo comparisons
}
\description{
Consider a Phase III trial in patients with rheumatoid arthritis that will be conducted to investigate the efficacy and safety of three doses of an experimental treatment versus placebo. The three doses will be labeled as the low, medium and high doses. The primary efficacy endpoint is a binary endpoint based on the American College of Rheumatology response criteria (ACR50). The endpoint represents the proportion of patients who experience a symptomatic improvement after 26 weeks of treatment and thus a higher value of the endpoint indicates a beneficial effect. 

A multiplicity adjustment will be applied to address Type I error rate inflation cased by the evaluation of three dose-placebo comparisons. The adjustment will rely on the fixed-sequence procedure with the following sequence of dose-placebo comparisons: 
\itemize{
	\item High dose versus placebo.
	\item Medium dose versus placebo. 
	\item Low dose versus placebo.
}
Each comparison in this sequence will be evaluated at a one-sided alpha of 0.025 (provided the preceding comparison was significant) to control the overall Type I error rate at a one-sided 2.5\% level.

The following design parameters will be used in this example:
\itemize{
    \item A total of 440 patients will be enrolled in the trial and randomly assigned to placebo and three doses of the experimental treatment using an equal randomization approach. 
    \item The patient dropout rate at the end of the 26-week treatment period is expected be 20\%.
}

Power calculations will be performed in this trial using the \code{\link{MultAdj}} function. A list of all design and multiplicity adjustment parameters (\code{parameters}) needs to be created and then passed to this function. A simulation report can be created by calling the \code{\link{GenerateReport}} function. Power calculations can also be performed in this trial using a Shiny-based application with a graphical user interface that can be launched by calling the \code{\link{MultAdjApp}} function.
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
parameters$endpoint_type = "Binary"

# Direction of beneficial effect
parameters$direction = "Higher"

# Number of dose-placebo comparisons
parameters$n_comparisons = 3

# Number of endpoints
parameters$n_endpoints = 1

# Number of enrolled patients (placebo, low dose, medium dose and high dose) 
parameters$sample_size = c(110, 110, 110, 110)

# Patient dropout rate
parameters$dropout_rate = 0.2

# Endpoint information 
parameters$control_rate = 0.3
parameters$treatment_rate = c(0.5, 0.55, 0.6)

# Multiple testing procedure
parameters$mult_test = "Fixed-sequence"

# Hypothesis testing sequence (the null hypotheses of no effect will be 
# tested in the following sequence: H3 [high dose versus placebo], H2 [medium dose 
# versus placebo] and H1 [low dose versus placebo])
parameters$sequence = c(3, 2, 1)

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
              tempfile("MultAdj Multiple dose-placebo comparisons.docx", 
              fileext=".docx"))
}