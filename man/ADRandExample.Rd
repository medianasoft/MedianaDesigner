\name{ADRandExample}
\alias{ADRandExample}
\title{
Simulation-based design of dose-finding Phase II trials with response-adaptive randomization (normally distributed endpoint)
}
\description{
Consider a Phase II trial to be conducted to evaluate the efficacy and safety profiles of several doses of an experimental treatment compared to placebo. The trial will be conducted for the treatment of cognitive impairment in patients with schizophrenia. The primary endpoint is a continuous, normally distributed endpoint (MCCB overall composite score) with a larger value indicating favorable effect. The treatment effect will be evaluated at 3 months after randomization, i.e., the treatment period will 3-month long. A multi-stage adaptive design will be employed in the trial to evaluate the dose-response relationship for the treatment while assigning most patients to the best performing doses. This will be accomplished using a response-adaptive randomization approach. 

The following design parameters will be used in this example:
\itemize{
    \item A total of 320 patients will be enrolled in the trial. 
    \item Patients will be randomly assigned to placebo and four doses of the experimental treatment (20 mg, 50 mg, 100 mg or 250 mg). 
    \item A four-stage design will be utilized with 80 patients enrolled in each stage. A balanced randomization approach will be used in the first stage of the trial and the randomization ratios in the dosing groups will be updated at the end of each stage. The randomization ratio for the placebo group will be fixed at 20\% throughout the trial. 
    \item Patients will be enrolled over the period of 24 months with the median enrollment time (time point by which 50\% of the patients will be enrolled into the trial) of 18 months.
    \item The patient dropout rate will be set to 10\% at the end of the treatment period.
}

Dose-response modeling at the end of each stage will be performed using the MCPMod approach based on four dose-response models (linear, exponential, Emax and logistic models). Non-linear parameters of these models will be specified as shown below (note that no parameter needs to be specified for the linear model). The rule for updating the randomization scheme will be based on the posterior probability of achieving a clinically meaningful improvement over placebo in each dosing group. The threshold for clinically meaningful improvement will be set at 1.5. 

Operating characteristics of the response-adaptive design will be evaluated using the \code{\link{ADRand}} function. As shown below, a list of all design and decision rule parameters (\code{parameters}) will be created and then passed to this function. It's important to note that, due to the use of Bayesian calculations, the simulation process will be rather slow. It is recommended to take advantage of parallel calculations by setting the number of cores (\code{ncores}) to a larger value as shown below. 

A simulation report for this adaptive design can be created by calling the \code{\link{GenerateReport}} function. The same functionality is also available via a Shiny-based application that can be launched by calling the \code{\link{ADRandApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{ADRand}}
}
\examples{
# List of all parameters

parameters = list()

# Endpoint type
parameters$endpoint_type = "Normal"

# Direction of favorable outcome
parameters$direction = "Higher"

# Dose levels in the trial 
parameters$dose_levels = c(0, 20, 50, 100, 250)

# Total number of enrolled patients in each trial stage
parameters$stage_sample_size = c(80, 80, 80, 80)

# Length of the treatment period (months)
parameters$treatment_period = 3

# Length of the patient enrollment period (months)
parameters$enrollment_period = 24

# Median enrollment time (months)
parameters$enrollment_parameter = 18

# Patient dropout rate
parameters$dropout_rate = 0.1

# Mean and SD for the primary endpoint in the placebo group
parameters$control_mean = 1.5
parameters$control_sd = 4

# Mean and SD for the primary endpoint in the dosing groups
parameters$treatment_mean = c(2, 2.5, 3, 3.5)
parameters$treatment_sd = c(4, 4, 4, 4)

# Fixed randomization ratio in the placebo group
parameters$ratio_placebo = 0.2

# Threshold for clinically meaningful improvement over placebo
parameters$delta = 1.5

# Balance parameter for adaptive randomization
parameters$balance = 2

# Non-linear parameters of the candidate dose-response models used in the MCPMod method 
# for modeling the dose-response relationship

# Non-linear parameter of the exponential model (delta)
parameters$exponential_model_parameter = 100

# Non-linear parameter of the Emax model (ED50)
parameters$emax_model_parameter = 200

# Non-linear parameters of the logistic model (ED50 and delta)
parameters$logistic_model_parameters = c(125, 25)

# One-sided Type I error rate
parameters$alpha = 0.025

# Number of cores for parallel calculations 
# (this is an example of how to run simulations on 3 cores)
# parameters$ncores = 3

# Number of simulations, you should prefer more
parameters$nsims = 30

# Remove this parameter in your code:
parameters$withoutCharts = TRUE

# Run simulations to compute operating characteristics
results = ADRand(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("ADRand Normally distributed endpoint.docx", fileext=".docx"))
}