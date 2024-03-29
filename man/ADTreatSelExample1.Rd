\name{ADTreatSelExample1}
\alias{ADTreatSelExample1}
\title{
Simulation-based design of an adaptive trial with treatment selection (normally distributed endpoint)
}
\description{
Consider a seamless Phase II/Phase III or Phase III trial with a normally distributed primary efficacy endpoint with a larger value indicating a more favorable outcome. The efficacy and safety profiles of two treatments, e.g., two doses of an experimental treatment, will be evaluated in this trial compared to a control. An adaptive design with two interim analyses will be employed in the trial (the first interim analysis may correspond to the end of the Phase II portion of the trial). The following decision rules will be applied at the interim looks:
\itemize{
	\item Futility stopping rule will be applied at the first interim analysis (IA1). A treatment arm may be dropped due to futility if the treatment-specific conditional power is too low. The trial will be terminated for futility if all treatment arms are dropped.
    \item Treatment selection rule will be applied at the second interim analysis (IA2). The best performing treatment will be selected for the final analysis (FA) and the remaining treatment will be dropped at this interim analysis.
} 
The following design parameters will be assumed:
\itemize{
    \item A balanced design with 200 enrolled patients per trial arm will be utilized in the trial.
    \item The patient dropout rate at the end of the treatment period is equal to 5\%.
    \item The information fractions at IA1 and IA2 are set to 40\% and 60\%, i.e., the first and second interim looks will be taken after 40\% and 60\% of the patients complete the treatment period or drop out of the trial before completing the treatment period. 
}
The futility threshold at IA1 is set to 25\%, which means that a treatment arm will be dropped at this interim analysis if conditional power does not exceed 25\%. The mean effect and standard deviation are assumed to be 0 and 1 in the control arm. The mean effects in the two treatment arms are assumed to be equal to 0.25 and 0.3, respectively, with a common standard deviation of 1.

To address multiplicity induced by the data-driven treatment selection, the Hochberg procedure will be applied at the final analysis.

Key operating characteristics of this adaptive design with a treatment selection option at the second interim look will be evaluated using the \code{\link{ADTreatSel}} function based on 10,000 simulation runs. To invoke this function, a list of all design and decision rule parameters (\code{parameters}) needs to be created as shown below. A comprehensive simulation report for this adaptive design can be generated by calling the \code{\link{GenerateReport}} function and a graphical user interface can be launched using the \code{\link{ADTreatSelApp}} function.
}
\arguments{No arguments}
\value{No return value}
\seealso{
  \code{\link{ADTreatSel}}, \code{\link{ADTreatSelExample2}}, \code{\link{ADTreatSelExample3}}
}
\examples{
# List of all parameters
parameters = list()

# Endpoint type
parameters$endpoint_type = "Normal"

# Direction of favorable outcome
parameters$direction = "Higher"

# Number of enrolled patients (control, two treatments)
parameters$sample_size = c(200, 200, 200)

# Patient dropout rate
parameters$dropout_rate = 0.05

# Mean and SD in the control arm 
parameters$control_mean = 0
parameters$control_sd = 1

# Means and SDs in the multiple treatment arms
parameters$treatment_mean = c(0.25, 0.30)
parameters$treatment_sd = c(1, 1)

# Information fractions at IA1, IA2, FA
parameters$info_frac = c(0.4, 0.6, 1)

# Futility threshold for conditional power at IA1
parameters$futility_threshold = 0.25

# Number of treatments to be selected at the second interim analysis
parameters$treatment_count = 1

# Multiple testing procedure to be used in the trial
parameters$mult_test = "Hochberg"

# One-sided alpha level
parameters$alpha = 0.025

# Number of simulations, you should prefer more
parameters$nsims = 100

# Number of cores for parallel calculations
parameters$ncores = 1

# Run simulations to compute operating characteristics
results = ADTreatSel(parameters)

# Generate a simulation report (remove tempfile)
GenerateReport(results, 
              tempfile("ADTreatSel Normally distributed endpoint.docx", fileext=".docx"))
}