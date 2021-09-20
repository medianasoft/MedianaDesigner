\name{MultAdj}
\alias{MultAdj}
\alias{MultAdjC}
\title{
Simulation-based design of traditional trials with multiple objectives
}
\description{
This function supports power calculations for a broad class of Phase III clinical trials with multiple objectives, including
\itemize{
    \item Trials with a single source of multiplicity, i.e., trials with a single endpoint and several dose-control comparisons, or two-arm trials with several endpoints.
    \item Trials with several sources of multiplicity, i.e., trials with several endpoints and several dose-placebo comparisons.
}
A fixed-sample design is assumed in this setting. Most commonly used multiplicity adjustments are supported for trials with a single source of multiplicity. In addition, global testing procedures can be applied in two-arm trials with several endpoints. Popular gatekeeping procedures are supported in the advanced multiplicity problems with several sources. For examples of the function call, see \link{MultAdjExample1}, \link{MultAdjExample2} or \link{MultAdjExample3}.
}
\usage{
MultAdj(parameters)
}
\arguments{
  \item{parameters}{List of the trial design and other parameters. The required elements are defined below: 
  \itemize{
    \item \code{endpoint_type}{: Character value defining the common type of trial endpoints. Possible values:  
     \itemize{
      \item \code{"Normal"}: Normally distributed endpoint.
      \item \code{"Binary"}: Binary endpoint.
    }
    }
    \item \code{direction}{: Character value defining the common direction of favorable outcome for all endpoints. Possible values: \code{"Higher"} (a higher value of each endpoint indicates a more favorable outcome) and \code{"Lower"} (a lower value of each endpoint indicates a more favorable outcome).}  
    \item \code{n_comparisons}{: Integer value defining the number of dose-control comparisons in the trial. This value must be positive.}
    \item \code{n_endpoints}{: Integer value defining the number of endpoints in the trial. This value must be positive. Either \code{n_comparisons} or \code{n_endpoints} must be greater than 1.}
    \item \code{sample_size}{: Integer vector defining the number of enrolled patients in each trial arm (control and experimental treatments). Each element must be positive.}
    \item \code{control_mean}{: Numeric vector defining the mean of each endpoint in the control arm. This parameter is required only with normally distributed endpoints (\code{endpoint_type="Normal"}).}
    \item \code{control_sd}{: Numeric vector defining the standard deviation of each endpoint in the control arm. Each element must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_mean}{: Numeric vector or matrix defining the mean of each endpoint in each experimental treatment arm. In clinical trials with several endpoints and several dose-placebo comparisons, the rows corresponds to the endpoints and the columns corresponds to the treatment-control comparisons. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_sd}{: Numeric vector or matrix defining the standard deviation of each endpoint in each experimental treatment arm. In clinical trials with several endpoints and several dose-placebo comparisons, the rows corresponds to the endpoints and the columns corresponds to the treatment-control comparisons.Each element must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{control_rate}{: Numeric vector defining the proportion or response rate for each endpoint in the control arm. Each element must be between 0 and 1. This parameter is required only with binary endpoints (\code{endpoint_type=} \code{"Binary"}).}
    \item \code{treatment_rate}{: Numeric vector or matrix defining the proportion or response rate for each endpoint in each experimental treatment arm. In clinical trials with several endpoints and several dose-placebo comparisons, the rows corresponds to the endpoints and the columns corresponds to the treatment-control comparisons. Each element must be between 0 and 1. This parameter is required only with binary endpoints.}
    \item \code{endpoint_correlation}{: Numeric matrix defining the pairwise correlations among the endpoint-specific test statistics. Each element must be between 0 and 1 and the matrix must be positive definite. This parameter is required only in trials with multiple endpoints.}
    \item \code{mult_test}{: Character value defining the multiple testing procedure, global testing procedure or gatekeeping procedure. Possible values:  
     \itemize{
      \item \code{"Bonferroni"}: Bonferroni multiple testing procedure.
      \item \code{"Holm"}: Holm multiple testing procedure.
      \item \code{"Fixed-sequence"}: Fixed-sequence multiple testing procedure.
      \item \code{"Chain"}: Chain multiple testing procedure.
      \item \code{"Hochberg"}: Hochberg multiple testing procedure in trials with a single source of multiplicity or Hochberg-based gatekeeping procedure in trials with several sources of multiplicity.
      \item \code{"Hommel"}: Hommel multiple testing procedure in trials with a single source of multiplicity or Hommel-based gatekeeping procedure in trials with several sources of multiplicity.
      \item \code{"O'Brien"}: O'Brien global testing procedure.
    }
    Note that the O'Brien procedure can be used only in two-arm trials with several endpoints, similarly gatekeeping procedures can be used only in trials with several endpoints and several dose-placebo comparisons.
    }
    \item \code{weights}{: Numeric vector defining the initial hypothesis weights. Each element must be between 0 and 1. This parameter is required only with multiple testing procedures.}
    \item \code{transition}{: Numeric matrix defining the hypothesis transition parameters. Each element must be between 0 and 1 and the sum of elements in each row must be less than or equal to 1. This parameter is required only with the chain multiple testing procedure.}
    \item \code{sequence}{: Integer vector defining the hypothesis testing sequence. This parameter is required only with the fixed-sequence multiple testing procedures.}
    \item \code{mult_method}{: Character value defining the mixture method for the gatekeeping procedure. Possible values:  
     \itemize{
      \item \code{"Standard"}: Standard mixture method.
      \item \code{"Modified"}: Modified mixture method.
      \item \code{"Enhanced"}: Enhanced mixture method.
    }
    This parameter is required only with gatekeeping procedures.
    }
    \item \code{mult_test_gamma}{: Numeric vector defining the truncation parameter for each endpoint-specific family of hypotheses. The vector's length must be equal to the number of endpoints. Each element must be between 0 and 1, the last element may be equal to 1 whereas the other elements must be strictly less than 1. This parameter is required only with gatekeeping procedures.}
    \item \code{dropout_rate}{: Numeric value defining the patient dropout rate. A uniform patient dropout process is assumed and thus this parameter defines the fraction of patients that will be excluded from the analysis. This value must be between 0 and 1.}
    \item{alpha}{: Numeric value defining the overall one-sided Type I error rate. The default value is 0.025.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    }
  }
}
\value{
The function returns an object of class \samp{MultAdjResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{sim_results}{Data frame containing the raw and adjusted p-values generated by the hypothesis tests for each simulation run. The first set of n columns correspond to the raw p-values for the n null hypotheses and the next set of n columns correspond to the adjusted p-values for the n null hypotheses.}
  \item{sim_summary}{List containing the power calculation results for the specified multiple testing procedure, global testing procedure or gatekeeping procedure.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{MultAdjApp}}, \code{\link{MultAdjExample1}}, \code{\link{MultAdjExample2}}, \code{\link{MultAdjExample3}}
}
