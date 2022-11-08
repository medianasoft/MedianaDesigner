\name{ClustRand}
\alias{ClustRand}
\alias{ClustRandC}
\alias{ClustRandSingleCore}
\alias{ClustRandGLMEMR}
\alias{ClustRandGEEC}
\alias{ExportRandomClusterSize}
\title{
Simulation-based design of cluster-randomized trials
}
\description{
This function performs a simulation-based evaluation of operating characteristics for cluster-randomized trials with normally distributed or binary endpoints. A multi-arm clinical trial is assumed (control arm and several treatment arms). For examples of the function call, see \link{ClustRandExample1} or \link{ClustRandExample2}.
}
\usage{
ClustRand(parameters)
}
\arguments{
  \item{parameters}{List of the trial design and other parameters. The required elements are defined below: 
  \itemize{
    \item \code{endpoint_type}{: Character value defining the primary endpoint's type. Possible values:  
     \itemize{
      \item \code{"Normal"}: Normally distributed endpoint.
      \item \code{"Binary"}: Binary endpoint.
    }
    }
    \item \code{direction}{: Character value defining the direction of favorable outcome. Possible values: \code{"Higher"} (a higher value of the endpoint indicates a more favorable outcome) and \code{"Lower"} (a lower value of the endpoint indicates a more favorable outcome).}    
    \item \code{sample_size}{: Integer vector defining the number of completers in the trial arms (control arm and several treatment arms). Completers are defined as patients who complete the trial and are included in the final analysis. Each element must be positive.}
    \item \code{method_type}{: Character value defining the data analysis method. Possible values:  
     \itemize{
      \item \code{"GEE"}: Generalized estimating equations.
      \item \code{"GLMEM"}: Generalized linear mixed effects model. The current implementation of the GLMEM method is much less efficient than that for the GEE method and results in longer simulation run times.
    }
    }
    \item \code{cluster_scheme}{: Character value defining the cluster scheme. Possible values:  
     \itemize{
      \item \code{"Fixed"}: Cluster sizes are pre-defined.
      \item \code{"Random"}: Cluster sizes are randomly generated from a generalized Bernoulli distribution. 
    }
    }
    \item \code{control_cluster_size}{: Numeric vector defining the number of cluster sizes in the control arm. This parameter is required only if the cluster sizes are pre-defined (\code{cluster_scheme="Fixed"}).}
    \item \code{treatment_cluster_size}{: Numeric matrix defining the number of cluster sizes in the treatment arms (each row corresponds to a treatment arm). This parameter is required only if the cluster sizes are pre-defined.}
    \item \code{control_cluster_proportion}{: Numeric vector of relative cluster sizes in the control arm. The cluster membership is determined using a generalized Bernoulli distribution based on the relative cluster sizes. Each element must be between 0 and 1. This parameter is required only if the cluster sizes are random (\code{cluster_scheme="Random"}).}
    \item \code{treatment_cluster_proportion}{: Numeric matrix of relative cluster sizes in the treatment arms (each row corresponds to a treatment arm). The cluster membership in each treatment arm is determined using a generalized Bernoulli distribution based on the relative cluster sizes. Each element must be between 0 and 1. This parameter is required only if the cluster sizes are random.}
    \item \code{control_mean}{: Numeric value defining the mean of the primary endpoint in the control arm. This parameter is required only with normally distributed endpoints (\code{endpoint_type="Normal"}).}
    \item \code{treatment_mean}{: Numeric vector defining the means of the primary endpoint in the treatment arms. This parameter is required only with normally distributed endpoints.}
    \item \code{control_rate}{: Numeric value defining the proportion or response rate for the primary endpoint in the control arm. This value must be between 0 and 1. This parameter is required only with binary endpoints \cr (\code{endpoint_type="Binary"}).}
    \item \code{treatment_rate}{: Numeric vector defining the proportion or response rate for the primary endpoint in the treatment arms. Each element must be between 0 and 1. This parameter is required only with binary endpoints.}
    \item \code{control_icc}{: Numeric value defining the intra-cluster correlation coefficient in the control arm. This value must be between 0 and 1.}
    \item \code{control_between_cluster_sd}{: Numeric value defining the between-cluster standard deviation in the control arm. This value must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{treatment_icc}{: Numeric vector defining the intra-cluster correlation coefficients in the treatment arms. Each element must be between 0 and 1.}
    \item \code{treatment_between_cluster_sd}{: Numeric vector defining the between-cluster standard deviation in the treatment arms. Each element must be positive. This parameter is required only with normally distributed endpoints.}
    \item \code{descriptive_statistics}{: Logical value requesting the computation of descriptive statistics (arm-specific effects and cluster sizes if the cluster sizes are random) from each simulation run. The default value is FALSE.}
    \item{alpha}{: Numeric value defining the overall one-sided Type I error rate. The default value is 0.025.}
    \item \code{random_seed}{: Integer value defining the random number generator seed. The default value is 49283.}
    \item \code{nsims}{: Integer value defining the number of simulation runs.}
    \item \code{ncores}{: Integer value defining the number of cores for parallel calculations. The number of cores cannot exceed the maximum available number of cores. The default value is 1.} 
    }
  }
}
\value{
The function returns an object of class \samp{ClustRandResults}. This object is a list with the following components:
  \item{parameters}{List containing the user-specified parameters.}
  \item{pval_results}{Data frame containing the convergence flags and p-values for each simulation run.}
  \item{coef_results}{Data frame containing the convergence flags and model parameter estimares for each simulation run. The data frame is created if the descriptive statistics are requested (\code{descriptive_statistics=TRUE}).}  
  \item{cluster_size_results}{Data frame containing the cluster sizes for each simulation run if the cluster sizes are random (\code{cluster_scheme="Random"}). The data frame is created if the descriptive statistics are requested (\code{descriptive_statistics=TRUE}).}  
  \item{sim_summary}{List containing the key operating characteristics of the cluster-randomized design.}
  A detailed summary of the simulation results can be created using the \code{\link{GenerateReport}} function.
}
\seealso{
  \code{\link{ClustRandApp}}
}
