\name{EventPredPriorDistribution}
\alias{EventPredPriorDistribution}
\title{
Calculation of the parameters of prior gamma distributions
}
\description{
This function can be used for calculating the parameters of prior gamma distributions for the event hazard rate, patient dropout hazard rate and patient enrollment rate in the \code{\link{EventPred}} function. For an example of the function call, see \link{EventPredExample}.
}
\usage{
EventPredPriorDistribution(expected, uncertainty)
}
\arguments{
  \item{expected}{Numeric value defining the expected value of the event hazard rate, patient dropout hazard rate or patient enrollment rate. This value must be positive.}
  \item{uncertainty}{Numeric value defining the amount of uncertainly around the expected value defined above. As a general guideline, the uncertainty parameters of 0.1, 0.3 and 0.5 define high-confidence, medium-confidence and low-confidence prior distributions. This value must be between 0 and 1.}

}
\value{
The function returns a numeric vector with the shape and rate parameters of the prior gamma distribution.}
\seealso{
  \code{\link{EventPred}}
}
