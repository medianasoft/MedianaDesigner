\name{MedianaDesigner-package}
\alias{MedianaDesigner-package}
\alias{MedianaDesigner}
\docType{package}
\title{
Efficient Simulation-Based Power and Sample Size Calculations for a Broad Class of Late-Stage Clinical Trials
}
\description{
The MedianaDesigner package provides efficient simulation-based tools for Phase II/III and Phase III clinical trials. The following five modules are included in the package:
\itemize{
\item Adaptive designs with data-driven sample size or event count re-estimation.
\item Adaptive designs with data-driven treatment selection.
\item Adaptive designs with data-driven population selection.
\item Optimal selection of a futility stopping rule.
\item Event prediction in event-driven trials.
}
For more information on this package, visit \url{https://mediana.us/free-software/}. The technical manuals with a detailed description of the statistical methodology implemented in each module can be found at \url{https://medianasoft.github.io/MedianaDesigner}.
}
\details{
\tabular{ll}{
Package: \tab MedianaDesigner\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2021-07-07\cr
License: \tab GPL-3
}
Key functions included in the package:
\itemize{
\item \code{\link{ADSSMod}}: Simulation-based design of adaptive trials with data-driven sample size or event count re-estimation.
\item \code{\link{ADSSModApp}}: Graphical user interface to design adaptive trials with data-driven sample size or event count re-estimation.
\item \code{\link{ADTreatSel}}: Simulation-based design of adaptive trials with data-driven treatment selection.
\item \code{\link{ADTreatSelApp}}: Graphical user interface to design adaptive trials with data-driven treatment selection.
\item \code{\link{ADPopSel}}: Simulation-based design of adaptive trials with data-driven population selection.
\item \code{\link{ADPopSelApp}}: Graphical user interface to design adaptive trials with data-driven population selection.
\item \code{\link{FutRule}}: Simulation-based selection of an optimal futility stopping rule at an interim analysis.
\item \code{\link{FutRuleApp}}: Graphical user interface for an optimal selection of a futility stopping rule.
\item \code{\link{EventPred}}: Simulation-based event prediction in trials with an event-driven design.
\item \code{\link{EventPredPriorDistribution}}: Calculation of the parameters of prior gamma distributions.
\item \code{\link{EventPredApp}}: Graphical user interface for event prediction.
\item \code{\link{GenerateReport}}: Simulation report for any module.
}
The package comes with one example data set:
\itemize{
\item \code{\link{EventPredData}}: Example data set to illustrate event prediction in event-driven trials.
}
Multiple case studies are included in this vignette to illustrate the use of this package for designing clinical trials with data-driven design elements:
\itemize{
\item \link{ADSSModExample1}: Simulation-based design of an adaptive trial with sample size re-estimation (normally distributed endpoint).
\item \link{ADSSModExample2}: Simulation-based design of an adaptive trial with sample size re-estimation (binary endpoint).
\item \link{ADSSModExample3}: Simulation-based design of an adaptive trial with event count re-estimation (time-to-event endpoint).
\item \link{ADTreatSelExample1}: Simulation-based design of an adaptive trial with treatment selection (normally distributed endpoint).
\item \link{ADTreatSelExample2}: Simulation-based design of an adaptive trials with treatment selection (binary endpoint).
\item \link{ADTreatSelExample3}: Simulation-based design of an adaptive trials with treatment selection (time-to-event endpoint).
\item \link{ADPopSelExample1}: Simulation-based design of an adaptive trial with population selection (normally distributed endpoint).
\item \link{ADPopSelExample2}: Simulation-based design of an adaptive trials with population selection (binary endpoint).
\item \link{ADPopSelExample3}: Simulation-based design of an adaptive trials with population selection (time-to-event endpoint).
\item \link{FutRuleExample1}: Simulation-based selection of an optimal futility stopping rule (normally distributed endpoint).
\item \link{FutRuleExample2}: Simulation-based selection of an optimal futility stopping rule (binary endpoint).
\item \link{FutRuleExample3}: Simulation-based selection of an optimal futility stopping rule (time-to-event endpoint).
\item \link{EventPredExample}: Simulation-based event prediction in trials with an event-driven design (time-to-event endpoint).

}

}
\author{
Mediana (www.mediana.us)
}
\references{
Bagiella, E., Heitjan, D.F. (2001). Predicting analysis times in randomized clinical trials. \emph{Statistics in Medicine}. 20, 2055-2063.

Chuang-Stein, C., Kirby, S., French, J., Kowalski, K., Marshall, S., Smith, M. K. (2011). A quantitative approach for making go/no-go decisions in drug development. \emph{Drug Information Journal}. 45, 187-202.

Dmitrienko, A., Paux, G. (2017). Subgroup analysis in clinical trials. \emph{Clinical Trial Optimization Using R}. Dmitrienko, A., Pulkstenis, E. (editors). Chapman and Hall/CRC Press, New York.

Dmitrienko, A., D'Agostino, R.B. (2018). Multiplicity considerations in clinical trials. \emph{New England Journal of Medicine}. 378, 2115-2122.

Herson, J. (1979). Predictive probability early termination plans for Phase II clinical trials. \emph{Biometrics}. 35, 775-783.

Millen, B., Dmitrienko, A., Ruberg, S., Shen, L. (2012). A statistical framework for decision making in confirmatory multi-population tailoring clinical trials. \emph{Drug Information Journal}. 46, 647-656.

Wang, D., Cui, L., Zhang, L., Yang, B. (2014). An ROC approach to evaluate interim go/no-go decision-making quality with application to futility stopping in the clinical trial designs. \emph{New Developments in Statistical Modeling, Inference and Application}. Jin, Z., Liu, M., Luo, X. (editors). Springer, New York. 121-147.

Wassmer, G., Brannath, W. (2016). \emph{Group Sequential and Confirmatory Adaptive Designs in Clinical Trials}. New York: Springer.

}
\keyword{ package }
