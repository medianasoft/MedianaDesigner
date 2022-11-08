\name{GenerateReport}
\alias{GenerateReport}
\alias{ReportDoc}
\alias{ADSSModReportDoc}
\alias{ADTreatSelReportDoc}
\alias{ADPopSelReportDoc}
\alias{EventPredReportDoc}
\alias{FutRuleReportDoc}
\alias{ADRandReportDoc}
\alias{MultAdjReportDoc}
\alias{ClustRandReportDoc}
\alias{print.ADSSModResults}
\alias{print.ADTreatSelResults}
\alias{print.ADPopSelResults}
\alias{print.FutRuleResults}
\alias{print.EventPredResults}
\alias{print.ADRandResults}
\alias{print.MultAdjResults}
\alias{print.ClustRandResults}
\title{
Simulation report  
}
\description{
This function creates a detailed simulation report in a Microsoft Word format. For examples of the function call, see \link{ADSSModExample1}, \link{ADTreatSelExample1}, \link{ADPopSelExample1}, \link{FutRuleExample1}, \link{EventPredExample} and \link{ADRandExample}.
}
\usage{
GenerateReport(results, report_filename)
}
\arguments{
  \item{results}{Object created by a simulation function, e.g., \code{\link{ADSSMod}}.  
  }
  \item{report_filename}{Character value defining the report's filename. The report is saved in the current working directory.
  }
}
\value{No return value, called for side effects}

