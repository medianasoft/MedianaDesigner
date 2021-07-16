\name{GenerateReport}
\alias{GenerateReport}
\alias{ReportDoc}
\alias{ADSSModReportDoc}
\alias{ADTreatSelReportDoc}
\alias{ADPopSelReportDoc}
\alias{EventPredReportDoc}
\alias{FutRuleReportDoc}
\alias{print.ADSSModResults}
\alias{print.ADTreatSelResults}
\alias{print.ADPopSelResults}
\alias{print.FutRuleResults}
\alias{print.EventPredResults}

\title{
Simulation report for adaptive trials  
}
\description{
This function creates a detailed simulation report in a Microsoft Word format. For examples of the function call, see \link{ADSSModExample1}, \link{ADTreatSelExample1}, \link{ADPopSelExample1}, \link{FutRuleExample1} and \link{EventPredExample}.
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

