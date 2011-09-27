\name{summary.dhglm}
\Rdversion{1.1}
\alias{summary.dhglm}
%\alias{dhglm}

\docType{package}

\title{Summary Method for dhglmfit Object}

\description{It provides standard summary statistics for the fitted dhglmfit objects.}

\usage{
\method{summary}{dhglm}(object, \dots)
}

\arguments{
  \item{object}{A \code{\link{dhglmfit}} object.}
  \item{...}{other arguments}
}

\details{
When applied to a dhglmfit object, summary statistics are produced if the underlying model estimation converges. 
Otherwise, an error message is generated.
}

\value{
  It returns a list of \code{\link{summary.dhglm}} object containing the following statistics and tables.
  \item{FixCoef}{A two-dimensional matrix of fixed effect parameter estimates, respective standard errors, and t-values.}
  \item{RandCoef}{A two-dimensional matrix of random effect estimates and respective standard errors.}
  \item{iter}{Number of iterations used in \code{\link{dhglmfit}} estimation.}
  \item{convergence}{A string, "converged" if the dhglmfit estimation converged, and "did not converge" otherwise.}
  \item{mtwolikelihood}{Vectors of -2 times maximum and restricted likelihood.}
}


