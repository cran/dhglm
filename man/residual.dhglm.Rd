\name{plotdhglm}
\Rdversion{1.7}
\alias{plotdhglm}

\docType{package}

\title{Produce Model-Checking Plots for a Double Hierarchical Generalized Linear Model Analysis}

\description{Plots residuals for the mean and dispersion models}

\usage{
plotdhglm(OUTPUT, type="mean", random=NULL)
}

\arguments{
  \item{OUTPUT}{The <\code{\link{dhglmfit}}> object to be plotted}
  \item{type}{Type of model required (mean, dispersion)}
  \item{random}{Random term whose residuals are to be plotted. Default (NULL) is the residuals from the full model}
}

\details{
Four types of plot are available (normal probability plot for residuals, histogram of residuals, 
residuals versus fitted values and absolute values of residuals versus fitted values).
}

\keyword{ package }
\seealso{
<\code{\link{dhglmfit}}>
}

\examples{
#### Model checking plot for crack-growth data
data(data_crack_growth)
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle+(1|specimen),RandDist="gaussian")
res_crack_dhglm1<-dhglmfit(RespDist="gamma", DataMain=data_crack_growth, 
MeanModel=model_mu, DispersionModel=model_phi,Maxiter=1)
plotdhglm(res_crack_dhglm1)
}

