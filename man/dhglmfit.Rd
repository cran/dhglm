\name{dhglmfit}
\Rdversion{1.1}
\alias{dhglmfit}

\docType{package}

\title{Fitting Double Hierarchical Generalized Linear Models using h-likelihood Approach}

\description{dhglmfit is used to fit a class of double hierarchical generalized linear models (DHGLMs) in which
fixed and random effects can be specified in both the mean and the dispersion components. A variety of distributions and link functions for both the 
response variables and the random effects are allowed. Fixed and random effects can also be fitted in both the mean and the dispersion components.
To call the fitting function dhglmfit, models for the mean and dispersion must be specified 
by DHGLMMODELING object preferably created by calling the DHGLMMODELING function.}

\usage{
dhglmfit(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel=,
DispersionModel=,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,Maxiter=200,
convergence=1e-06)
}

\arguments{
  \item{RespDist}{The distribution of the response is set by
     the option RespDist. The user can set it to: "gaussian" (default), "binomial", "poisson", or "gamma".}
  \item{BinomialDen}{When RespDist="binomial", one should use the option BinomialDen to specify the
    denominator for the binomial distribution. This should be "NULL" (default) or a numeric vector of length equal to the length of DataMain. 
    When specified as BinomialDen=NULL and RespDist="binomial", the denominator is 1.}
  \item{DataMain}{The option DataMain determines the data frame to be used (non-optional).}
  \item{MeanModel}{For the mean model, this option requries DGHLMMODELING object which should specified by the option Model="mean".}
  \item{DispersionModel}{For the overdispersion model, this option requries DGHLMMODELING object which should be specified by the option Model="dispersion".}
  \item{PhiFix, LamFix}{Two options that determine whether the overdispersion parameters (phi) and random-effect variance (lambda) are 
     to be estimated or maintaned constant. Specifying defaults such as PhiFix =NULL (or LamFix =NULL) implies that phi (or lambda) is to be estimated. 
     If not, phi (or lambda) is fixed at a value specified by PhiFix (or LamFix).}
  \item{mord}{The option mord specifies the order of Laplace approximation to the marginal likelihood for fitting the mean parameters.
     The choice is either 0 or 1 (default).}
  \item{dord}{The option dord specifies the order of Laplace approximation to the adjusted profile likelihood for fitting the dispersion
      parameters. The choice is either 1 (default) or 2.}
  \item{Maxiter}{Maxiter determines the maximum number of iterations for estimating all parameters. The default number is 200 iterations.}
  \item{convergence}{Setting this option determines the criterion for convergence, which is computed as the absolute
      difference between the values of all the estimated parameters in the previous and current iterations. The default criterion is 1e-06.}
}

\examples{
#### Analysis of crack-growth data ####
data(data_crack_growth)
## GLM ##
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0)
model_phi<-DHGLMMODELING(Model="dispersion")
res_glm<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## JGLM ##
model_mu<-DHGLMMODELING(Model="mean", Link="log",
LinPred=y~crack0)
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle)
res_jglm<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## HGLM I ##
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion")
res_hglm1<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## HGLM II ##
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle)
res_hglm2<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## DHGLM I ##
model_mu<-DHGLMMODELING(Model="mean", Link="log",
 LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle+(1|specimen),RandDist="gaussian")
res_dhglm1<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## DHGLM II ##
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma",
LinkRandVariance="log",
LinPredRandVariance=lambda~1+(1|specimen),
RandDistRandVariance="gaussian")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle)
res_dhglm2<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
## DHGLM III ##
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~crack0+(1|specimen),RandDist="inverse-gamma",
         LinkRandVariance="log",
LinPredRandVariance=lambda~1+(1|specimen),
RandDistRandVariance="gaussian")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle+(1|specimen),RandDist="gaussian")
res_dhglm3<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)


#### Analysis of epilepsy data ####
data(data_epilepsy)
## HGLM 
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~B+T+A+T*B+V+(1|patient)+(1|id),
RandDist=c("gaussian","gamma"))
model_phi<-DHGLMMODELING(Model="dispersion")
res_hglm<-dhglmfit(RespDist="poisson",DataMain=data_epilepsy,
MeanModel=model_mu,DispersionModel=model_phi)
## DHGLM
model_mu<-DHGLMMODELING(Model="mean", Link="log", 
LinPred=y~B+T+A+T*B+V+(1|patient)+(1|id),
RandDist=c("gaussian","gaussian"),LinkRandVariance="log",
LinPredRandVariance=lambda~1+(1|id),
RandDistRandVariance="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion")
res_dhglm<-dhglmfit(RespDist="poisson",DataMain=data_epilepsy,
MeanModel=model_mu,DispersionModel=model_phi)
}

