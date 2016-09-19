\name{dhglm-package}
\Rdversion{1.6}
\alias{dhglm-package}
%\alias{dhglm}

\docType{package}

\title{Double Hierarchical Genearlized Linear Models}

\description{The dhglm package is used to fit double hierarchical generalized linear models (DHGLMs) in which random effects 
can be specified in both the mean and the dispersion components (Lee and Nelder, 2006; Lee, Nelder, and Pawitan, 2006).
It can also be used to fit generalized linear models (GLMs) of Nedler and Wedderburn (1972), joint GLMs of Nelder and Lee (1991),
and hierarchical GLMs (HGLMs) of Lee and Nelder (1996, 2001). Dispersion parameters of the random effects in the mean model can also be modeled 
with random effects (Noh, Lee and Pawitan, 2005). The response variable is allowed to follow a Gaussain, binomial, Poisson, or gamma distribution. The distribution of random
effects can be specified as Gaussian, gamma, inverse-gamma or beta. It can handle complex structures such as crossed or nested designs
in which various combinations of different distributions for random effects can be specified. Fixed effects in the mean can be estimated by maximizing
the h-likelihood or a first-order Laplace approximation to the marginal likelihood. Dispersion parameters are estimated by using first-order
adjusted profile likelihood, an extension of the restricted maximum likelihood; alternatively, 
these parameters can be assigned fixed values. The dhglm package also produces model-checking plots for various component of the model.}
\details{
\tabular{ll}{
Package: \tab dhglm\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2016-09-19\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
}
This is version 1.6 of the dhglm package.
}
\author{
Manegseok Noh, Youngjo Lee

Maintainer: Maengseok Noh <msnoh@pknu.ac.kr>
}
\references{
Lee, Y. and Nelder, J. A. (1996). Hierarchical generalised linear models (with discussion), Journal of the Royal Statistical Society B,
58, 619--678.

Lee, Y. and Nelder, J. A. (2001). Hierarchical generalised linear models : A synthesis of generalised linear
models, random-effect model and structured dispersion, Biometrika, 88, 987--1006.

Lee, Y. and Nelder, J. A. (2006). Double hierarchical generalized linear models (with discussion), Applied Statistics 55, 139--185.

Lee, Y. Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Nelder, J. A. and Lee, Y. (1991). Generalised linear models for the analysis of Taguchi-type experiments, 
Applied Stochastic Models and Data Analysis, 7, 107--120.

Nelder, J. A. and Wedderburn, R. W. M. (1972). Generalised linear models, Journal of the Royal Statistical Society A, 135, 370--384.

Noh, M., Lee, Y. and Pawitan, Y. (2005). Robust ascertainment-adjusted parameter estimation, Genetic Epidemiology, 29, 68--75.
}
\keyword{ package }
\seealso{
<\code{\link{dhglmfit}}>
}
\examples{
### DHGLM introducing random effects in the overdispersion for crack growth data
data(data_crack_growth)
model_mu<-DHGLMMODELING(Model="mean", Link="log", LinPred=y~crack0+(1|specimen),
RandDist="inverse-gamma")
model_phi<-DHGLMMODELING(Model="dispersion", Link="log", 
LinPred=phi~cycle+(1|specimen), RandDist="gaussian")
res_crack<-dhglmfit(RespDist="gamma",DataMain=data_crack_growth,
MeanModel=model_mu,DispersionModel=model_phi)
}

