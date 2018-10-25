DHGLMMODELING <-
function(Model="mean",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL,LMatrix=NULL,LinkRandVariance=NULL,LinPredRandVariance=NULL,
RandDistRandVariance="gaussian",LinkRandVariance2=NULL,LinPredRandVariance2=NULL,corr=NULL,spatial=NULL,longitude=NULL,latitude=NULL,spline=NULL,Neighbor=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    if (is.null(corr)) corr="IND"
    if (is.null(spatial)) spatial="IND"
    if (is.null(RandDist)) RandDist="gaussian"
    res<-list(Model,Link,LinPred,RandDist,Offset,LMatrix,LinkRandVariance,LinPredRandVariance,
              RandDistRandVariance,LinkRandVariance2,LinPredRandVariance2,corr,spatial,longitude,latitude,spline,Neighbor)
    return(res)
}
