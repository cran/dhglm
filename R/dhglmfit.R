dhglmfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,EstCorr=TRUE,Dmethod="deviance") {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
#    date<-matrix(c(1:n),n,1)
#    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau,date))
    if (RespDist=="gaussian") MeanModel[[2]][1] <- "identity"
    if (is.null(DispersionModel[12][[1]])==FALSE) {
        if(DispersionModel[12][[1]]=="AR") {
             res<-dhglmfit_run_sv(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=0,dord=1,REML=REML,
             Maxiter=200,convergence=1e-02,Iter_mean=3)
        }
    }
    if (is.null(DispersionModel[12][[1]])==FALSE) {
    if (DispersionModel[12][[1]]=="GARCH") {
             DispersionModel[3][[1]]<-"phi~yt12"
             res<-dhglmfit_run(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="GARCH")
    }
    }
    else res<-dhglmfit_run(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
   return(res)
}
