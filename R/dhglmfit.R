dhglmfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,EstCorr=FALSE,Dmethod="deviance") {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
#    date<-matrix(c(1:n),n,1)
#    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau,date))
    if (RespDist=="gaussian" && is.null(MeanModel[[2]][1])) MeanModel[[2]][1] <- "identity"
  if (is.null(MeanModel[[6]][1])==TRUE) {
    if (is.null(DispersionModel[[12]][1])==FALSE) {
        if(DispersionModel[[12]][1]=="AR") {
             MeanModel[[2]][1] <- "identity"
             res<-dhglmfit_run_sv(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=0,dord=1,REML=REML,
             Maxiter=200,convergence=1e-02,Iter_mean=3)
        }
    }
    if (is.null(DispersionModel[[12]][1])==FALSE) {
    if (DispersionModel[[12]][1]=="GARCH") {
             MeanModel[[2]][1] <- "identity"
             DispersionModel[3][[1]]<-"phi~yt12"
             res<-dhglmfit_run_GARCH(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="GARCH")
    }
    if(is.null(MeanModel[[16]][1])==TRUE) MeanModel[[16]][1]=="temp"
    if (DispersionModel[[12]][1]=="IND" && MeanModel[[13]][1]=="IND" && is.null(MeanModel[[16]][1])==TRUE) {
            res<-dhglmfit_run(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
    }
    }
    if (MeanModel[[13]][1]=="Matern" && is.null(MeanModel[[14]][1])==FALSE && is.null(MeanModel[[15]][1])==FALSE) 
    if (MeanModel[[13]][1]=="MRF" || MeanModel[[13]][1]=="IAR" && is.null(MeanModel[[17]][1])==FALSE) {
         res<-HGLMREML.h(formulaMain=MeanModel[[3]],DataMain=DataMain,Offset=MeanModel[[5]],RespDist=RespDist,
             RespLink=MeanModel[[2]],RandDist="normal", spatial=MeanModel[[13]],Neighbor=MeanModel[[17]]) 
    }
    if (is.null(MeanModel[[16]][1])==FALSE && is.null(DispersionModel[[16]][1])==TRUE) res<-dhglmfit_spline(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
#    if (is.null(MeanModel[[16]][1])==FALSE && is.null(DispersionModel[[16]][1])==FALSE) res<-dhglmfit_spline2(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
#             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
#             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
  } else {
             res<-dhglmfit_run_Lmatrix(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5)
   }
   return(res)
}
