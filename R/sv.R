sv <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1e-02,Iter_mean=3) {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    date<-matrix(c(1:n),n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau,date))
    res<-dhglmfit_run_sv(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=Iter_mean)
   return(res)
}
