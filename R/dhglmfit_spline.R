dhglmfit_spline <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=EstCorr,Dmethod=Dmethod) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL

      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
    spline=MeanModel[[16]]
    xx=matrix(0,nrow(x),p-1)
    for (i in 2:p) {
          if (i==2) xx[,1]<-x[,2]
          else {
             xx[,i-1]<-x[,i]
          }
    }
    if (p==2) colnames(xx)<-namesX[2]
    else colnames(xx)<-namesX[2:p]
    xx <- terms(nobars(formulaMean))[[3]]
    formulaMean1=formulaMean
    pp<-p-1
    print(spline)
    ppp<-0
    xx=matrix(0,nrow(x),pp)
    if (length(spline)==1) {
         xx[,1]=x[,2]
         colnames(xx)<-namesX[2]
         ppp<-1
    }
    else {
       xx<-NULL
       casecase<-1
       ppp<-0
       for (i in 1:pp) {
          if (spline[i]=="cubic" && casecase==1) {
               ppp<-ppp+1
               xx<-x[,i+1]
               casecase<-2
          } else if (casecase==2) {
               xx<-cbind(xx,x[,i+1])  
               ppp<-ppp+1
          }
       }
    }
    print(RespDist) 
    v1=NULL
    v2=NULL
    v3=NULL
    v4=NULL
    if (ppp==1) {
          if (RespDist=="gaussian") m<-lm(y~xx, data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx, family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx, family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx, family=binomial, data=DataMain)
          xx=matrix(xx,length(xx),1)
          v1=(xx[,1]-m$fitted.values)
          res1<-crPlotsHGLM(m)
    }
    if (ppp==2) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
          res1<-crPlotsHGLM(m)
    }
    if (ppp==3) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2]+xx[,3], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
          v3=(xx[,3]-m$fitted.values)
          v3=(v3-mean(v3))/sqrt(var(v3))
          res1<-crPlotsHGLM(m)
    }
    if (ppp==4) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
          v3=(xx[,4]-m$fitted.values)
          v3=(v4-mean(v4))/sqrt(var(v4))
          v4=(xx[,4]-m$fitted.values)
          v4=(v4-mean(v4))/sqrt(var(v4))
          res1<-crPlotsHGLM(m)
    }
    colnames(x)<-namesX
    ystar=y-m$fitted.values
    if (RespDist=="gaussian") res<-lm(ystar~x,data=DataMain)
    if (RespDist=="poisson") {
         ystar=round(abs(y-m$fitted.values))
         res<-glm(ystar~x, family=poisson,data=DataMain)
    }
    if (RespDist=="gamma") {
         ystar=abs(y-m$fitted.values)
         res<-glm(ystar~x, family=Gamma, data=DataMain)
    }
    if (RespDist=="binomial") res<-glm(y~x, family=binomial, data=DataMain)
    aa<-summary(res)
    beta_coeff=NULL
    if (RespDist=="gaussian") {
         beta_coeff=aa[[4]]
         print(beta_coeff,4)
    }
    else {
          beta_coeff=summary(res)[[12]]
          print(beta_coeff,4)
    }
    v_h=NULL
    if (ppp==1) v_h=v1
    if (ppp==2) v_h=cbind(v1,v2)
    if (ppp==3) v_h=cbind(v1,v2,v3)
    if (ppp==4) v_h=cbind(v1,v2,v3,v4)
    res<-list(res,res1,dstar=ystar,v_h=v_h,beta_coeff=beta_coeff)
    return(res)
}
