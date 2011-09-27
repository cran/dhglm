dhglmfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=1,dord=1,Maxiter=200,convergence=1e-06) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrast=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
      z<-zz<-z[[1]]
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
    x_lambda <- fr_lambda$X
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
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
    DispersionModel_1<-DispersionModel[3][[1]]
    if (DispersionModel[3][[1]]=="constant") DispersionModel[3][[1]]<-phi~1
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrast=NULL)
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
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
##    print(model_number)
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if (RespDist=="poisson" || RespDist=="binomail") PhiFix<-1
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=inv_disp,offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       oq_disp<-matrix(1,qcum_disp[nrand+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
    if (nrand>1) {
        II<-diag(rep(1,n))
        z<-cbind(z,II)
    }
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
       if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=inv_disp,offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-5
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    beta_h<-beta_mu
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%(y-mu)-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="gamma") {
            dhdv<-t(z)%*%(y-mu)+1/lambda-exp(v_h)/lambda
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda)+exp(-v_h)/lambda
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    v_h_old<-v_h
    v_h<-v_h-solve(d2hdv2)%*%dhdv
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
    if (RandDist=="gaussian") u_h <- v_h
    if (RandDist=="gamma") u_h <-exp(v_h)
    if (RandDist=="inverse-gamma") u_h <-exp(v_h)
    if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    a<-matrix(0,n,1)
    if (q[1]>0 && mord==1 && RespDist!="gamma") {
    T<-t(cbind(t(z),I))
##    Null<-matrix(0,n,n)
##    W<-matrix(0,(2*n),(2*n))
##    W[c(1:n),]<-cbind(W1,Null)
##    W[c((n+1):(2*n)),]<-cbind(Null,W2)   
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--z%*%solve(t(T)%*%W%*%T)%*%t(z)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qq in 1:n){
            d2[i]<-d2[i]+P[qq,qq]*K1[qq,i]
        }
        if (RandDist=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+z%*%solve(W2)%*%t(z))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    beta_h<-solve(t(x)%*%invSig%*%x)%*%(t(x)%*%invSig%*%(z1-a))
    se_beta<-sqrt(diag(solve(t(x)%*%invSig%*%x)))
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
       TT<-rbind(cbind(x,z),cbind(OO1,I))
       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW[c(1:n),]<-cbind(W1,Null1)
       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
       PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
       leverage<-rep(0,n)
       for (kk in 1:n) leverage[kk]<-PP[kk]
    }
    resp_disp<-deviance_residual/(1-leverage)
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2
##############################################################
######### GLM fit for phi #####################
##############################################################
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma") {
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weight=weight_disp,offset=Offset_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est<-1/inv_disp
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    if (q_disp[1]>0) {
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(resp_disp,x_disp,disp_rand)
       reshglm_disp<-hglmfit_corr(resp_disp~x_disp-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=1)
       disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
              temp17<-2*(1-log(u_h)-(1-u_h))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    }
##############################################################
######### GLM fit for lambda            #####################
##############################################################
  if (is.null(LamFix)) {
    if (q[1]>0 && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
##       print(resp_lambda)
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
   model_number1<-0
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && q_lambda[1]>0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
###       lambda<-lambda_est1
       lambda_est<-lambda
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_err<-convergence3
    names(print_i) <- "iteration : "
    print(print_i)
    names(print_err) <- "convergence : "
    print(print_err)
    max_iter<-max_iter+1
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    md<-RespDist
    names(md)<-"Distribution of Main Response : "
    print(md)
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        z_beta<-beta_h/se_beta
        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (q[1]>0 && q_lambda[1]==0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        res3<-summary(resglm_lambda,dispersion=2)
        p_lambda<-nrand
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]*1.1
#        lambda_h[2]<-lambda_h[2]*1.7
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]*sqrt(1.1)
#        lambda_se[2]<-lambda_se[2]*sqrt(1.7)
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(lambda_h,lambda_se)
        colnames(lambda_coeff) <- c("Estimate", "Std. Error")
        rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }        
    if (q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        res3<-summary(resglm_lambda,dispersion=2)
        if (nrand==1) {
           lambda_coeff<-cbind(beta_lambda,se_lambda)
           colnames(lambda_coeff) <- c("Estimate", "Std. Error")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(lambda_h,lambda_se)
           colnames(lambda_coeff) <- c("Estimate", "Std. Error")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        alpha_coeff<-cbind(beta_alpha,se_alpha[1,1])
        colnames(alpha_coeff) <- c("Estimate", "Std. Error")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           phi_coeff<-cbind(beta_phi,se_phi)
           colnames(phi_coeff) <- c("Estimate", "Std. Error")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           phi_coeff<-cbind(beta_phi,se_phi)
           colnames(phi_coeff) <- c("Estimate", "Std. Error")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           tau_coeff<-cbind(beta_tau,se_tau[1,1])
           colnames(tau_coeff) <- c("Estimate", "Std. Error")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%W1%*%x
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
    }
    if (model_number >=3) {
       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%W2%*%v_h-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))  
       if (RandDist=="gamma") hv<-log(u_h-u_h)/lambda_est-log(lambda_est)/lambda_est-lgamma(1/lambda_est)
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) hv12<-reshglm_disp[[8]][3]
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       AA<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(-1*d2hdv2)))
       BB<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(t(z)%*%W1%*%z)))
       rl<--2*hlikeli-2*sum(hv)+hv10+hv20+log(abs(det(AA/(2*pi))))
       pd<- sum(diag(solve(AA) %*% BB))
       caic<- -2*hlikeli + 2*pd
    }
    likeli_coeff<-rbind(ml,rl,caic)
    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ")
    }
    if (model_number<4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ")
    }
    print(likeli_coeff)
    res<-list(mean_residual,mu)
    return(res)
}

