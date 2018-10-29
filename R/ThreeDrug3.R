#' Drawing contour plots of the index surface
#' 
#' This function draws contour plots showing the additive/synergistic/antagonistic effect
#' @param location1 location of Drug 1's dataset on your hard disk
#' @param location2 location of Drug 2's dataset on your hard disk
#' @param location3 location of Drug 3's dataset on your hard disk
#' @param location4 location of combination dataset on your hard disk
#' @param method input the method you received in the previous step
#' @return Contour plots showing the additive/synergistic/antagonistic effect
#' @examples
#' # Try this function with the following example:
#' SYN.analysis1 <-ThreeDrug_analysis(
#' location1="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug1.txt", 
#' location2="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug2.txt", 
#' location3="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug3.txt", 
#' location4 = "C:/Users/58412/Desktop/SRA/Three drugs/data2/Combination.txt",
#' method="lilili")
#' @export


ThreeDrug_analysis <- function(location1, location2, location3, location4,method)
{
  # Purpose: to fit combination index surface of three drugs the single dose responses are known
  # Input: 
  
  #location1: location of Drug 1's dataset on your hard disk
  #Drug1: data set of Drug 1, column 1-- dose, column 2--response  
  #location2: location of Drug 2's dataset on your hard disk
  #Drug2: data set of Drug 2, column 1-- dose, column 2--response
  #location3: tion of Drug 3's dataset on your hard disk
  #Drug3: data set of Drug 3, column 1-- dose, column 2--response
  #location4: location of Combination dataset on your hard disk
  # r: number of replications at each dose-mixture, and usually 3<=r<=8
  # method: input the method you received in the previous step
  # levels: levels of contour plot, for example, c(0, 0.3, 0.5, 1, 1.5, 2)
  
  # Output:
  # Contour Plot of the interaction index surface; 
  # Red lines are 95% confidence intervals of the additive action (combination index=1)  
  # synergistic, antagonistic and additive multifold
  
  ## this programm fit the three drugs combination data by using additive model and B-spline approximation by using
  ## m times combination and r times replications for r times
  
  #' @import Matrix
  #' @import splines
  #' @import scatterplot3d
  #library(quantreg)
  #library(ucminf)
  set.seed(10^7) 
  Drug1 <-read.table(file= location1, header = T)
  Drug2 <-read.table(file= location2, header = T)
  Drug3 <-read.table(file= location3, header = T)
  Combination <- read.table(file= location4, header = T)
  #id0 <-read.table("Drug3.txt", header = T)
  r<-8
  levels <- c(0, 0.3, 0.5, 1, 1.5, 2)
  level1 <- c(0, 10, 20, 30, 40, 50)
  tau=0.5
  Comb <- Combination
  if ( method != "logloglog" & method != "lilili" & method != "loglili" & method != "lililog" & method !="lilogli" & method !="loglogli" & method !="loglilog" & method!= "liloglog") 
  { 
    stop( "method = unknown regression model." ) 
  }
  
  ##--- estimation of the single dose-responses
  if ( method == "logloglog") {
    coef1 <- as.numeric(lm(Drug1[,2] ~ log(Drug1[,1]))$coef)      #intercept and slope of drug 1 
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)      #intercept and slope of drug 2 	
    coef3 <- as.numeric(lm(Drug3[,2] ~ log(Drug3[,1]))$coef)
    SD1 <- exp((Comb[, 4] - coef1[1])/coef1[2])
    SD2 <- exp((Comb[, 4] - coef2[1])/coef2[2])
    SD3 <- exp((Comb[, 4] - coef3[1])/coef3[2])
  }
  if ( method == "lilili") {
    coef1 <- as.numeric((lm(Drug1[,2] ~ Drug1[,1]))$coef)
    coef2 <- as.numeric((lm(Drug2[,2] ~ Drug2[,1]))$coef)
    coef3 <- as.numeric((lm(Drug3[,2] ~ Drug3[,1]))$coef)
    SD1 <- ((Comb[, 4] - coef1[1])/coef1[2])
    SD2 <- ((Comb[, 4] - coef2[1])/coef2[2])
    SD3 <- ((Comb[, 4] - coef3[1])/coef3[2])
  }
  
  if ( method == "lilogli") {
    coef1 <- as.numeric((lm(Drug1[,2] ~ Drug1[,1]))$coef)
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)
    coef3 <- as.numeric(lm(Drug3[,2] ~ (Drug3[,1]))$coef)
    SD1 <- ((Comb[, 4] - coef1[1])/coef1[2])
    SD2 <- exp((Comb[, 4] - coef2[1])/coef2[2])
    SD3 <- ((Comb[, 4] - coef3[1])/coef3[2])
  }
  
  if ( method == "loglogli"){
    coef1 <- as.numeric(lm(Drug1[,2] ~ log(Drug1[,1]))$coef)     
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)
    coef3 <- as.numeric((lm(Drug3[,2] ~ Drug3[,1]))$coef)
    SD1 <- exp((Comb[, 4] - coef1[1])/coef1[2])
    SD2 <- exp((Comb[, 4] - coef2[1])/coef2[2])
    SD3 <- ((Comb[, 4] - coef3[1])/coef3[2])
  }
  
  if ( method == "loglili" ) {
    method <- "lilogli"
    Drug1 <- Drug2
    Drug2 <- Drug1
  }
  
  if ( method == "lililog" ) {
    method <- "lilogli"
    Drug.temp <- Drug2
    Drug2 <- Drug3
    Drug3 <- Drug.temp
  }  
  
  if ( method == "loglilog"){
    method <- "loglogli"
    Drug.temp <- Drug2
    Drug2 <- Drug3
    Drug3 <- Drug.temp
  }
  
  if ( method == "liloglog"){
    method <- "loglogli"
    Drug.temp <- Drug1
    Drug1 <- Drug3
    Drug3 <- Drug.temp 
  }
  
  
  Index <- Comb[, 1]/SD1 + Comb[, 2]/SD2+Comb[, 3]/SD3
  WW <- cbind(Comb, Index)
  
  
  n.Comb <- dim(Comb)[1]
  m <- n.Comb/r
  
  
  ###-------------------function of select optimal kn-----------------------
  opkn <- function(y, time1,time2,time3, N, degree)
  {
    
    Okn=NULL
    for (i in 1:r)
    {
      
      
      SIC = NULL
      for(kn in 0:3)       
      {
        q = kn+degree
        u = seq(0, 1, length=kn+2)[-c(1,kn+2)]
        Knots1 = as.numeric(quantile(time1[,i], u))
        Knots2 = as.numeric(quantile(time2[,i], u))
        Knots3 = as.numeric(quantile(time3[,i], u))
        pi.t1 = bs(time1[,i], knots=Knots1, intercept=TRUE, degree=degree)
        pi.t2 = bs(time2[,i], knots=Knots2, intercept=TRUE, degree=degree)
        pi.t3 = bs(time3[,i], knots=Knots3, intercept=TRUE, degree=degree)
        G = cbind(pi.t1,pi.t2, pi.t3)
        pkn=ncol(G)
        lse=(lm(y[,i] ~ G))
        res=lse$resid
        SIC = c(SIC, log(sum(res^2)) + 0.5*log(N)*pkn/N) 
      }
      okn = max(which(SIC == min(SIC))-1,0)
      Okn=c(Okn,okn)
    }
    ################# Fit the model again with the optimal kn
    fiokn=floor(mean(Okn))
    return(list(fiokn))
  }
  ##-----------function of compute estimator--------------------------------------------------
  
  estpar<-function(Y, X1,X2,X3, okn, tau, degree)
  {
    Ecoef=matrix(0,r,(3*(okn+degree)+1))##estimated coefficient
    Est=matrix(0,r,m)##estimated index
    Estd=matrix(0,r,(3*(okn+degree)+1))
    for(i in 1:r){
      
      
      q = okn+degree
      u = seq(0, 1, length=okn+2)[-c(1,okn+2)]
      Knots1 = as.numeric(quantile(X1[, i], u))
      Knots2 = as.numeric(quantile(X2[, i], u))
      Knots3 = as.numeric(quantile(X3[, i], u))
      pi.t1 = bs(X1[, i], knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
      pi.t2 = bs(X2[, i], knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
      pi.t3 = bs(X3[, i], knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
      G = cbind(pi.t1,pi.t2, pi.t3)
      dm=matrix(cbind(rep(1,m),G),ncol=q*3+1,byrow=F)
      #est=rq(Y[,i]~G, tau=tau)$coeff#quantile regression
      est= lm(Y[,i]~ G)$coef#LS---slightly better than qr
      Idi=dm%*%matrix(est)
      Est[i,]=Idi#estimated index
      Ecoef[i ,]=est# estimated parameter
      
    }
    return(list(eindex=Est,ecoeff=Ecoef,sd=Estd))
  }
  ###--------------------------------------------------------------------------------
  
  ##------------------variance estimation---------------------------------------
  varpar<-function(Y,X1,X2,X3,hatpar,kn,tau,degree,m)
  {
    X10=X1[, 1]
    X20=X2[, 1]
    X30=X3[, 1]
    q = kn+degree
    u = seq(0, 1, length=kn+2)[-c(1,kn+2)]
    Knots1 = as.numeric(quantile(X10, u))
    Knots2 = as.numeric(quantile(X20, u))
    Knots3 = as.numeric(quantile(X30, u))
    pi.t1 = bs(X10, knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
    pi.t2 = bs(X20, knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
    pi.t3 = bs(X30, knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
    G = cbind(pi.t1,pi.t2, pi.t3)
    dm=matrix(cbind(rep(1,m),G),ncol=q*3+1,byrow=F)# design matrix
    Y0=matrix(0,m,1)
    for(i in 1:m)
    {Y0[i,]=mean(Y[i,])}
    res=Y0-dm%*%hatpar
    sige=sum(res^2)-(mean(res))^2
    sigma=solve(t(dm)%*%dm)*sige/(m-(3*q+1))
    return(sigma)}
  ##----------------------------------------------------------------------------------
  XX1=matrix(rep(WW[,1],1),ncol=r,byrow=T)
  XX2=matrix(rep(WW[,2],1),ncol=r,byrow=T)
  XX3=matrix(rep(WW[,3],1),ncol=r,byrow=T)
  YY=matrix(rep(WW[,5],1),ncol=r,byrow=T)
  YYY=matrix(rep((WW[,4])^(1/10),1),ncol=r,byrow=T)
  
  # YYY=matrix(rep(WW[,4],1),ncol=r,byrow=T)
  
  #okn=0#number of interior knots for interaction index
  #okn1=0
  degree=3
  NB=2
  okn <- as.numeric(opkn(YY, XX1,XX2,XX3,m,degree))
  okn1 <- as.numeric(opkn(YYY, XX1,XX2,XX3,m,degree))
  hatY<-hatid<-matrix(0,m,1)
  estvar=matrix(0,3*(okn+degree)+1,3*(okn+degree)+1)
  Par<-matrix(0,3*(okn+degree)+1,NB)
  ParY<-matrix(0,3*(okn1+degree)+1,NB)
  estvarY=matrix(0,3*(okn1+degree)+1,3*(okn1+degree)+1)
  
  for (j in 1:NB){
    
    YYY0<-YY0<-matrix(0,m,r)
    bootf=matrix(0,m,r)
    for (i in 1:m){bfoot=sample(r, r, replace = F)
    YYi=YY[i,]
    YY0[i,]=YYi[bfoot]
    YYYi=YYY[i,]
    YYY0[i,]=YYYi[bfoot]
    bootf[i,]=bfoot}
    
    
    ####----------estimate the parameter for interaction index----------------------
    
    
    
    Est=estpar(YY0, XX1,XX2,XX3, okn, tau, degree)$eindex
    Ecoef=estpar(YY0, XX1,XX2,XX3, okn, tau, degree)$ecoeff
    
    ####----------estimate the parameter for response-----------------------
    
    Est1=estpar(YYY0, XX1,XX2,XX3, okn1, tau, degree)$eindex
    Ecoef1=estpar(YYY0, XX1,XX2,XX3, okn1, tau, degree)$ecoeff
    
    estid=matrix(rep(apply(Est,2,mean),1),ncol=1,byrow=T)
    esty=matrix(rep(apply(Est1,2,mean),1),ncol=1,byrow=T)
    hatid=hatid+estid/NB
    hatY=hatY+esty/NB
    
    pari=matrix(rep(apply(Ecoef,2,mean),1),ncol=1,byrow=T)
    Par[,j]=pari
    ParY[,j]=matrix(rep(apply(Ecoef1,2,mean),1),ncol=1,byrow=T)
    
    estvari=varpar(YY0,XX1,XX2,XX3,pari,okn,tau,degree,m)
    estvar=estvar+estvari/NB	
    
    estvari1=varpar(YYY0,XX1,XX2,XX3,pari,okn,tau,degree,m)
    estvarY=estvarY+estvari1/NB
    
  }##end for bootstrap
  
  W0=matrix(rep(Index,1),ncol=r,byrow=T)
  id=matrix(rep(apply(W0,1,mean),1),ncol=1,byrow=T)
  hatid-id
  y=matrix(rep(apply(YYY0,1,mean),1),ncol=1,byrow=T)
  (hatY)^10-y
  
  
  
  par=apply(Par,1,mean)
  pary=apply(ParY,1,mean)
  
  sdpar=diag(estvar)^(1/2)
  
  sdpary=diag(estvarY)^(1/2)
  
  #######---Contour plot for fixed x1(dose of Drug1)--------------------------------
  opar <- par(no.readonly = TRUE)
  par(mar = rep(4.5,4)) 
  n.grid <- 50
  
  x1 <- seq(0, max(WW[, 1]), length = n.grid)
  x2 <- seq(min(WW[, 2]), max(WW[, 2]), length = n.grid)
  x3 <- seq(min(WW[, 3]), max(WW[, 3]), length = n.grid)
  
  
  W <- W1 <- W2 <-matrix(0, n.grid, n.grid)
  x23<-as.matrix(expand.grid(x2, x3))
  XX10=XX1[, 1]
  XX20=XX2[, 1]
  XX30=XX3[, 1]
  q = okn+degree
  u = seq(0, 1, length=okn+2)[-c(1,okn+2)]
  Knots1 = as.numeric(quantile(XX10, u))
  Knots2 = as.numeric(quantile(XX20, u))
  Knots3 = as.numeric(quantile(XX30, u))
  pi.t1 = bs(x1, knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t2 = bs(x23[,1], knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t3 = bs(x23[,2], knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
  i=1##--dose of Drug1 is 0
  pit1fi=pi.t1[i,]
  px1=matrix(rep(pit1fi,n.grid^2),ncol=q,byrow=T)
  dm=matrix(cbind(rep(1,n.grid^2),px1,pi.t2, pi.t3),ncol=q*3+1,byrow=F)# design matrix
  Wf1=matrix(dm %*% par,ncol=n.grid,byrow=F)
  
  
  
  ######-----------------------confidence interval--------------------------------------------
  
  std_f1=matrix((diag(dm%*%estvar%*%t(dm)))^(1/2),ncol=n.grid,byrow=F)
  W11=Wf1+qnorm(0.975, 0, 1)*std_f1
  W12=Wf1-qnorm(0.975, 0, 1)*std_f1
  
  WWf1=dm %*% par
  sstd_f1=(diag(dm%*%estvar%*%t(dm)))^(1/2)
  
  lwf1=(1-qnorm(0.975, 0, 1)*sstd_f1) 
  sdose1=x23[WWf1<lwf1 ,]
  write.table(sdose1, file="syn23.txt") 
  
  
  write.table(matrix(cbind(W,W1,W2),ncol=3,byrow=F), file="index-plot.txt")
  
  cont.lineposi1 <- contourLines(x2, x3, W11)
  cont.lineneg1 <- contourLines(x2, x3, W12)
  cont.linecent1 <- contourLines(x2, x3, Wf1)
  
  
  
  color_temp.palette <- palette( c("#00FF00", "#00A000", "yellow", "#FF0000", "#800000" ))
  filled.contour(x2, x3, Wf1, levels = levels,
                 col = color_temp.palette, xlim = range(x2), ylim = range(x3),
                 plot.title = title(main = "Contour plot of Drug2 and 3 (fixed Drug1=0)", xlab = "Drug 2 Dose", ylab = "Drug 3 Dose"))
  
  #Plotting the confidence intervals
  mar.orig <- par("mar")
  w <- (3 + mar.orig[2]) * par("csi") * 2.54
  layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
  
  ind <- 0
  for (i in 1:length(cont.lineposi1)) {
    if ( abs(cont.lineposi1[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineposi1[ind][[1]]$x, cont.lineposi1[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  ind <- 0
  for (i in 1:length(cont.lineneg1)) {
    if ( abs(cont.lineneg1[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineneg1[ind][[1]]$x, cont.lineneg1[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  
  
  #####---------------------------------------------------------------------------------------
  #######---Contour plot for fixed x2(dose of Drug2)--------------------------------
  
  n.grid <- 50
  
  x1 <- seq(min(WW[, 1]), max(WW[, 1]), length = n.grid)
  x2 <- seq(0, max(WW[, 2]), length = n.grid)
  x3 <- seq(min(WW[, 3]), max(WW[, 3]), length = n.grid)
  
  
  W <- W1 <- W2 <-matrix(0, n.grid, n.grid)
  x13<-as.matrix(expand.grid(x1, x3))
  XX10=XX1[, 1]
  XX20=XX2[, 1]
  XX30=XX3[, 1]
  q = okn+degree
  u = seq(0, 1, length=okn+2)[-c(1,okn+2)]
  Knots1 = as.numeric(quantile(XX10, u))
  Knots2 = as.numeric(quantile(XX20, u))
  Knots3 = as.numeric(quantile(XX30, u))
  pi.t1 = bs(x13[,1], knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t2 = bs(x2, knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t3 = bs(x13[,2], knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
  i=1
  pit2fi=pi.t2[i,]
  px2=matrix(rep(pit2fi,n.grid^2),ncol=q,byrow=T)
  dm=matrix(cbind(rep(1,n.grid^2),pi.t1,px2, pi.t3),ncol=q*3+1,byrow=F)# design matrix
  Wf2=matrix(dm %*% par,ncol=n.grid,byrow=F)
  
  ###-----------------------confidence interval--------------------------------------------
  #        W1=matrix(dm %*% parhi,ncol=n.grid,byrow=F)
  #        W2=matrix(dm %*% parlw,ncol=n.grid,byrow=F)
  std_f2=matrix((diag(dm%*%estvar%*%t(dm)))^(1/2),ncol=n.grid,byrow=F)
  W1=Wf2+qnorm(0.975, 0, 1)*std_f2
  W2=Wf2-qnorm(0.975, 0, 1)*std_f2
  WWf2=dm %*% par
  sstd_f2=(diag(dm%*%estvar%*%t(dm)))^(1/2)
  lwf2=(1-qnorm(0.975, 0, 1)*sstd_f2) 
  sdose2=x13[WWf2<lwf2 ,]
  
  write.table(sdose2, file="syn13.txt")      
  
  
  
  
  write.table(matrix(cbind(W,W1,W2),ncol=3,byrow=F), file="index-plot.txt")
  
  cont.lineposi <- contourLines(x1, x3, W1)
  cont.lineneg <- contourLines(x1, x3, W2)
  cont.linecent <- contourLines(x1, x3, Wf2)
  
  
  
  color_temp.palette <- palette( c("#00FF00", "#00A000", "yellow", "#FF0000", "#800000" ))
  filled.contour(x1, x3, Wf2, levels = levels, 
                 col = color_temp.palette, xlim = range(x1), ylim = range(x3),
                 plot.title = title(main = "Contour plot of Drug1 and 3 (fixed Drug2=0)", xlab = "Drug 1 Dose", ylab = "Drug 3 Dose"))
  
  #Plotting the confidence intervals
  mar.orig <- par("mar")
  w <- (3 + mar.orig[2]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
  
  ind <- 0
  for (i in 1:length(cont.lineposi)) {
    if ( abs(cont.lineposi[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineposi[ind][[1]]$x, cont.lineposi[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  ind <- 0
  for (i in 1:length(cont.lineneg)) {
    if ( abs(cont.lineneg[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineneg[ind][[1]]$x, cont.lineneg[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  
  
  
  #######---Contour plot for fixed x3(dose of Drug3)--------------------------------
  
  n.grid <- 50
  
  x1 <- seq(min(WW[, 1]), max(WW[, 1]), length = n.grid)
  x2 <- seq(min(WW[, 2]), max(WW[, 2]), length = n.grid)
  x3 <- seq(0, max(WW[, 3]), length = n.grid)
  
  W <- W1 <- W2 <-matrix(0, n.grid, n.grid)
  x12<-as.matrix(expand.grid(x1, x2))
  XX10=XX1[, 1]
  XX20=XX2[, 1]
  XX30=XX3[, 1]
  q = okn+degree
  u = seq(0, 1, length=okn+2)[-c(1,okn+2)]
  Knots1 = as.numeric(quantile(XX10, u))
  Knots2 = as.numeric(quantile(XX20, u))
  Knots3 = as.numeric(quantile(XX30, u))
  pi.t1 = bs(x12[,1], knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t2 = bs(x12[,2], knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t3= bs(x3, knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
  
  i=1
  pit3fi=pi.t3[i,]
  px3=matrix(rep(pit3fi,n.grid^2),ncol=q,byrow=T)
  dm=matrix(cbind(rep(1,n.grid^2),pi.t1,pi.t2, px3),ncol=q*3+1,byrow=F)# design matrix
  Wf3=matrix(dm %*% par,ncol=n.grid,byrow=F)
  
  ###-----------------------confidence interval--------------------------------------------
  #        W1=matrix(dm %*% parhi,ncol=n.grid,byrow=F)
  #        W2=matrix(dm %*% parlw,ncol=n.grid,byrow=F)
  std_f3=matrix((diag(dm%*%estvar%*%t(dm)))^(1/2),ncol=n.grid,byrow=F)
  W1=Wf3+qnorm(0.975, 0, 1)*std_f3
  W2=Wf3-qnorm(0.975, 0, 1)*std_f3
  
  WWf3=dm %*% par
  sstd_f3=(diag(dm%*%estvar%*%t(dm)))^(1/2)
  
  lwf1=(1-qnorm(0.975, 0, 1)*sstd_f3) 
  sdose1=x12[WWf3<lwf1 ,]
  
  
  write.table(sdose1, file="syn12.txt")
  
  
  
  write.table(matrix(cbind(W,W1,W2),ncol=3,byrow=F), file="index-plot.txt")
  
  cont.lineposi <- contourLines(x1, x2, W1)
  cont.lineneg <- contourLines(x1, x2, W2)
  cont.linecent <- contourLines(x1, x2, Wf3)
  
  
  
  color_temp.palette <- palette( c("#00FF00", "#00A000", "yellow", "#FF0000", "#800000" ))
  filled.contour(x1, x2, Wf3, levels = levels, 
                 col = color_temp.palette, xlim = range(x1), ylim = range(x2), 
                 plot.title = title(main = "Contour plot of Drug1 and 2 (fixed Drug3=0)", xlab = "Drug 1 Dose", ylab = "Drug 2 Dose"))
  
  #Plotting the confidence intervals
  mar.orig <- par("mar")
  w <- (3 + mar.orig[2]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
  
  ind <- 0
  for (i in 1:length(cont.lineposi)) {
    if ( abs(cont.lineposi[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineposi[ind][[1]]$x, cont.lineposi[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  ind <- 0
  for (i in 1:length(cont.lineneg)) {
    if ( abs(cont.lineneg[i][[1]]$level - 1) < 1e-5)	{
      ind <- i
      break;
    }
  }
  if (ind != 0)	{
    lines(cont.lineneg[ind][[1]]$x, cont.lineneg[ind][[1]]$y, type = "l", lty = 2, lwd = 1, col=rgb(0,0,0))	
  }
  
  par(opar)
  
  
  
  
  
  ###-------------------------------------------------------------------------------------------------------------------------------------------
  ##------------------------find the synegestic surface----------------------------------
  n.grid=50
  nc=100
  xx1 <- seq(min(WW[,1]), max(WW[, 1]), length = n.grid)
  x2 <- seq(min(WW[, 2]), max(WW[, 2]), length = n.grid)
  x3 <- seq(min(WW[, 3]), max(WW[, 3]), length = n.grid)
  
  
  
  x123<-as.matrix(expand.grid(xx1,x2, x3))
  Knots1 = as.numeric(quantile(XX10, u))
  Knots2 = as.numeric(quantile(XX20, u))
  Knots3 = as.numeric(quantile(XX30, u))
  pi.t4 = bs(x123[,1], knots=Knots1, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t5 = bs(x123[,2], knots=Knots2, intercept=TRUE, degree=degree)[,1:(q)]
  pi.t6 = bs(x123[,3], knots=Knots3, intercept=TRUE, degree=degree)[,1:(q)]
  dm1=matrix(cbind(rep(1,dim(pi.t4)[1]),pi.t4,pi.t5, pi.t6),ncol=q*3+1,byrow=F)# design matrix
  W3=dm1 %*% par
  W4=(dm1 %*% pary)
  ##----------------standard deviation at thest point--------------------
  mm=dim(matrix(W3))[1]
  std2=std1=matrix(0,mm,1)
  
  mm1=mm/nc
  for(i in 0:(nc-1)){
    
    n1=mm1*i+1
    n2=mm1*(i+1)
    dmi=dm1[n1:n2,]
    stdi=matrix((diag(dmi%*%estvar%*%t(dmi)))^(1/2),ncol=1,byrow=T)
    std1[n1:n2,]=stdi
    stdi2=matrix((diag(dmi%*%estvarY%*%t(dmi)))^(1/2),ncol=1,byrow=T)
    std2[n1:n2,]=stdi2
  }
  
  
  ###---------synergistic multifold-------------------------
  lw=(1-qnorm(0.975, 0, 1)*std1) 
  uw=(1+qnorm(0.975, 0, 1)*std1)    
  xn=x123[W3>lw ,]
  
  
  #color <- rep("yellow",length( xn[,2]))
  #scatterplot3d(xn[,2],xn[,3],xn[,1], color, 
  #col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
  #cex.lab = 1.1, main = "synergistic multifold", pch = 20, 
  #xlab="Drug 1", ylab="Drug 2", zlab="Drug 3",)
  
  par(mfrow=c(1,1))
  
  scatterplot3d(xn[,2],xn[,3],xn[,1], highlight.3d = T, 
                col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                cex.lab = 1.1, main = "synergistic multifold", pch = 20, 
                xlab="Drug 1", ylab="Drug 2", zlab="Drug 3")
  
  write.table(xn, file="synergistic.txt")
  par(opar)
  
  
  
  
  
  
  ###-----------antagonistic multifold-----------------
  xn1=x123[W3>uw ,]
  par(mfrow=c(1,1))
  
  scatterplot3d(xn1[,2],xn1[,3],xn1[,1], highlight.3d = TRUE, 
                col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                cex.lab = 1.1, main = "antagonistic multifold", pch = 20, 
                xlab="Drug 2", ylab="Drug 3", zlab="Drug 1")
  
  write.table(xn1, file="antagonistic.txt")
  par(opar)
  
  
  
  ###-----------additive multifold-----------------
  xn2=x123[(W3<=uw)&(W3>=lw) ,]
  par(mfrow=c(1,1))
  
  scatterplot3d(xn2[,2],xn2[,3],xn2[,1], highlight.3d = TRUE, 
                col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                cex.lab = 1.1, main = "additive multifold", pch = 20, 
                xlab="Drug 2", ylab="Drug 3", zlab="Drug 1")
  
  write.table(xn2, file="addtive.txt")
  par(opar)
  
  
  ###########################################################################
  
  ###----------IC50 of response------------------------------
  # yic50=50^{1/10}
  lw2=(50-3) 
  uw2=(50+3)
  xn4=x123[(W4^10<=uw2)&(W4^10>=lw2) ,]
  xn5=x123[W4^10<lw2 ,]
  write.table(xn4, file="ic50.txt")
  
  
  #scatterplot3d(xn4[,2],xn4[,3],xn4[,1], highlight.3d = TRUE, 
  #col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
  #cex.lab = 1.1, main = "Dose region with viability smaller than 50%", pch = 20, 
  #xlab="Drug 1", ylab="Drug 2", zlab="Drug 3",)
  
  par(mfrow=c(1,2))
  
  scatterplot3d(xn4[,2],xn4[,3],xn4[,1], highlight.3d = TRUE, 
                col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                cex.lab = 1.1, main = "IC50 of the combination", pch = 20, 
                xlab="Drug 2", ylab="Drug 3", zlab="Drug 1")
  
  scatterplot3d(xn5[,2],xn5[,3],xn5[,1], highlight.3d = TRUE, 
                col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                cex.lab = 1.1, main = "Dose region with viability smaller than IC50", pch = 20, 
                xlab="Drug2", ylab="Drug3", zlab="Drug1")
  
  par(opar)
  
}


 