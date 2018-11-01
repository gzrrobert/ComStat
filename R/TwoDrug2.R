#' Testing the additive effect and drawing a contour plot
#' 
#' This function implements an F test to find out the additive/synergistic/antagonistic effect under each two
#' drug combination level, and draws a contour plot showing this effect
#' @param location1 location of Drug 1's dataset on your hard disk
#' @param location2 location of Drug 2's dataset on your hard disk
#' @param location4 location of combination dataset on your hard disk
#' @param r number of replications at each dose-mixture, and usually 3<=r<=8
#' @param method input the method you received in the previous step
#' @return F test of additive/synergistic/antagonistic effect under drug combination,
#'         a contour plot showing this effect
#' @examples
#' # Try this function with the following example:
#' TwoDrug_Test_Analysis_out <- TwoDrug_Test_Analysis(
#'  location1="C:/Users/58412/Desktop/SRA/Two drugs/Drug1.txt",
#'  location2="C:/Users/58412/Desktop/SRA/Two drugs/Drug2.txt",
#'  location4 = "C:/Users/58412/Desktop/SRA/Two drugs/Combination.txt", 
#'                                               r=6,method="logli")
#' @export

#' @import fields


TwoDrug_Test_Analysis <- function(location1,location2,location4, r=5,method="loglog"){
  
  #location1: location of Drug 1's dataset on your hard disk
  #Drug1: data set of Drug 1, column 1-- dose, column 2--response  
  #location2: location of Drug 2's dataset on your hard disk
  #Drug2: data set of Drug 2, column 1-- dose, column 2--response
  #location4: location of Combination dataset on your hard disk
  # c0: the smallest meaningful difference to be detected
  # r: number of replications at each dose-mixture, and usually 3<=r<=5
  # method: input the method you received in the previous step
  # All dataset should be in txt file
  
  #1st Drug
  Drug1<-read.table(file = location1, header=T)
  #if response value is larger than 100%, convert it to 100%
  for (i in 1:length(Drug1[,2])) {
    if (Drug1[i,2]>100) {
      Drug1[i,2]<-100
    }
  }
  #remove missing observations
  Drug1<-Drug1[complete.cases(Drug1), ]
  
  #2nd Drug
  Drug2<-read.table(file = location2, header=T)
  #if response value is larger than 100%, convert it to 100%
  for (i in 1:length(Drug2[,2])) {
    if (Drug2[i,2]>100) {
      Drug2[i,2]<-100
    }
  }
  #remove missing observations
  Drug2<-Drug2[complete.cases(Drug2), ]
  
  ##---Calculate combination index---------------------
  Combination<-read.table(file = location4, header=T)
  #If response > 100, convert it to 100
  for (i in 1:length(Combination[,3])) {
    if (Combination[i,3]>100) {
      Combination[i,3]<-100
    }
  }
  #Eliminate missing observations in the Combination Dataset
  Comb <- Combination[complete.cases(Combination), ]
  
  
  SYN.test <- function(Drug1=Drug1, Drug2=Drug2, r=r)
  {
    # Purpose: to test the combinations of two drugs are additive 
    #
    # Input: 
    # Drug1: dataset of Drug 1, column 1-- dose, colummn 2--response     
    # Drug2: dataset of Drug 2, column 1-- dose, colummn 2--response  
    # Combination: dataset of Combination response of drug 1 and 2;
    #              column1--dose of drug1, column 2--dose of drug2, column 3--response 
    # r: number of replication at each dose-mixture 	
    #
    # Output:
    # Degrees of freedom of the F test
    # F-statistic
    # P-value
    
    coef1 <- as.numeric(lm(Drug1[,2] ~ log(Drug1[,1]))$coef)      #intercept and slope of drug 1
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)      #intercept and slope of drug 2
    a1 <- coef1[1]
    b1 <- coef1[2]
    a2 <- coef2[1]
    b2 <- coef2[2]
    if (coef1[2]<coef2[2]){		
      a1 <- coef2[1]
      b1 <- coef2[2]
      a2 <- coef1[1]
      b2 <- coef1[2]
    }
    
    #In the example Combination Dataset, r=6
    #Combination Dataset
    #If response > 100, convert it to 100

    #Eliminate missing observations in the Combination Dataset
    Comb <- Comb[complete.cases(Comb), ]
    n.Comb <- dim(Comb)[1]
    n <- n.Comb/r           # n is the number of experimental units 
    
    M <- Comb[order(Comb[, 1], Comb[, 2]), 1:3]
    x1 <- M[, 1]
    x2 <- M[, 2]
    y <- M[, 3]
    
    k <- b2/b1
    
    rho <- exp((a2 - a1)/b1)
    Z1 <- x1 + 0.5 * x2^k * (1 + sqrt(1 + (4 * (k - 1) * x1)/(rho * x2^k)))
    Z2 <- x1/Z1
    Z <- matrix(c(a1 + b1 * log(Z1), b1 * log(Z2)), ncol = 2, byrow = F)
    
    if (b1 < b2){
      rho <- exp((a1 - a2)/b2)
      Z1 <- x2 + 0.5 * x1^(1/k) * (1 + sqrt(1 + (4 * (1/k - 1) * x2)/(rho * x1^(1/k))))
      Z2 <- x2/Z1
      Z <- matrix(c(a2 + b2 * log(Z1), b2 * log(Z2)), ncol = 2, byrow = F)
    }
    
    U <- kronecker(diag(rep(1, n)), rep(1, r))
    J <- U %*% solve(t(U) %*% U) %*% t(U)
    V <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    
    A <- diag(rep(1, (n * r)))
    DF <- c(n - 2, (n * r - n))
    F.value <- ((n * r - n) * (t(y) %*% (J - V) %*% y))/((n - 2) * t(y) %*% (A - J) %*% y)
    P.value <- 1 - pf(F.value, (n - 2), (n * r - n), ncp = 0)
    
    F.statistic <- as.numeric(F.value)
    return(list(DF=DF, F.statistic=F.statistic, P.value=P.value))
  }
  
  SYN.test<-SYN.test(Drug1=Drug1, Drug2=Drug2, r=r)
  #In the example Combination Dataset, r=6
  
  
  SYN.analysis <- function(Drug1=Drug1, Drug2=Drug2, r=r, method = "loglog")
  {
    # Purpose: to fit combination index surface of two drugs the single dose responses are known
    # Input: 
    # Drug1: dataset of Drug 1, column 1-- dose, colummn 2--response     
    # Drug2: dataset of Drug 2, column 1-- dose, colummn 2--response  
    # Combination: dataset of Combination response of drug 1 and 2;
    #              column1--dose of drug1, column 2--dose of drug2, column 3--response 
    # r: number of replication at each dose-mixture 	
    # levels: levels of contour plot, for example, c(0, 0.3, 0.5, 1, 1.5, 2)
    # method: specifying the fitted models for the drug data;
    #         "loglog" for linear log models on both drugs,
    #         "lili" for linear models on both drugs.
    #
    # Output:
    # Contour Plot of the interaction index surface; 
    # Red lines are 95% confidence intervals of the additive action (combination index=1)
    # fit the drug data with linear log & linear models
    if ( method == "logli") {
      coef1 <- as.numeric( lm( Drug1[,2] ~ log(Drug1[,1]))$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ Drug2[,1] )$coeff )
    }
    
    # fit the drug data with linear & linear log models
    if ( method == "lilog") {
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1]))$coeff )
    }
    
    # fit the drug data with linear models
    if ( method == "lili" ) {
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ Drug2[,1] )$coeff )
    }
    
    # fit the drug data with linear log models
    if ( method == "loglog") {
      ##---Calculate parameters---------------------------	
      coef1 <- as.numeric(lm(Drug1[,2] ~ log(Drug1[,1]))$coef)      #intercept and slope of drug 1
      coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)      #intercept and slope of drug 2
    }	
    
    ##---Calculate combination index---------------------

    
    SD1 <- exp((Comb[, 3] - coef1[1])/coef1[2])
    SD2 <- exp((Comb[, 3] - coef2[1])/coef2[2])
    Index <- Comb[, 1]/SD1 + Comb[, 2]/SD2
    W <- cbind(Comb, Index)
    
    M <- as.matrix(W[order(W[, 1], W[, 2]), 1:4])
    y <- M[, 4]
    
    n.Comb <- dim(Comb)[1]
    m <- n.Comb/r
    
    A0 <- matrix(0, m, 2)
    for(i in 1:m) {A0[i,] <- M[((i - 1) * r + 1),1:2]}
    qrstr <- qr(cbind(rep(1, m), A0))
    A <- qr.Q(qrstr, complete = T)[,4:m]
    
    B0 <- matrix(0, m, m)
    for(j in 1:m) B0[,j] <- (matrix(rep(A0[j,],m),ncol=2,byrow=T)-A0)^2 %*% c(1,1)	
    
    Pert <- 1E-7
    ind <- abs(B0) < Pert
    B0[ind] <- Pert
    temp <- (B0 * log(B0 + diag(rep(1, m))))/(16 * pi)
    B <- t(A) %*% temp %*% A
    BC <- chol((B + t(B))/2)
    inver.BC <- solve(BC)
    
    Z <- kronecker(temp, rep(1, r)) %*% A %*% inver.BC
    C <- cbind(rep(1,n.Comb),M[,1:2],Z)
    D <- diag(c(0, 0, 0, rep(1, m-3)))
    
    ##---Estimate parameters----------------------------------
    si2u <- si2e <- 2
    n.iteration <- 150
    
    y<-as.matrix(y,nrow=ncol(J),ncol=1)
    for(k in 1:n.iteration) {		 
      J <- solve(t(C) %*% C + (si2e/si2u) * D)
      G <- J %*% t(C) %*% y
      U <- G[4:m]
      si2u <- c(t(U) %*% U + si2e * sum(diag(J)[4:m]))/(m-3)
      si2e <- c(t(y) %*% y - 2 * t(y) %*% C %*% G + t(G) %*% t(C) %*% C %*% G + si2e * sum(diag(t(C) %*% C %*% J)))/n.Comb
    }
    
    ##---Contour plot-------------------------------------------
    Lam <- si2e/si2u
    J <- solve(t(C) %*% C + Lam * D)
    G <- J %*% t(C) %*% y
    
    n.grid <- 100
    x1 <- seq(min(M[, 1]), max(M[, 1]), length = n.grid)
    x2 <- seq(min(M[, 2]), max(M[, 2]), length = n.grid)
    W <- W1 <- W2 <-matrix(0, n.grid, n.grid)
    
    Temp1 <- (matrix(rep(x1,m),ncol=n.grid,byrow=T)-matrix(rep(A0[,1],n.grid),ncol=n.grid,byrow=F))^2
    Temp2 <- (matrix(rep(x2,m),ncol=n.grid,byrow=T)-matrix(rep(A0[,2],n.grid),ncol=n.grid,byrow=F))^2
    
    ind <- 1
    Mat <- matrix(0,n.grid^2,m)
    for(j in 1:n.grid) {
      for(i in 1:n.grid) {
        Mat[ind,] <- Temp1[,i]+Temp2[,j]
        ind <- ind+1
      }
    }
    ttemp <- (Mat * log(Mat))/(16 * pi)
    z0 <- ttemp %*% A %*% inver.BC
    Cd <- cbind(rep(1,n.grid^2), as.matrix(expand.grid(x1, x2)), z0)
    W <- matrix(Cd %*% G,ncol=n.grid,byrow=F) 
    
    ind <- 1 
    for (j in 1:n.grid){
      CC <- Cd[ind:(ind+n.grid-1),]
      W1[, j] <- W[,j] + qnorm(0.975, 0, 1) * sqrt(si2e) * sqrt(diag(CC %*% J %*% t(CC)))
      W2[, j] <- W[,j] + qnorm(0.025, 0, 1) * sqrt(si2e) * sqrt(diag(CC %*% J %*% t(CC)))
      ind <- ind+n.grid
    }
    
    
    cont.lineposi <- contourLines(x1, x2, W1)
    cont.lineneg <- contourLines(x1, x2, W2)
    cont.linecent <- contourLines(x1, x2, W)
    
    #contour(x1, x2, W, levels = levels, xlab = "Drug 1 Dose", ylab = "Drug 2 Dose", 
    #xlim = range(x1), ylim = range(x2))	
    #x11()
    
    color_temp.palette <- palette( c("#00FF00", "#00A000", "yellow", "#FF0000", "#800000" ))
    filled.contour(x=x1, y=x2, z=W, xlab = "Drug 1 Dose", ylab = "Drug 2 Dose", 
                   col = color_temp.palette, color.palette = cm.colors, xlim = range(x1), ylim = range(x2), zlim = range(W), 
                   nlevels = 6)
    
    #Another way to plot the contour data: Use plotly package
    #install.packages("plotly")
    #library(plotly)
    #packageVersion('plotly')
    #p <- plot_ly(x=x1,y=x2,z = W, col = color_temp.palette, type = "contour")
    #plot the response surface
    Comb1 <-Comb[,1]
    Comb2 <- Comb[,2]
    resp <- Tps(Comb[,1:2], Comb[,3], m=2)
    Drug1.plot <- seq(min(Comb1), max(Comb1), length=30)
    Drug2.plot <- seq(min(Comb2), max(Comb2), length=30)
    plotgrid <- expand.grid(Comb1=Drug1.plot, Comb2=Drug2.plot)
    resp.fit <- matrix(predict(resp, plotgrid), 30, 30)
    resp.fit[resp.fit <= 0] <- 0;
    x11()
    persp(Drug1.plot, Drug2.plot, resp.fit, theta = 120, phi = 30, expand = 0.5, col = "pink",
          ltheta = 120, shade = .75, ticktype = "detailed",
          xlab = "Drug1", ylab = "Drug2", zlab = "Viability(%)", zlim = c(0, 1.1*max(resp.fit)),
          main = "The response surface of the combination of Drug1 and Drug2")
    
    #Plotting the confidence intervals
    mar.orig <- par("mar")
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
    
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
  }
  SYN.analysis <- SYN.analysis(Drug1=Drug1, Drug2=Drug2, r=6, method = method)
  
  #Output the results of SYN.test
  SYN.test.Dataframe<-as.data.frame(SYN.test)
  
  print(paste0("The Degrees of Freedom of the F distribution are "," ",SYN.test.Dataframe[1,1]," ","and"," ", SYN.test.Dataframe[2,1] ))
  print(paste0("The F-statistic is"," ", round(SYN.test$F.statistic,3)))
  print(paste0("The P-value is"," ", round(SYN.test$P.value,3)))
}






