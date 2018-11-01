#' Testing the additive effect
#' 
#' This function implements an F test to find out the additive/synergistic/antagonistic effect under each three
#' drug combination level
#' @param location1 location of Drug 1's dataset on your hard disk
#' @param location2 location of Drug 2's dataset on your hard disk
#' @param location3 location of Drug 3's dataset on your hard disk
#' @param location4 location of combination dataset on your hard disk
#' @param r number of replications at each dose-mixture, and usually 3<=r<=8
#' @return F test of additive/synergistic/antagonistic effect under drug combination
#' @examples
#' # Try this function with the following example:
#In the example Combination Dataset, r=6:
#' SY.test3<-ThreeDrug_test(
#' location1="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug1.txt", 
#' location2="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug2.txt", 
#' location3="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug3.txt", 
#' location4="C:/Users/58412/Desktop/SRA/Three drugs/data2/Combination.txt", 
#' r=6)
#' @export


ThreeDrug_test <- function(location1, location2, location3,location4, r)
{
  # Purpose: to test the combinations of three drugs are additive 
  #(equivalent to test each pair of three drugs are additive) 
  #
  # Input: 
  # Drug1: dataset of Drug 1, column 1-- dose, colummn 2--response     
  # Drug2: dataset of Drug 2, column 1-- dose, colummn 2--response  
  # Drug3: dataset of Drug 2, column 1-- dose, colummn 2--response  
  # Combination(location4): dataset of Combination response of drug 1 and 2;
  #              column1--dose of drug1, column 2--dose of drug2, column 3--dose of drug3, column 4--response 
  # r: number of replication at each dose-mixture   
  #
  # Output:
  # Degrees of freedom of the F test
  # F-statistic
  # P-value
  
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
  
  #3rd Drug
  Drug3<-read.table(file = location3, header=T)
  #if response value is larger than 100%, convert it to 100%
  for (i in 1:length(Drug3[,2])) {
    if (Drug3[i,2]>100) {
      Drug3[i,2]<-100
    }
  }
  #remove missing observations
  Drug3<-Drug3[complete.cases(Drug3), ]
  
  Combination<-read.table(file = location4, header=T)
  
  coef1 <- as.numeric(lm(Drug1[,2] ~ (Drug1[,1]))$coef)      #intercept and slope of drug 1
  coef2 <- as.numeric(lm(Drug2[,2] ~ (Drug2[,1]))$coef)      #intercept and slope of drug 2
  coef3 <- as.numeric(lm(Drug3[,2] ~ (Drug3[,1]))$coef)      #intercept and slope of drug 3
  
  a<-matrix(c(coef1,coef2,coef3),nrow=2,ncol=3,byrow=FALSE)
  b<-matrix(c(coef1[2],coef2[2],coef3[2]),nrow=1,ncol=3)
  maxbeta=max(b);
  minbeta=min(b);
  coef1<-a[, which(b==maxbeta)]
  coef3<-a[, which(b==minbeta)]
  coef2<-a[, which(b<maxbeta & b>minbeta)]
  
  a1 <- coef1[1]
  b1 <- coef1[2]
  a2 <- coef2[1]
  b2 <- coef2[2]
  a3 <- coef3[1]
  b3 <- coef3[2]
  
  
  Comb <- Combination
  n.Comb <- dim(Comb)[1]
  n <- n.Comb/r                          # n is the number of experimental units
  
  M <- Comb[order(Comb[, 1], Comb[, 2],Comb[, 3]), 1:4]
  x1 <- M[, 1]
  x2 <- M[, 2]
  x3 <- M[, 3]
  y <- M[, 4]
  
  k1 <- b1/b2
  k2 <- b1/b3
  rho0 <- exp((a2 - a1)/b1)
  rho1 <- exp((a3-a1)/b1)
  rho2 <- rho0/rho1
  h1<-x1/rho1+rho2^(k1)*k1*(b3-b2)*x2/(b3-b1)
  h2<-rho2^(k1)*b3*(b2-b1)*x2/(b2*(b3-b1))+x3
  psi1<-(1-k2)*h2-h1^(k2)
  psi2<-(((1-k2)*h2-h1^(k2))^2+2*k2*(1-k2)*h2^2)^(0.5)
  psi<-b3^2*h1/(b1*(b3-b1)*h2^2)*(psi1+psi2)
  Z1<-x1+rho2^(k1-1)*psi^(k1*(b3-b2)/(b3-b1))*x2+psi*x3
  dz2<-x1+rho2^(k1-1)*psi^(k1*(b3-b2)/(b3-b1))*x2+psi*x3*rho1/rho0
  Z2<-x1/dz2
  Z3<-psi*x3/Z1
  gz1<-a1+b1*log(Z1)
  gz2<-b1*log((1-rho0)*Z2+rho0)
  gz3<-b1*log((1-1/rho2)*(1-Z3)+1/rho2)
  Z <- matrix(c(gz1,gz2,gz3), ncol = 3, byrow = F)
  
  U <- kronecker(diag(rep(1, n)), rep(1, r))
  J <- U %*% solve(t(U) %*% U) %*% t(U)
  V <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  
  A <- diag(rep(1, (n * r)))
  DF <- c(n - 3, (n * r - n))
  F.value <- ((n * r - n) * (t(y) %*% (J - V) %*% y))/((n - 3) * t(y) %*% (A - J) %*% y)
  P.value <- 1 - pf(F.value, (n - 3), (n * r - n), ncp = 0)
  
  F.statistic <- as.numeric(F.value)
  return(list(DF=DF, F.statistic=round(F.statistic,3), P.value=round(P.value,3)))
}







