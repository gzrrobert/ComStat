#' Identifying the linear/log-linear model of each dose-response curve
#' 
#' This function utilized nls() function to model the HILL model of each relationship and substitute each relationship
#' with linear or log-linear relationship, it then designs a dose-mixture combination dataset for the study by using Uniform design.
#' @param location1 location of Drug 1's dataset on your hard disk
#' @param location2 location of Drug 2's dataset on your hard disk
#' @param location3 location of UniformDesign on your hard disk
#' @param c0 the smallest meaningful difference to be detected
#' @param r number of replications at each dose-mixture, and usually 3<=r<=8, default=3
#' @return The fitted models and method for two drugs, 
#'         the number of dose-mixtures for combination, 
#'         total sample size,
#'         fitted plots for two drugs,
#'         dose mixtures for experimental design
#' @examples 
#' # Try this function with the following example:
#' TwoDrug_identify_design1<-TwoDrug_identify_design(
#'  location1="C:/Users/58412/Desktop/SRA/Two drugs/Drug1.txt",
#'  location2="C:/Users/58412/Desktop/SRA/Two drugs/Drug2.txt",
#'  location3="C:/Users/58412/Desktop/SRA/Two drugs/UniformDesign.txt",
#'  c0=0.01, r=6)
#' @export


#Econ=100
#HILL MODEL can be transformed into a linear regression model
#mlogD-mlogIC50=log[(Y-b)/(Econ-b)]
#Denote IC50 as IC_50
#viability range 20%-80% will be used in order to fit linear & log model
#data format: 1st column should be named as 'DOSE', 2nd column should be named as 'Viability'
#All missing data will be eliminated when using this program
#Attention: All datasets should be in TXT format

#' @import BB
#' @import nleqslv
#' @import ktsolve
#' @import drc


TwoDrug_identify_design<- function (location1, location2, location3, c0, r) {
  
  #location1: location of Drug 1's dataset on your hard disk
  #Drug1: data set of Drug 1, column 1-- dose, column 2--response  
  #location2: location of Drug 2's dataset on your hard disk
  #Drug2: data set of Drug 2, column 1-- dose, column 2--response
  #location3: location of UniformDesign on your hard disk
  # c0: the smallest meaningful difference to be detected
  # r: number of replications at each dose-mixture, and usually 3<=r<=5
  #All dataset should  be in txt file
  location1<-location1
  location2<-location2
  location3<-location3
  c0<-c0
  r<-r
  

  
  #1st Drug
  data<-read.table(file = location1, header=T)
  #if response value is larger than 100%, convert it to 100%
  for (i in 1:length(data$Viability)) {
    if (data[i,"Viability"]>100) {
      data[i,"Viability"]<-100
    }
  }
  
  #remove missing observations and reorder the data
  data<-data[complete.cases(data), ]
  data<-data[order(data[,1],decreasing=FALSE),]
  Dose<-data$Dose
  Viability<-data$Viability
  
  #Intercept b is the minimum viability oberserved
  b<-min(Viability)
  
  #Reformulate the Hill model into a linear model
  #Now the response is:
  y<-matrix(NA,length(Viability),1)
  Viability.matrix<-as.matrix(Viability)
  for (i in 1:length(Viability.matrix)) {
    if (Viability.matrix[i,1]>b) {
      y[i,1]<-log((Viability.matrix[i,1]-b)/(100-b))
    } else {
      break
    }
  }
  y<-y[!is.na(y),]
  y<-data.frame((y))
  
  data.length<-dim(y)[1]
  Dose<-as.matrix(Dose)
  Dose<-Dose[1:data.length,]
  Dose<-as.numeric(Dose)
  
  data<-data[1:data.length,]
  
  data.transform<-as.data.frame(cbind(Dose,y$X.y.))
  
  #fit the Hill model  
  #Estimate the unknown coefficients in the HILL model
  model<-nls(V2~m*log(Dose)+m*log(IC_50),start=list(m=-20, IC_50=1.1),data=data.transform,trace=F)
  
  #the coefficients of optimized model
  coef<-summary(model)$coefficients
  m<-coef[1,1]
  IC_50<-coef[2,1]
  
  #Build the optimized Hill Model
  #y=(((Econ-b)(Dose/IC_50)^m)/(1+(Dose/IC_50)^m))+b
  #And calculate the fitted response value
  Viability_hat<-matrix(NA,length(Dose),1)
  Econ=100
  for (i in 1:length(Dose)) {
    Viability_hat[i]<-(((Econ-b)*((data[i,1]/IC_50)**m))/(1+(data[i,1]/IC_50)**m))+b
  }
  Viability_hat<-data.frame(Dose,Viability_hat)
  
  #visual check of the hill model
  par(mfrow=c(1,2))
  fit.ll <- drm(Viability~Dose, data=data, fct=LL.5(), type="continuous")
  plot(fit.ll,type = "none",col='red',xlab = "Dose",ylab = "Viability(%)",main = "Drug 1")
  points(x=data$Dose,y=data$Viability)
  #Can use ggplot() to achieve a more smooth curve on the scattered plot
  
  #Create a dataset with residuals that can be used to fit linear or loglinear model
  #Viability range from 20% to 80% will be used;
  residual_data_2080_range<-data[data$Viability>=20 & data$Viability<=80,]
  residual_data_2080_range_dose<-residual_data_2080_range$Dose
  residual<-as.data.frame(Viability_hat$Viability_hat)-as.matrix(data$Viability)
  residual<-data.frame(data$Dose,residual)
  residual_data<-residual[max(residual_data_2080_range_dose)>=residual$data.Dose & residual$data.Dose>=min(residual_data_2080_range_dose),]
  
  
  #Decide either Linear or Log model suffice
  #1.fit the data with linear model
  #Note: Only the response range from 20% to 80% in the dataset with residuals will be used to fit both linear and log model
  linear<-lm(Viability_hat.Viability_hat~data.Dose,data=residual_data)
  if (summary(linear)$coefficients[2,4]>0.05) {
    print("Warning: the linear model may not suffice")
  }
  #Calculate fitted responses under linear model
  linear_hat<-matrix(NA,length(residual_data$data.Dose),1)
  for (i in 1:length(residual_data$data.Dose)){
    linear_hat[i]<-(summary(linear)$coefficients[1,1])+(summary(linear)$coefficients[2,1])*residual_data[i,1]
  }
  #2.fit the data with log-linear model
  log<-lm(Viability_hat.Viability_hat~log(data.Dose),data=residual_data)
  if (summary(log)$coefficients[2,4]>0.05) {
    print("Warning: the log model may not suffice")
  }
  #Calculate fitted responses under log model
  log_hat<-matrix(NA,length(residual_data$data.Dose),1)
  for (i in 1:length(residual_data$data.Dose)){
    log_hat[i]<-(summary(log)$coefficients[1,1])+(summary(log)$coefficients[2,1])*log(residual_data[i,1])
  }
  
  #The fitted response from Hill Model
  #range Viability 20%-80%
  hill_hat<-Viability_hat[max(residual_data_2080_range_dose)>=residual$data.Dose & residual$data.Dose>=min(residual_data_2080_range_dose),]
  hill_hat<-as.matrix(hill_hat[ ,2])
  
  #linear v.s Hill
  diff_linear_hill<-matrix(NA,length(linear_hat),1)
  for (i in 1:length(linear_hat)) {
    diff_linear_hill[i]<-(linear_hat[i]-hill_hat[i])**2
  }
  #Sum of square of Difference between linear model and Hill model
  SSD_linear_Hill<-0
  for (i in 1:length(diff_linear_hill)) {
    SSD_linear_Hill<-diff_linear_hill[i]+SSD_linear_Hill
  }
  #log v.s Hill
  diff_log_hill<-matrix(NA,length(log_hat),1)
  for (i in 1:length(log_hat)) {
    diff_log_hill[i]<-(log_hat[i]-hill_hat[i])**2
  }
  #Sum of square of Difference between linear model and Hill model
  SSD_log_Hill<-0
  for (i in 1:length(diff_log_hill)) {
    SSD_log_Hill<-diff_log_hill[i]+SSD_log_Hill
  }
  #Compare the linear model with log model
  #Can use the following ifelse() function as well.
  #ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Linear Model should suffice"), ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Log Model should suffice"), print("Either Log model or lieanr model suffices")))
  
  
  
  #2nd Drug
  data2<-read.table(file = location2, header=T)
  #if response value is larger than 100%, convert it to 100%
  for (i in 1:length(data2$Viability)) {
    if (data2[i,"Viability"]>100) {
      data2[i,"Viability"]<-100
    }
  }
  
  #remove missing values and reorder the data
  data2<-data2[complete.cases(data2), ]
  data2<-data2[order(data2[,1],decreasing=FALSE),]
  Dose2<-data2$Dose
  Viability2<-data2$Viability
  
  #Intercept b is the minimum viability oberserved
  b2<-min(Viability2)
  
  #Reformulate the Hill model into a linear model
  #The response is now:
  y2<-matrix(NA,length(Viability2),1)
  Viability.matrix2<-as.matrix(Viability2)
  for (i in 1:length(Viability.matrix2)) {
    if (Viability.matrix2[i,1]>b2) {
      y2[i,1]<-log((Viability.matrix2[i,1]-b2)/(100-b2))
    } else {
      break
    }
  }
  y2<-y2[!is.na(y2),]
  y2<-data.frame((y2))
  
  data.length2<-dim(y2)[1]
  Dose2<-as.matrix(Dose2)
  Dose2<-Dose2[1:data.length2,]
  Dose2<-as.numeric(Dose2)
  
  data2<-data2[1:data.length2,]
  
  data.transform2<-as.data.frame(cbind(Dose2,y2$X.y2.))
  
  #fit the Hill model  
  #Estimate the coefficients in the HILL model
  model2<-nls(V2~m*log(Dose2)+m*log(IC_50),start=list(m=-20, IC_50=1.1),data=data.transform2,trace=F)
  
  #the coefficients of optimized model
  coef2<-summary(model2)$coefficients
  m2<-coef2[1,1]
  IC_50_2<-coef2[2,1]
  
  #Build the optimized Hill Model
  #y=(((Econ-b)(Dose/IC_50)^m)/(1+(Dose/IC_50)^m))+b
  #And calculated the fitted response
  Viability_hat2<-matrix(NA,length(Dose2),1)
  Econ=100
  for (i in 1:length(Dose2)) {
    Viability_hat2[i]<-(((Econ-b2)*((data2[i,1]/IC_50_2)**m2))/(1+(data2[i,1]/IC_50_2)**m2))+b2
  }
  Viability_hat2<-data.frame(Dose2,Viability_hat2)
  
  #visual check of the hill model
  fit.ll2 <- drm(Viability~Dose, data=data2, fct=LL.5(), type="continuous")
  plot(fit.ll2,type = "none",col='red',xlab = "Dose",ylab = "Viability(%)",main = "Drug 2")
  points(x=data2$Dose,y=data2$Viability)
  #Can use ggplot to achieve a more smooth curve on the scattered plot
  
  #Create a dataset with residuals
  #To be used to fit linear|log model
  #Viability range from 20% to 80% will be used;
  residual_data_2080_range2<-data2[data2$Viability>=20 & data2$Viability<=80,]
  residual_data_2080_range_dose2<-residual_data_2080_range2$Dose
  residual2<-as.data.frame(Viability_hat2$Viability_hat)-as.matrix(data2$Viability)
  residual2<-data.frame(data2$Dose,residual2)
  residual_data2<-residual2[max(residual_data_2080_range_dose2)>=residual2$data2.Dose & residual2$data2.Dose>=min(residual_data_2080_range_dose2),]
  
  
  #Decide either Linear or Log model suffice
  #1.fit the data with linear model
  #Note: Only the response range from 20% to 80% in the dataset with residuals will be used to fit both linear and log model
  linear2<-lm(Viability_hat2.Viability_hat~data2.Dose,data=residual_data2)
  if (summary(linear2)$coefficients[2,4]>0.05) {
    print("Warning: the linear model may not suffice")
  }
  #Calculate fitted responses under linear model
  linear_hat2<-matrix(NA,length(residual_data2$data2.Dose),1)
  for (i in 1:length(residual_data2$data2.Dose)){
    linear_hat2[i]<-(summary(linear2)$coefficients[1,1])+(summary(linear2)$coefficients[2,1])*residual_data2[i,1]
  }
  #2.fit the data with log-linear model
  log2<-lm(Viability_hat2.Viability_hat~log(data2.Dose),data=residual_data2)
  if (summary(log2)$coefficients[2,4]>0.05) {
    print("Warning: the log model may not suffice")
  }
  #Calculate fitted responses under log model
  log_hat2<-matrix(NA,length(residual_data2$data2.Dose),1)
  for (i in 1:length(residual_data2$data2.Dose)){
    log_hat2[i]<-(summary(log2)$coefficients[1,1])+(summary(log2)$coefficients[2,1])*log(residual_data2[i,1])
  }
  
  #The fitted response from Hill Model
  #range Viability 20%-80%
  hill_hat2<-Viability_hat2[max(residual_data_2080_range_dose2)>=residual2$data2.Dose & residual2$data2.Dose>=min(residual_data_2080_range_dose2),]
  hill_hat2<-as.matrix(hill_hat2[ ,2])
  
  #linear v.s Hill
  diff_linear_hill2<-matrix(NA,length(linear_hat2),1)
  for (i in 1:length(linear_hat2)) {
    diff_linear_hill2[i]<-(linear_hat2[i]-hill_hat2[i])**2
  }
  #Sum of square of Difference between linear model and Hill model
  SSD_linear_Hill2<-0
  for (i in 1:length(diff_linear_hill2)) {
    SSD_linear_Hill2<-diff_linear_hill2[i]+SSD_linear_Hill2
  }
  
  #log v.s Hill
  diff_log_hill2<-matrix(NA,length(log_hat2),1)
  for (i in 1:length(log_hat2)) {
    diff_log_hill2[i]<-(log_hat2[i]-hill_hat2[i])**2
  }
  #Sum of square of Difference between linear model and Hill model
  SSD_log_Hill2<-0
  for (i in 1:length(diff_log_hill2)) {
    SSD_log_Hill2<-diff_log_hill2[i]+SSD_log_Hill2
  }
  #Compare the linear model with log model
  #Can use the following ifelse() function as well.
  #ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Linear Model should suffice"), ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Log Model should suffice"), print("Either Log model or lieanr model suffices")))
  
  if (SSD_log_Hill>SSD_linear_Hill) {
    print("The Linear Model should suffice for Drug 1")
    model_drug1="linear"
    drug1_a1<-summary(linear)$coefficients[1,1]
    drug1_b1<-summary(linear)$coefficients[2,1]
  } else if (SSD_log_Hill<SSD_linear_Hill) {
    print("The Log Model should suffice for Drug 1")
    model_drug1="log"
    drug1_a1<-summary(log)$coefficients[1,1]
    drug1_b1<-summary(log)$coefficients[2,1]
  } else {
    print("Either Log model or lieanr model suffices for Drug 1")
    model_drug1="linear or log"
  }
  
  if (SSD_log_Hill2>SSD_linear_Hill2) {
    print("The Linear Model should suffice for Drug 2")
    model_drug2<-"linear"
    drug2_a2<-summary(linear2)$coefficients[1,1]
    drug2_b2<-summary(linear2)$coefficients[2,1]
  } else if (SSD_log_Hill2<SSD_linear_Hill2) {
    print("The Log Model should suffice for Drug 2")
    model_drug2<-"log"
    drug2_a2<-summary(log2)$coefficients[1,1]
    drug2_b2<-summary(log2)$coefficients[2,1]
  } else {
    print("Either Log model or lieanr model suffices for Drug 2")
    model_drug2<-"linear or log"
  }
  
  
  # method: specifying the fitted models for the drug data;
  #"loglog" for linear log models on both drugs,
  #"lili" for linear models on both drugs,
  #"lilog" for linear model on Drug1 and log model on Drug2
  #"logli" for log model on Drug1 and linear model on Drug2
  if (model_drug1=="linear" & model_drug2=="linear") {
    method<-"lili"
  } else if (model_drug1=="linear" & model_drug2=="log") {
    method<-"lilog"
  } else if (model_drug1=="log" & model_drug2=="linear") {
    method<-"logli"
  } else if (model_drug1=="log" & model_drug2=="log") {
    method<-"loglog"
  }
  #stop("Please enter the correct location of your dataset on the disk. Note:Use left slash sign after every folder names.")

  Drug1<-data
  Drug2<-data2
  method<-method

  UniformDesign<-read.table(file = location3, sep=',', header=T)

#example
#Try this the function by specifying the data location on your hard disk
#TwoDrug<-TwoDrug(location1= "C:/Users/yinru/Dropbox/SRA 2016Fall/SynergyProgram/RProgram/Drug1.txt",
#                location2= "C:/Users/yinru/Dropbox/SRA 2016Fall/SynergyProgram/RProgram/Drug2.txt" )


SYN.design <- function(Drug1=Drug1, Drug2=Drug2, h=c(20,80), c0, r, method = "loglog")
{
  # Input:
  # data: single-dose experiments of two drugs
  #       Drug1: data set of Drug 1, column 1-- dose, column 2--response  
  #       Drug2: data set of Drug 2, column 1-- dose, column 2--response  
  # h=c(h1, h2): experimental domain for given low (h1) and high (h2) response boundary 
  #              for example: h1=20, h2=80 
  # c0: the smallest meaningful difference to be detected
  # r: number of replications at each dose-mixture, and usually 3<=r<=5
  # Note: m*r<=30
  # method: specifying the fitted models for the drug data;
  #         "loglog" for linear log models on both drugs,
  #         "lili" for linear models on both drugs,
  #         "lilog" for linear model on Drug1 and log model on Drug2
  #         "logli" for log model on Drug1 and linear model on Drug2
  #
  # Output: 
  # Mixture number (m): number of dose-mixtures for combination study
  # Experiment number (m*r): total sample size
  # Discrepancy: the central L2-discrepancy of m mixtures
  # Dose.mixtures: dose mixtures for combination experiments 
  #
  # Default:
  # data set "UniformDesign" is needed
  # power: 80%
  # significance level: 5%
 
  
  if ( method != "loglog" & method != "lili" & method != "lilog" & method != "logli") { 
    stop( "method = unknown regression model." ) 
  }
  
  ##--- estimation of the single dose-responses
  if ( method == "loglog") {
    coef1 <- as.numeric(lm(Drug1[,2] ~ log(Drug1[,1]))$coef)      #intercept and slope of drug 1 
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)      #intercept and slope of drug 2 	
  }
  if ( method == "lili") {
    coef1 <- as.numeric((lm(Drug1[,2] ~ Drug1[,1]))$coef)
    coef2 <- as.numeric((lm(Drug2[,2] ~ Drug2[,1]))$coef)
  }
  if ( method == "logli" ) {
    method <- "lilog"
    Drug.temp <- Drug1
    Drug1 <- Drug2
    Drug2 <- Drug.temp
  }
  if ( method == "lilog") {
    coef1 <- as.numeric((lm(Drug1[,2] ~ Drug1[,1]))$coef)
    coef2 <- as.numeric(lm(Drug2[,2] ~ log(Drug2[,1]))$coef)
  }
  
  si2 <- var(c(Drug1[, 2], Drug2[, 2]))     #si2: variance estimated from the single experiment data
  ##--- calculate the sample size for uniform design
  #m is the number of mixtures
  #k is the number of drugs, in this two drug combination case, k=2
  #3<=r<=8
  #3<=m<=30
  #k is automatically set to 2 since this is a two drug combination study
  #d <- c0^2/si2
  #delta0: non-central parameter of the F-distribution
  #delta0 <- (m * (r+1)) * d
  #To obtain m, Step 1: Obtain the quantile under F-distribution
  #Note: alpha (significance level)=0.05
  #Substitute quantile with the following equation, otherwise an error occurs since m is the unknown parameter
  #and has to be estimated
  #quantile<-qf(0.95,(m-2),(m*r))
  #Step 2: Calculate the non-central parameter as delta0
  #Substitute delta0 with the following equation, otherwise an error occurs since m is the unknown parameter
  #and has to be estimated
  #Since d <- c0^2/si2 and delta0 <- (m * (r+1)) * d, then delta0 can be written into:
  #delta0 <- (m * (r+1)) * (c0^2/si2)
  #Step 3: Use pre-determined power (80%) to determine m
  #Constructing essential elements to be used in 'ktsolve' function
  r<-r
  yfunc <- function(x) {
    y <- vector()
    y[1] <- exp(-(m * (r+1)) * (c0^2/si2)/2)*pf(qf(0.95,(m-2),(m*r)),m-2,m*r)+
      exp(-(m * (r+1)) * (c0^2/si2)/2)*((m * (r+1)) * (c0^2/si2)/2)*pf(((m-2)/m)*qf(0.95,(m-2),(m*r)), m, m*r)+
      ((exp(-(m * (r+1)) * (c0^2/si2)/2)*((m * (r+1)) * (c0^2/si2)/2)^2)/2)*pf(((m-2)/(m+2))*qf(0.95,(m-2),(m*r)), (m+2), m*r)-0.2
  }
  known=list(r=r,c0=c0,si2=si2)
  guess=list(m=3)
  #ktsolve function
  solvm <- ktsolve(yfunc, known=known, guess=guess, tool=c("BB","nleqslv"), show=TRUE)
  if(ceiling(solvm$results$par)>30 | ceiling(solvm$results$par)<3) { m<-NULL
  print("Please choose a larger r (Iteration) and rerun this program")} else 
  {m<-ceiling(solvm$results$par)}
  
  
  ##--- Specify the U-type matrix location on your hard disk
  UD <- UniformDesign
  Dis <- UD[1, m]
  MU <- matrix(c(UD[2:(m + 1), 1], UD[2:(m + 1), m]), ncol = 2, byrow = F)
  
  ##--- choose the experimental domain
  # the total dose is from e1 to e2
  
  a1 <- coef1[1]
  b1 <- coef1[2]
  a2 <- coef2[1]
  b2 <- coef2[2]
  Variables <- c("Drug1", "Drug2")
  
  if (coef1[2]<coef2[2]){		
    a1 <- coef2[1]
    b1 <- coef2[2]
    a2 <- coef1[1]
    b2 <- coef1[2]
    Variables <- c("Drug2", "Drug1")
  }
  
  rho <- exp((a2 - a1)/b1)
  e1 <- exp((h[2] - a1)/b1)
  e2 <- exp((h[1] - a1)/b1)
  
  ##--- find the m dose-mixtures
  if ( method == "lili" | method == "lilog") {
    A <- matrix(0, 2, 4)
    A[, 1] <- c(c((h[2] - a1)/b1), 0)
    A[, 2] <- c(c((h[1] - a1)/b1), 0)
    A[, 3] <- c(0, c((h[1] - a1)/b2))
    A[, 4] <- c(0, c((h[2] - a1)/b2))
    d0 <- det(cbind(A[, 2] - A[, 1], A[, 4] - A[, 1]))
    d1 <- det(cbind(A[, 2] - A[, 1], A[, 3] - A[, 4]))
    d2 <- det(cbind(A[, 3] - A[, 2], A[, 4] - A[, 1]))
    y1 <- ( - d0 + sqrt(d0^2 + 2 * d1 * c((MU[, 1] - 0.5)/m) * (d0 + 0.5 * d1)))/d1
    y2 <- c((MU[, 2] - 0.5)/m)
    Z <- matrix(0, m, 2)
    for(i in 1:m) {
      Z[i,  ] <- A %*% c((1 - y1[i]) * (1 - y2[i]), y1[i] * (1 - y2[i]), y1[i] * y2[i], (1 - y1[i]) * y2[i])
    }
    z1 <- c(Z[, 1])
    z2 <- c(Z[, 2])
    X <- matrix(0, m, 2)
    X[, 1] <- z1
    if ( method == "lili") {
      X[, 2] <- z2 - (a2 - a1)/b2
    }
    if ( method == "lilog") {
      X[, 2] <- (z2 * exp((a1 - a2)/b2) * exp(z1 * (b1/b2) + z2))/(z1 * (b1/b2) + z2)
    }
    Dose.mixtures <- X
  }
  
  if ( method == "loglog") {
    y1 <- ((MU[, 1] - 0.5)/m)*(e2-e1)+e1
    y2 <- (MU[, 2] - 0.5)/m
    d1 <- y1*y2
    d2 <- rep(0, m)
    D <- seq(0.001, (2 * e2), 0.001)
    k <- b2/b1
    for(j in 1:m) {
      f1 <- (2 * (y1[j] - d1[j]) - D^k * (1 + sqrt(1 + (4 * (k - 1) * d1[j])/(D^k * rho))))^2
      Ind <- sort.list(f1)[1]
      d2[j] <- ((D[Ind])^(1/k))* exp((a1 - a2)/b2)
    }
    Dose.mixtures <- matrix(c(d1, d2), ncol = 2, byrow = F)
  }
  dimnames(Dose.mixtures) <- list(NULL, Variables)
  
  return(list(Mixture.number=m, Experiment.number=m*r, Discrepancy=Dis, Dose.mixtures=round(Dose.mixtures,3)))	
}

SYN.design <- SYN.design(
Drug1=Drug1, Drug2=Drug2, h=c(20,80), c0=c0, r=r, method = "loglog")

print(paste0("The method used is"," ",method))
print(paste0("The number of combinations (mixtures) needed is m="," ", SYN.design$Mixture.number))
print(paste0("Total sample size at least is m*r="," ", SYN.design$Experiment.number))
print(paste0("The central L2-discrepancy of m mixtures is"," ", round(SYN.design$Discrepancy,3)))
print(paste0("Please use the following scheme of dose mixtures to measure the response, input the results into a txt file and name as 'Combination' "))
print(SYN.design$Dose.mixtures)
print("The Combination dataset should have the pattern as: column1--dose of drug1, column 2--dose of drug2, column 3--response")   
}

