#' Identifying the linear/log-linear model of each dose-response curve
#' 
#' This function utilized nls() function to model the HILL model of each relationship and substitute each relationship
#' with linear or log-linear relationship, it then designs a dose-mixture combination dataset for the study by using Uniform design.
#' @param location1 location of Drug 1's dataset on your hard disk
#' @param location2 location of Drug 2's dataset on your hard disk
#' @param location3 location of Drug 3's dataset on your hard disk
#' @param c0 the smallest meaningful difference to be detected
#' @param r number of replications at each dose-mixture, and usually 3<=r<=8, default=3
#' @return The fitted models and method for three drugs, 
#'         the number of dose-mixtures for combination, 
#'         total sample size,
#'         fitted plots for three drugs,
#'         dose mixtures for experimental design
#' @examples 
#' # Try this function with the following examples:
#' #(please be patient when the method is"logloglog"):
#' ThreeDrug_identify_design_1<- ThreeDrug_identify_design(
#' location1="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug1.txt", 
#' location2="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug2.txt", 
#' location3="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug3.txt", 
#' c0=13, r=6)

#' ThreeDrug_identify_design_2<- ThreeDrug_identify_design(
#' location1="C:/Users/58412/Desktop/SRA/Three drugs/data1/ARAC03.txt", 
#' location2="C:/Users/58412/Desktop/SRA/Three drugs/data1/SAHA03.txt", 
#' location3="C:/Users/58412/Desktop/SRA/Three drugs/data1/VP1603.txt", 
#' c0=15, r=6)
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


#Econ=100
#HILL MODEL can be transformed into a linear regression model
#mlogD-mlogIC50=log[(Y-b)/(Econ-b)]
#Denote IC50 as IC_50
#viability range 20%-80% will be used in order to fit linear & log model
#data format: 1st column should be named as 'Dose', 2nd column should be named as 'Viability'
#All missing data will be eliminated when using this program
#Attention: All datasets should be in TXT format

ThreeDrug_identify_design<- function (location1, location2, location3, c0, r) {
  
  #location1: location of Drug 1's dataset on your hard disk
  #Drug1: data set of Drug 1, column 1-- dose, column 2--response  
  #location2: location of Drug 2's dataset on your hard disk
  #Drug2: data set of Drug 2, column 1-- dose, column 2--response
  #location3: tion of Drug 3's dataset on your hard disk
  #Drug3: data set of Drug 3, column 1-- dose, column 2--response
  # c0: the smallest meaningful difference to be detected
  # r: number of replications at each dose-mixture, and usually 3<=r<=5,sometimes can be 6
  #All dataset should  be in txt file
  
  ThreeDrugIdentify<-function (location1,location2, location3) {
    
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
    
    
    
    ####3rd Drug
    
    data3<-read.table(file = location3, header=T)
    #if response value is larger than 100%, convert it to 100%
    for (i in 1:length(data3$Viability)) {
      if (data3[i,2]>100) {
        data3[i,2]<-100
      }
    }
    
    #remove missing values and reorder the data
    data3<-data3[complete.cases(data3), ]
    data3<-data3[order(data3[,1],decreasing=FALSE),]
    Dose3<-data3$Dose
    Viability3<-data3$Viability
    
    #Intercept b is the minimum viability oberserved
    b3<-min(Viability3)
    
    #Reformulate the Hill model into a linear model
    #The response is now:
    y3<-matrix(NA,length(Viability3),1)
    Viability.matrix3<-as.matrix(Viability3)
    for (i in 1:length(Viability.matrix3)) {
      if (Viability.matrix3[i,1]>b3) {
        y3[i,1]<-log((Viability.matrix3[i,1]-b3)/(100-b3))
      } else {
        break
      }
    }
    y3<-y3[!is.na(y3),]
    y3<-data.frame((y3))
    
    data.length3<-dim(y3)[1]
    Dose3<-as.matrix(Dose3)
    Dose3<-Dose3[1:data.length3,]
    Dose3<-as.numeric(Dose3)
    
    data3<-data3[1:data.length3,]
    
    data.transform3<-as.data.frame(cbind(Dose3,y3$X.y3.))
    
    #fit the Hill model  
    #Estimate the coefficients in the HILL model
    model3<-nls(V2~m*log(Dose3)+m*log(IC_50),start=list(m=-20, IC_50=1.1),data=data.transform3,trace=F)
    
    #the coefficients of optimized model
    coef3<-summary(model3)$coefficients
    m3<-coef3[1,1]
    IC_50_3<-coef3[2,1]
    
    #Build the optimized Hill Model
    #y=(((Econ-b)(Dose/IC_50)^m)/(1+(Dose/IC_50)^m))+b
    #And calculated the fitted response
    Viability_hat3<-matrix(NA,length(Dose3),1)
    Econ=100
    for (i in 1:length(Dose3)) {
      Viability_hat3[i]<-(((Econ-b3)*((data3[i,1]/IC_50_3)**m3))/(1+(data3[i,1]/IC_50_3)**m3))+b3
    }
    Viability_hat3<-data.frame(Dose3,Viability_hat3)
    

    
    #Create a dataset with residuals
    #To be used to fit linear|log model
    #Viability range from 20% to 80% will be used;
    residual_data_2080_range3<-data3[data3[,2]>=20 & data3[,2]<=80,]
    residual_data_2080_range_dose3<-residual_data_2080_range3[,1]
    residual3<-as.data.frame(Viability_hat3$Viability_hat)-as.matrix(data3$Viability)
    residual3<-data.frame(data3$Dose,residual3)
    residual_data3<-residual3[max(residual_data_2080_range_dose3)>=residual3$data3.Dose & residual3$data3.Dose>=min(residual_data_2080_range_dose3),]
    
    
    #Decide either Linear or Log model suffice
    #1.fit the data with linear model
    #Note: Only the response range from 20% to 80% in the dataset with residuals will be used to fit both linear and log model
    linear3<-lm(Viability_hat3.Viability_hat~data3.Dose,data=residual_data3)
    if (summary(linear3)$coefficients[2,4]>0.05) {
      print("Warning: the linear model may not suffice")
    }
    #Calculate fitted responses under linear model
    linear_hat3<-matrix(NA,length(residual_data3$data3.Dose),1)
    for (i in 1:length(residual_data3$data3.Dose)){
      linear_hat3[i]<-(summary(linear3)$coefficients[1,1])+(summary(linear3)$coefficients[2,1])*residual_data3[i,1]
    }
    #2.fit the data with log-linear model
    log3<-lm(Viability_hat3.Viability_hat~log(data3.Dose),data=residual_data3)
    if (summary(log3)$coefficients[2,4]>0.05) {
      print("Warning: the log model may not suffice")
    }
    #Calculate fitted responses under log model
    log_hat3<-matrix(NA,length(residual_data3$data3.Dose),1)
    for (i in 1:length(residual_data3$data3.Dose)){
      log_hat3[i]<-(summary(log3)$coefficients[1,1])+(summary(log3)$coefficients[2,1])*log(residual_data3[i,1])
    }
    
    #The fitted response from Hill Model
    #range Viability 20%-80%
    hill_hat3<-Viability_hat3[max(residual_data_2080_range_dose3)>=residual3$data3.Dose & residual3$data3.Dose>=min(residual_data_2080_range_dose3),]
    hill_hat3<-as.matrix(hill_hat3[ ,2])
    
    #linear v.s Hill
    diff_linear_hill3<-matrix(NA,length(linear_hat3),1)
    for (i in 1:length(linear_hat3)) {
      diff_linear_hill3[i]<-(linear_hat3[i]-hill_hat3[i])**2
    }
    #Sum of square of Difference between linear model and Hill model
    SSD_linear_Hill3<-0
    for (i in 1:length(diff_linear_hill3)) {
      SSD_linear_Hill3<-diff_linear_hill3[i]+SSD_linear_Hill3
    }
    
    #log v.s Hill
    diff_log_hill3<-matrix(NA,length(log_hat3),1)
    for (i in 1:length(log_hat3)) {
      diff_log_hill3[i]<-(log_hat3[i]-hill_hat3[i])**2
    }
    #Sum of square of Difference between linear model and Hill model
    SSD_log_Hill3<-0
    for (i in 1:length(diff_log_hill3)) {
      SSD_log_Hill3<-diff_log_hill3[i]+SSD_log_Hill3
    }
    #Compare the linear model with log model
    #Can use the following ifelse() function as well.
    #ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Linear Model should suffice"), ifelse(SSD_log_Hill>SSD_linear_Hill, print("The Log Model should suffice"), print("Either Log model or lieanr model suffices")))
    
    
    
    #Compare the linear model with log model
    #identify the model for each drug-response relationship
    
    #1st drug
    if (SSD_log_Hill>SSD_linear_Hill) {
      print("The Linear Model should suffice for Drug 1")
      model_drug1="linear"
      drug1_a1<<-summary(linear)$coefficients[1,1]
      drug1_b1<<-summary(linear)$coefficients[2,1]
    } else if (SSD_log_Hill<SSD_linear_Hill) {
      print("The Log Model should suffice for Drug 1")
      model_drug1="log"
      drug1_a1<<-summary(log)$coefficients[1,1]
      drug1_b1<<-summary(log)$coefficients[2,1]
    } else {
      print("Either Log model or lieanr model suffices for Drug 1")
      model_drug1="linear"
    }
    
    #2nd drug  
    if (SSD_log_Hill2>SSD_linear_Hill2) {
      print("The Linear Model should suffice for Drug 2")
      model_drug2<-"linear"
      drug2_a2<<-summary(linear2)$coefficients[1,1]
      drug2_b2<<-summary(linear2)$coefficients[2,1]
    } else if (SSD_log_Hill2<SSD_linear_Hill2) {
      print("The Log Model should suffice for Drug 2")
      model_drug2<-"log"
      drug2_a2<<-summary(log2)$coefficients[1,1]
      drug2_b2<<-summary(log2)$coefficients[2,1]
    } else {
      print("Either Log model or lieanr model suffices for Drug 2")
      model_drug2<-"linear"
    }
    
    #3rd drug
    if (SSD_log_Hill3>SSD_linear_Hill3) {
      print("The Linear Model should suffice for Drug 3")
      model_drug3<-"linear"
      drug3_a3<<-summary(linear3)$coefficients[1,1]
      drug3_b3<<-summary(linear3)$coefficients[2,1]
    } else if (SSD_log_Hill3<SSD_linear_Hill3) {
      print("The Log Model should suffice for Drug 3")
      model_drug3<-"log"
      drug3_a3<<-summary(log3)$coefficients[1,1]
      drug3_b3<<-summary(log3)$coefficients[2,1]
    } else {
      print("Either Log model or lieanr model suffices for Drug 3")
      model_drug3<-"linear"
    }
    
    #Create new variables to be further used in the following if-else function
    Drug1<<-data
    Drug2<<-data2
    Drug3<<-data3
    
    
    # method: specifying the fitted models for the drug data;
    #"liloglog" 
    #"lilili" 
    #"lililog" 
    #"logloglog" 
    if (model_drug1=="linear" & model_drug2=="linear" & model_drug3=="linear") {
      # method="lilili"
      # Drug1: y=a1+b1*x, linear response
      # Drug2: y=a2+b2*x, linear response
      # Drug3: y=a3+b3*x, linear response
      method<<-"lilili"
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ Drug2[,1] )$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ Drug3[,1] )$coeff )
      a1<<-coef1[1]
      b1<<-coef1[2]
      a2<<-coef2[1]
      b2<<-coef2[2]
      a3<<-coef3[1]
      b3<<-coef3[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="linear" & model_drug2=="log" & model_drug3=="log") {
      # method="liloglog"
      # Drug1: y=a1+b1*x, linear response
      # Drug2: y=a2+b2*log(x), log response
      # Drug3: y=a3+b3*log(x), log response
      method<<-"liloglog"
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1] ))$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ log(Drug3[,1] ))$coeff )
      a1<<-coef1[1]
      b1<<-coef1[2]
      a2<<-coef2[1]
      b2<<-coef2[2]
      a3<<-coef3[1]
      b3<<-coef3[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="log" & model_drug2=="linear" & model_drug3=="log") {
      # method="loglilog"
      # Drug1: y=a1+b1*log(x), log response
      # Drug2: y=a2+b2*x, linear response
      # Drug3: y=a3+b3*log(x), log response
      method<<-"loglilog"
      coef1 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1] ))$coeff )
      coef2 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ log(Drug3[,1] ))$coeff )
      a1<<-coef2[1]
      b1<<-coef2[2]
      a2<<-coef1[1]
      b2<<-coef1[2]
      a3<<-coef3[1]
      b3<<-coef3[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="log" & model_drug2=="log" & model_drug3=="linear") {
      # method="loglogli"
      # Drug1: y=a1+b1*log(x), log response
      # Drug2: y=a2+b2*log(x), log response
      # Drug3: y=a3+b3*x, linear response
      method<<-"loglogli"
      coef1 <- as.numeric( lm( Drug1[,2] ~ log(Drug1[,1] ))$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1] ))$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ Drug3[,1] )$coeff )
      a1<<-coef3[1]
      b1<<-coef3[2]
      a2<<-coef1[1]
      b2<<-coef1[2]
      a3<<-coef2[1]
      b3<<-coef2[2]
      print(paste0("Method should be", " ", method))
    }  else if (model_drug1=="linear" & model_drug2=="linear" & model_drug3=="log") {
      # method="lililog"
      # Drug1: y=a1+b1*x, linear response
      # Drug2: y=a2+b2*x, linear response
      # Drug3: y=a3+b3*log(x), log response
      method<<-"lililog"
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ Drug2[,1] )$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ log(Drug3[,1] ))$coeff )
      a1<<-coef1[1]
      b1<<-coef1[2]
      a2<<-coef2[1]
      b2<<-coef2[2]
      a3<<-coef3[1]
      b3<<-coef3[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="linear" & model_drug2=="log" & model_drug3=="linear") {
      # method="lilogli"
      # Drug1: y=a1+b1*x, linear response
      # Drug2: y=a2+b2*log(x), log response
      # Drug3: y=a3+b3*x, linear response
      method<<-"lilogli"
      coef1 <- as.numeric( lm( Drug1[,2] ~ Drug1[,1] )$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1] ))$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ Drug3[,1] )$coeff )
      a1<<-coef1[1]
      b1<<-coef1[2]
      a2<<-coef3[1]
      b2<<-coef3[2]
      a3<<-coef2[1]
      b3<<-coef2[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="log" & model_drug2=="linear" & model_drug3=="linear") {
      # method="loglili"
      # Drug1: y=a1+b1*log(x), log response
      # Drug2: y=a2+b2*x, linear response
      # Drug3: y=a3+b3*x, linear response
      method<<-"loglili"
      coef1 <- as.numeric( lm( Drug1[,2] ~ log(Drug1[,1] ))$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ Drug2[,1] )$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ Drug3[,1] )$coeff )
      a1<<-coef2[1]
      b1<<-coef2[2]
      a2<<-coef3[1]
      b2<<-coef3[2]
      a3<<-coef1[1]
      b3<<-coef1[2]
      print(paste0("Method should be", " ", method))
    } else if (model_drug1=="log" & model_drug2=="log" & model_drug3=="log") {
      # method="logloglog"
      # Drug1: y=a1+b1*log(x), log response
      # Drug2: y=a2+b2*log(x), log response
      # Drug3: y=a3+b3*log(x), log response
      method<<-"logloglog"
      coef1 <- as.numeric( lm( Drug1[,2] ~ log(Drug1[,1] ))$coeff )
      coef2 <- as.numeric( lm( Drug2[,2] ~ log(Drug2[,1] ))$coeff )
      coef3 <- as.numeric( lm( Drug3[,2] ~ log(Drug3[,1] ))$coeff )
      a1<<-coef1[1]
      b1<<-coef1[2]
      a2<<-coef2[1]
      b2<<-coef2[2]
      a3<<-coef3[1]
      b3<<-coef3[2]
      print(paste0("Method should be", " ", method))
    } 
    
    
    #Create variables to be further used in SYN.design, SYN.test and SYN.analysis programs
    
    method<<-method
    
    #stop("Please enter the correct location of your dataset on the disk. Note:Use left slash sign after every folder names.")
    
    #visual check of the fitness: Drug1
    par(mfrow=c(1,3))
    fit.ll <- drm(Viability~Dose, data=data, fct=LL.5(), type="continuous")
    plot(fit.ll,type = "none",col='red',xlab = "Dose",ylab = "Viability(%)",main = "Drug 1")
    points(x=data$Dose,y=data$Viability)
    #Can use ggplot() to achieve a more smooth curve on the scattered plot
    
    #visual check of the fitness: Drug2
    fit.ll2 <- drm(Viability~Dose, data=data2, fct=LL.5(), type="continuous")
    plot(fit.ll2,type = "none",col='red',xlab = "Dose",ylab = "Viability(%)",main = "Drug 2")
    points(x=data2$Dose,y=data2$Viability)
    
    #visual check of the fitness: Drug3
    #visual check of the fitness: Drug3
    fit.ll3 <- drm(Viability~Dose, data=data3, fct=LL.5(), type="continuous")
    plot(fit.ll3,type = "none",col='red',xlab = "Dose",ylab = "Viability(%)",main = "Drug 3")
    points(x=data3$Dose,y=data3$Viability)
    
  }
  
  ThreeDrugIdentify<-ThreeDrugIdentify(location1,location2,location3)
  
  #Try this ThreeDrugIdentify function with the following example:
  #ThreeDrugIdentify<-ThreeDrugIdentify(location1 = "ARAC03.txt",location2 = "SAHA03.txt",location3 = "VP1603.txt")
  
  SY.design3<-function (c0, r,method)
  {
    if ( method != "liloglog" & method != "loglilog" & method != "loglogli" &  method != "lilili" & method != "loglili" & method != "lilogli" & method != "lililog" & method != "logloglog") { 
      stop( "method = unknown regression model." ) 
    }
    # function name: SY.design3(c0, r,method="lilili",UDFactor3=location)
    # data: from the single-dose experiments of three drugs
    # cited: UDFactor3
    # output: dose mixtures for combination experiments
    # parameters: 
    # a1 and b1: intercept and slope of drug1
    # a2 and b2: intercept and slope of drug2
    # a3 and b3: intercept and slope of drug3
    ##### b3 <= b2 <= b1 ????
    # si2: variance estimated from the pooled data of the single experiments
    # m: number of dose-mixtures for combination study
    # r: number of replications at each dose-mixture
    # c0: the smallest meaningful difference to be detected
    # totla sample size: m*r
    # power: 80%
    # significance level: 5%
    # Dis: the central L2-discrepancy of m mixtures
    # h1: to choose the low dose boundary from drug1: eg. 80%
    # h2: to choose the high dose boundary from drug1: eg. 20%
    ##--- input parameters of the single dose-responses
    si2 <- var(c(Drug1[, 2], Drug2[, 2], Drug3[,2]))   
    a1 <- a1
    b1 <- b1
    a2 <- a2
    b2 <- b2
    a3 <- a3
    b3 <- b3
    #si2 <- 988.422
    #a1 <- 4.8
    #b1 <- -12.76
    #a2 <- 41.52
    #b2 <- -13.02
    #a3 <- 54.55
    #b3 <- -23.98
    h1<-80
    h2<-20
    rho0 <- exp((a2 - a1)/b1)
    rho1 <- exp((a3 - a1)/b1)
    ##--- calculate the sample size for uniform design
    #m is the number of mixtures
    #k is the number of drugs, in this two drug combination case, k=3
    #3<=r<=8
    #3<=m<=30
    #k is automatically set to 3 since this is a three-drug combination study
    #d <- c0^2/si2
    #delta0: non-central parameter of the F-distribution
    #delta0 <- (m * (r+1)) * d
    #To obtain m, Step 1: Obtain the quantile under F-distribution
    #Note: alpha (significance level)=0.05
    #Substitute quantile with the following equation, otherwise an error occurs since m is the unknown parameter
    #and has to be estimated
    #quantile<<-qf(0.95,(m-2),(m*r))
    #Step 2: Calculate the non-central parameter as delta0
    #Substitute delta0 with the following equation, otherwise an error occurs since m is the unknown parameter
    #and has to be estimated
    #Since d <<- c0^2/si2 and delta0 <<- (m * (r+1)) * d, then delta0 can be written into:
    #delta0 <<- (m * (r+1)) * (c0^2/si2)
    #Step 3: Use pre-determined power (80%) to determine m
    #Constructing essential elements to be used in 'ktsolve' function
    r<<-r
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
    ##--- find the U-type matrix
    UDFactor3<-structure(.Data = list(c(0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
                                        17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30.)
                                      , c(4., 3., 1., 4., 2., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.039476, 1., 4., 3., 2., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(5., 2., 5., 1., 4., 3., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.026331, 4., 2., 1., 5., 3., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(6., 2., 4., 6., 1., 3., 5., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.018637, 3., 6., 2., 5., 1., 4., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(7., 5., 2., 7., 3., 6., 1., 4., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.01425, 4., 2., 6., 7., 1., 5., 3., NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(8., 3., 7., 5., 1., 8., 4., 2., 6., NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.0105, 4., 7., 1., 6., 3., 8., 2., 5., NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(9., 6., 2., 9., 3., 5., 7., 1., 8., 4., NA, NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.008758, 3., 8., 6., 1., 5., 9., 4., 2., 7., NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(10., 5., 8., 2., 10., 3., 7., 1., 9., 6., 4., NA, NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.007398, 4., 8., 2., 6., 10., 1., 7., 3., 9., 5., NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(11., 7., 3., 9., 4., 10., 1., 6., 11., 5., 2., 8., NA, NA, NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.006193, 6., 3., 10., 8., 2., 5., 11., 7., 1., 9., 4., NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(12., 6., 10., 2., 8., 4., 12., 3., 9., 7., 1., 11., 5., NA, NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.005332, 4., 8., 11., 1., 7., 10., 2., 6., 12., 5., 3., 9., NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(13., 8., 3., 12., 6., 10., 1., 4., 13., 9., 5., 2., 11., 7., NA, NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.004629, 7., 11., 4., 2., 13., 8., 5., 10., 1., 12., 3., 6., 9., NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(14., 7., 12., 4., 10., 1., 13., 5., 8., 2., 14., 9., 3., 11., 6., NA, NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.004064, 5., 9., 13., 1., 7., 12., 3., 8., 11., 4., 14., 2., 6., 10., NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(15., 5., 14., 10., 3., 7., 12., 2., 9., 15., 6., 1., 13., 8., 11., 4., NA,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.00359, 9., 5., 12., 2., 15., 7., 13., 1., 10., 4., 6., 14., 8., 3., 11.,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(16., 12., 6., 4., 14., 9., 1., 16., 8., 5., 11., 2., 15., 7., 13., 3., 10.,
                                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.003194, 7., 12., 4., 15., 1., 10., 5., 16., 8., 13., 2., 11., 6., 3.,
                                          14., 9., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(17., 6., 13., 11., 4., 16., 2., 8., 14., 10., 1., 17., 7., 5., 15., 9.,
                                          3., 12., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.002846, 7., 13., 3., 16., 10., 5., 11., 1., 17., 9., 6., 14., 2., 15.,
                                          4., 12., 8., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(18., 13., 3., 7., 18., 11., 5., 15., 9., 1., 16., 10., 4., 12., 6., 17.,
                                          2., 14., 8., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.002558, 8., 15., 2., 11., 18., 6., 3., 13., 9., 16., 5., 12., 1., 17.,
                                          7., 4., 14., 10., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(19., 12., 5., 17., 3., 15., 10., 7., 19., 1., 8., 11., 14., 4., 18., 9.,
                                          2., 16., 6., 13., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.002313, 6., 14., 11., 2., 18., 9., 16., 4., 7., 12., 19., 1., 10., 15.,
                                          3., 17., 8., 5., 13., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(20., 9., 17., 3., 15., 7., 12., 2., 20., 5., 13., 19., 6., 10., 14., 1.,
                                          18., 8., 11., 16., 4., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.002109, 15., 8., 4., 19., 11., 1., 17., 12., 6., 14., 3., 20., 7., 10.,
                                          9., 16., 2., 18., 5., 13., NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(21., 6., 16., 12., 1., 20., 9., 18., 4., 11., 15., 3., 21., 8., 14., 10.,
                                          5., 19., 2., 17., 7., 13., NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.001939, 9., 18., 2., 14., 5., 21., 11., 4., 16., 7., 19., 13., 1., 15.,
                                          8., 12., 20., 6., 3., 17., 10., NA, NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(22., 11., 4., 18., 14., 7., 21., 1., 15., 9., 17., 3., 22., 6., 13., 10.,
                                          19., 5., 16., 2., 20., 8., 12., NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.001785, 14., 8., 18., 3., 21., 6., 12., 16., 1., 9., 19., 13., 5., 22.,
                                          10., 2., 17., 11., 4., 20., 15., 7., NA, NA, NA, NA, NA, NA, NA, NA)
                                      , c(23., 14., 8., 21., 3., 17., 10., 6., 19., 13., 1., 23., 5., 16., 12., 9.,
                                          20., 2., 22., 11., 4., 18., 15., 7., NA, NA, NA, NA, NA, NA, NA)
                                      , c(0.001648, 12., 20., 4., 8., 17., 14., 2., 22., 6., 18., 10., 15., 1., 23.,
                                          9., 13., 5., 19., 3., 21., 7., 16., 11., NA, NA, NA, NA, NA, NA, NA)
                                      , c(24., 15., 6., 22., 4., 19., 12., 9., 24., 2., 13., 17., 7., 21., 1., 10.,
                                          18., 16., 5., 11., 23., 3., 14., 20., 8., NA, NA, NA, NA, NA, NA)
                                      , c(0.001527, 13., 19., 4., 8., 23., 16., 2., 11., 21., 5., 18., 6., 15., 12.,
                                          24., 1., 9., 17., 10., 20., 3., 22., 7., 14., NA, NA, NA, NA, NA,
                                          NA)
                                      , c(25., 14., 6., 22., 11., 18., 3., 20., 9., 24., 5., 16., 2., 13., 25., 8.,
                                          12., 21., 1., 17., 10., 23., 7., 15., 4., 19., NA, NA, NA, NA, NA)
                                      , c(0.001421, 16., 5., 23., 12., 2., 21., 10., 19., 7., 14., 25., 3., 8., 18.,
                                          11., 22., 4., 17., 15., 1., 13., 24., 6., 9., 20., NA, NA, NA, NA,
                                          NA)
                                      , c(26., 14., 7., 25., 4., 20., 17., 11., 1., 23., 9., 18., 13., 3., 26., 15.,
                                          6., 22., 8., 21., 12., 2., 16., 24., 10., 5., 19., NA, NA, NA, NA)
                                      , c(0.001321, 10., 22., 14., 3., 25., 7., 18., 16., 5., 12., 20., 1., 24., 9.,
                                          21., 6., 17., 13., 2., 26., 8., 15., 23., 4., 19., 11., NA, NA, NA,
                                          NA)
                                      , c(27., 17., 8., 25., 4., 20., 10., 14., 22., 1., 12., 27., 7., 18., 3., 24.,
                                          16., 9., 21., 5., 11., 26., 13., 2., 23., 15., 6., 19., NA, NA, NA)
                                      , c(0.001215, 20., 5., 13., 24., 2., 16., 9., 26., 11., 22., 7., 18., 15., 3.,
                                          23., 6., 27., 10., 14., 1., 19., 12., 21., 4., 25., 8., 17., NA, NA,
                                          NA)
                                      , c(28., 22., 4., 16., 11., 26., 8., 18., 1., 20., 13., 24., 6., 14., 28., 3.,
                                          21., 10., 7., 25., 15., 19., 2., 27., 9., 17., 5., 23., 12., NA, NA)
                                      , c(0.001157, 19., 6., 24., 12., 3., 27., 9., 17., 15., 8., 26., 2., 21., 11.,
                                          23., 5., 20., 14., 16., 1., 28., 10., 22., 4., 13., 25., 7., 18.,
                                          NA, NA)
                                      , c(29., 15., 24., 3., 9., 19., 28., 7., 13., 22., 2., 16., 27., 10., 18., 5.,
                                          26., 11., 21., 6., 20., 1., 25., 14., 12., 29., 8., 23., 4., 17.,
                                          NA)
                                      , c(0.001092, 21., 11., 8., 28., 4., 25., 15., 1., 18., 23., 13., 6., 17., 29.,
                                          3., 20., 10., 9., 26., 22., 12., 2., 16., 24., 14., 5., 27., 19.,
                                          7., NA)
                                      , c(30., 18., 5., 24., 12., 28., 7., 20., 2., 15., 25., 13., 10., 30., 4., 21.,
                                          8., 27., 16., 1., 22., 9., 17., 29., 14., 3., 19., 26., 6., 23., 11.)
                                      , c(0.001026, 9., 27., 19., 4., 13., 16., 19., 6., 21., 2., 23., 14., 25., 11.,
                                          7., 30., 17., 1., 18., 22., 8., 26., 5., 12., 24., 15., 28., 3., 10.,
                                          20.)
    )
    , names = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
                "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22",
                "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30", "C31", "C32",
                "C33", "C34", "C35", "C36", "C37", "C38", "C39", "C40", "C41", "C42",
                "C43", "C44", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52",
                "C53", "C54", "C55")
    , row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                    "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
                    "26", "27", "28", "29", "30", "31")
    , class = "data.frame"
    )
    
    UD <- UDFactor3
    Dis <- UD[1, (2 * m - 5)]
    MU <- matrix(c(UD[2:(m + 1), 1], UD[2:(m + 1), (2 * m - 6)], UD[2:
                                                                      (m + 1), (2 * m - 5)]), ncol = 3, byrow = F)
    
    if ( method == "liloglog" | method == "loglilog" | method == "loglogli" |  method == "lilili" | method == "loglili" | method == "lilogli" | method == "lililog") {
      # the total dose is from d01 to d02
      h01 <- h1 - a1
      h02 <- h2 - a1
      d01 <- (h1 - a1)/b1
      d02 <- (h2 - a1)/b1
      h001 <- h01/b1
      h002 <- h02/b1
      ##--- find the m dose-mixtures
      u1 <- (as.numeric(MU[, 1]) - 0.5)/m
      u2 <- (as.numeric(MU[, 2]) - 0.5)/m
      u3 <- (as.numeric(MU[, 3]) - 0.5)/m
      y1 <- 1 - sqrt(1 - u1)
      y2 <- (1 - y1) * u2
      y03 <- u3 * (h01^3 - h02^3) + h02^3
      y3 <- ((abs(y03))^(1/3) * sign(y03) - h02)/(h01 - h02)
      A <- matrix(c(h002, 0, 0, 0, (h02/b2), 0, 0, 0, (h02/b3), h001, 0,
                    0, 0, (h01/b2), 0, 0, 0, (h01/b3)), ncol = m , byrow = F)
      B <- matrix(c((1 - y1 - y2) * (1 - y3), y1 * (1 - y3), y2 * (1 - y3),
                    (1 - y1 - y2) * y3, y1 * y3, y2 * y3), nrow = m, byrow = T)
      Z <- A%*% B
      z1 <- c(Z[1,  ])
      z2 <- c(Z[2,  ])
      z3 <- c(Z[3,  ])
      ##--- d031 <- (h1 - a2)/b2
      ##--- d032 <- (h2 - a2)/b2
      ##--- d041 <- (exp((h1 - a3)/b3))
      ##--- d042 <- (exp((h2 - a3)/b3))
      if ( method == "lilili"){
        # method="lilili"
        # Drug1: y=a1+b1*x, linear response
        # Drug2: y=a2+b2*x, linear response
        # Drug3: y=a3+b3*x, linear response
        d1 <- abs(z1)
        d2 <- abs(z2 + (a1 - a2)/b2)
        d3 <- abs(z3 + (a1 - a3)/b3)
        method <- "lilili"
      }
      
      if ( method == "lililog"){
        # method="lililog"
        # Drug1: y=a1+b1*x, linear response
        # Drug2: y=a2+b2*x, linear response
        # Drug3: y=a3+b3*log(x), log-linear response
        z01 <- exp((b1 * z1 + b2 * z2 + z3 + a1 - a3)/b3)
        d1 <- abs(z1)
        d2 <- abs(z2 + (a1 - a2)/b2)
        d3 <- abs((z3 * z01)/(b1 * z1 + b2 * z2 + z3))
        method <-"lililog"
      }
      
      if ( method == "lilogli"){
        # method="lilogli"
        # Drug1: y=a1+b1*x, linear response
        # Drug2: y=a2+b2*log(x), log-linear response
        # Drug3: y=a3+b3*x, linear response
        z01 <- exp((b1 * z1 + b3 * z3 + z2 + a1 - a2)/b2)
        d1 <- abs(z1)
        d2 <- abs((z2 * z01)/(b1 * z1 + b3 * z3 + z2))
        d3 <- abs(z3 + (a1 - a3)/b3)
        method <-"lilogli"
      }
      
      if ( method == "loglili"){
        # method="loglili"
        # Drug1: y=a1+b1*log(x), log-linear response
        # Drug2: y=a2+b2*x, linear response
        # Drug3: y=a3+b3*x, linear response
        z01 <- exp((b2 * z2 + b3 * z3 + z1 + a2 - a1)/b1)
        d1 <- abs((z1 * z01)/(b2 * z2 + b3 * z3 + z1))
        d2 <- abs(z2)
        d3 <- abs(z3 + (a1 - a3)/b3)
        method <-"loglili"
      }
      
      if ( method == "liloglog"){
        # method="liloglog"
        # Drug1: y=a1+b1*x, linear response
        # Drug2: y=a2+b2*log(x), log-linear response
        # Drug3: y=a3+b3*log(x), log-linear response
        z02 <- exp((a2-a3)/b2)
        z03 <-(b1 * z1 + z2 + z3 + a1 - a3)/b3
        d1 <- abs(z1)
        d2 <- abs((z2*z02*exp((z03)^(b3/b2)))/(b1 * z1 + z2 + z3))
        d3 <- abs((z3*exp(z03))/(b1 * z1 + z2 + z3))
        method <- "liloglog"
      }
      
      if ( method == "loglilog"){
        # method="loglilog"
        # Drug1: y=a1+b1*log(x), log-linear response
        # Drug2: y=a2+b2*x, linear response
        # Drug3: y=a3+b3*log(x), log-linear response
        z02 <- exp((a1-a3)/b1)
        z03 <-(b2 * z2 + z1 + z3 + a2 - a3)/b3
        d1 <- abs((z1*z02*exp((z03)^(b3/b1)))/(b2 * z2 + z1 + z3))
        d2 <- abs(z2)
        d3 <- abs((z3*exp(z03))/(b2 * z2 + z1 + z3))
        method <- "loglilog"
      }
      
      if ( method == "loglogli"){
        # method="loglogli"
        # Drug1: y=a1+b1*log(x), log-linear response
        # Drug2: y=a2+b2*log(x), log-linear response
        # Drug3: y=a3+b3*x, linear response
        z02 <- exp((a1-a2)/b1)
        z03 <-(b3 * z3 + z1 + z2 + a3 - a2)/b2
        d1 <- abs((z1*z02*exp((z03)^(b2/b1)))/(b3 * z3 + z1 + z2))
        d2 <- abs((z2*exp(z03))/(b3 * z3 + z1 + z2))
        d3 <- abs(z3)
        method <- "loglogli"
      }
      
      Dose.mixtures <- matrix(c(d1, d2, d3), ncol = 3, byrow = F)
      Mixture.number <- m+1
      Discrepancy <- Dis
      Experiment.number <- (m+1) * r
      Variables <- c("Drug1", "Drug2", "Drug3")
      dimnames(Dose.mixtures) <- list(NULL, Variables)
      return(list(Mixture.number, Experiment.number, Discrepancy, Dose.mixtures,method))
    }
    
    if ( method == "logloglog"){
      # method="logloglog"
      # Drug1: y=a1+b1*log(x), log-linear response
      # Drug2: y=a2+b2*log(x), log-linear response
      # Drug3: y=a3+b3*log(x), log-linear response
      ##--- choose the experimental domain
      # the total dose is from d01 to d02
      h1 <- 80
      h2 <- 20
      d01 <- exp((h1 - a1)/b1)
      d02 <- exp((h2 - a1)/b1)
      d03 <- (exp((h2 - a2)/b2))
      d04 <- (exp((h2 - a3)/b3))
      ##--- find the m dose-mixtures
      y1 <- (as.numeric(MU[, 1]) - 0.5)/m
      y2 <- (as.numeric(MU[, 2]) - 0.5)/m
      y3 <- (as.numeric(MU[, 3]) - 0.5)/m
      z1 <- y1 * (d02 - d01) + d01
      z2 <- y2 * sqrt(y3)
      z3 <- (1 - y2) * sqrt(y3)
      d1 <- z1 * z2 * (1 - (1 - rho1/rho0) * z3)
      d2 <- rep(0, m)
      d3 <- rep(0, m)
      D2 <- seq(0.001, d03, 0.002)
      D3 <- seq(0.001, d04, 0.002)
      r1 <- length(D2)
      r2 <- length(D3)
      D4 <- c(kronecker(D2, rep(1, r2)))
      D5 <- c(kronecker(rep(1, r1), D3))
      D0 <- rep(0, m)
      dd2 <- rep(0, m)
      dd3 <- rep(0, m)
      w0 <- (rho0/rho1)^(b1/b2)
      w1 <- (b3 * (b2 - b1))/(b2 * (b3 - b1))
      w2 <- b1/b3
      #return(w1,w0,rho0,rho1)
      f1 <- w0 * w1 * D4 + D5
      for(j in 1:m) {
        f2 <- d1[j]/rho1 + w0 * (1 - w1) * D4
        f3 <- ((1 - w2) * f2 * f1 - f2^(1 + w2) + f2 * sqrt(((1 - w2) *
                                                               f1 - f2^w2)^2 + 2 * w2 * (1 - w2) * f1^2))/(w2 * (
                                                                 1 - w2) * f1^2)
        f4 <- (z1[j] * (1 - z2[j] - z3[j]) + (1 - rho1/rho0) * z1[
          j] * z2[j] * z3[j] - (rho0/rho1)^(b1/b2 - 1) * D4 *
            f3^w1)^2 + (z1[j] * z3[j] - D5 * f3^w1)^2
        Ind <- sort.list(f4)[1]
        d2[j] <- (D4[Ind])^(b1/b2) * exp((a1 - a2)/b2)
        d3[j] <- (D5[Ind])^(b1/b3) * exp((a1 - a3)/b3)
        D0[j] <- f4[Ind]
        dd2[j] <- D4[Ind]
        dd3[j] <- D5[Ind]
        method <-"logloglog"
      }
      Dose.mixtures <- matrix(c(d1, d2, d3), ncol = 3, byrow = F)
      Mixture.number <- m+1
      Discrepancy <- Dis
      Experiment.number <- (m+1) * r
      Variables <- c("Drug1", "Drug2", "Drug3")
      dimnames(Dose.mixtures) <- list(NULL, Variables)
      return(list(Mixture.number, Experiment.number, round(Discrepancy,3), round(Dose.mixtures,3),method,
                  D0, dd2, dd3))
    }
  }
  
  # Try in this example:
  #SYN.design3 <- SY.design3(0.01, 6,method="logloglog")
  SYN.design3<-SY.design3(c0=c0, r=r,method)
  
  print(paste0("The method used is"," ",SYN.design3[[5]]))
  print(paste0("The number of combinations (mixtures) needed is m="," ", SYN.design3[[1]]))
  print(paste0("Total sample size at least is m*r="," ", SYN.design3[[2]]))
  print(paste0("The central L2-discrepancy of m mixtures is"," ", SYN.design3[[3]]))
  print(paste0("Please use the following scheme of dose mixtures to measure the response, input the results into a txt file and name as 'Combination' "))
  print(SYN.design3[[4]])
  print("The Combination dataset should have the pattern as: column1--dose of drug1, column 2--dose of drug2, column 3--dose of drug3,column 4--response")   
}
  

  
  