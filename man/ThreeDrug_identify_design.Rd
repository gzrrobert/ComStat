% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ThreeDrug1.R
\name{ThreeDrug_identify_design}
\alias{ThreeDrug_identify_design}
\title{Identifying the linear/log-linear model of each dose-response curve}
\usage{
ThreeDrug_identify_design(location1, location2, location3, c0, r)
}
\arguments{
\item{location1}{location of Drug 1's dataset on your hard disk}

\item{location2}{location of Drug 2's dataset on your hard disk}

\item{location3}{location of Drug 3's dataset on your hard disk}

\item{c0}{the smallest meaningful difference to be detected}

\item{r}{number of replications at each dose-mixture, and usually 3<=r<=8, default=3}
}
\value{
The fitted models and method for three drugs, 
        the number of dose-mixtures for combination, 
        total sample size,
        fitted plots for three drugs,
        dose mixtures for experimental design
}
\description{
This function utilized nls() function to model the HILL model of each relationship and substitute each relationship
with linear or log-linear relationship, it then designs a dose-mixture combination dataset for the study by using Uniform design.
}
\examples{
# Try this function with the following examples:
#(please be patient when the method is"logloglog"):
ThreeDrug_identify_design_1<- ThreeDrug_identify_design(
location1="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug1.txt", 
location2="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug2.txt", 
location3="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug3.txt", 
c0=13, r=6)
ThreeDrug_identify_design_2<- ThreeDrug_identify_design(
location1="C:/Users/58412/Desktop/SRA/Three drugs/data1/ARAC03.txt", 
location2="C:/Users/58412/Desktop/SRA/Three drugs/data1/SAHA03.txt", 
location3="C:/Users/58412/Desktop/SRA/Three drugs/data1/VP1603.txt", 
c0=15, r=6)
}
