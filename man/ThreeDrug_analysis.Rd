% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ThreeDrug3.R
\name{ThreeDrug_analysis}
\alias{ThreeDrug_analysis}
\title{Drawing contour plots of the index surface}
\usage{
ThreeDrug_analysis(location1, location2, location3, location4, method)
}
\arguments{
\item{location1}{location of Drug 1's dataset on your hard disk}

\item{location2}{location of Drug 2's dataset on your hard disk}

\item{location3}{location of Drug 3's dataset on your hard disk}

\item{location4}{location of combination dataset on your hard disk}

\item{method}{input the method you received in the previous step}
}
\value{
Contour plots showing the additive/synergistic/antagonistic effect
}
\description{
This function draws contour plots showing the additive/synergistic/antagonistic effect
}
\examples{
# Try this function with the following example:
SYN.analysis1 <-ThreeDrug_analysis(
location1="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug1.txt", 
location2="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug2.txt", 
location3="C:/Users/58412/Desktop/SRA/Three drugs/data2/Drug3.txt", 
location4 = "C:/Users/58412/Desktop/SRA/Three drugs/data2/Combination.txt",
method="lilili")
}
