# ComStat
Version 1.beta, Oct 8, 2018
Statistical Design and Analysis of Synergy in Combination Studies 
Ming T. Tan, Hong-Bin Fang, with assistance of Robert Gu
Georgetown University Medical Center

Contact: Dr. Ming Tan (mtt34@georgetown.edu ) or Dr. Hong-Bin Fang (hf183@georgetown.edu).
For comments/questions/corrections regarding this Readme document and the ComStat Tutorial, contact Robert Gu (zg83@georgetown.edu) 

This R package program may be freely copied and noncommercially distributed as long as the source is acknowledged.

GETTING STARTED: 

Installation of the Software: 
1.	You are looking at the “Readme” file included with the ComStat software. That means you probably have already opened the ComStat.zip file from the github. If not please go to: http://www.umgcc.org/research/biostat_software.htm and click on “download COMSTAT.” 
2.	Unzip the files in the ComStat package and save them to a directory on your computer (eg., “My Documents” if using Windows). You will have 8 files to unzip and save to the directory: 
            One adobe file (Readme.pdf – this document), and one text file 
            Three R program files (SYN.design.R, SYN.test.R, SYN.analysis.R) 
            Three Excel files (contain sample data sets for training purposes)  
3.	Install the R software environment. To use ComStat, you will need to have the R computing language/environment installed on your computer. If you do not already have R on your computer, download R for free from the Internet. Go to http://lib.stat.cmu.edu/R/CRAN/ and download. For example, if you run Windows, click Windows, then on the subdirectory “base,” then on the file R-2.7.2-win32.exe to download the Windows version. The setup process will automatically start and you will end up with the R icon on your desktop. Double click the R icon to start the R program. Once in the program (R Graphics Users Interface, or RGui), under the Packages menu, go to a CRAN mirror site nearest you, and then in the same menu, Install packages “xlsReadWrite” and “fields.” The package “spam” will also install automatically. If you are unsure how to download R, or confused about its use, your friendly neighborhood statistician will probably be able to help. 
4.	While in R (RGui/RConsole), next engage the ComStat R function programs into the R system. First, under the RGui field, point to the directory where the ComStat files are, using File/Change dir…, then select the directory containing the files (eg., C:\Documents and Settings\\MyDocuments\200801COMSTAT). Then load the SysStat program files into R by using the following commands in the RConsole field after the red prompt (suggest you “copy-and-paste” each command below, excluding >); hit “Enter” after copying each command: 
> source("SYN.design.R") 
> source("SYN.test.R") 
> source("SYN.analysis.R") 

5.	You are now ready to input data files. We recommend that you use the sample data files (in Excel) in conjunction with the tutorial on page 6 to become familiar with how to use ComStat for the design and analysis of drug combination experiments. 


THEORETICAL BACKGROUND AND METHODS: 

This suite of computer programs provides the entire statistical design and analysis modules for two-drug combination studies. It consists of three R functions: SYN.design, SYN.test, SYN.analysis. It represents an efficient way to select doses in the combinations and calculate the sample sizes to produce data to characterize the joint action of two agents, as well as an efficient way to analyze the data so produced.

The dose-responses of individual drugs are assumed to be log-linear or linear, for a goodnees-of-fit of individual dose-responses, sometimes we can consider only the responses in a dose-range of interest, such as ED20 ~ ED80, or IC10 ~IC90, etc.

For a common Hill model (nonlinear sigmoid model) compared with log-linear model or linear model to fit the dose-response curves, the log-linear individual dose-response model is used for experimental design and interaction analysis in combination studies.


The experimental design (dose-finding and sample size determination) is derived by means of uniform measures that maximize the minimum power of the F-test to detect any departures from additive action, and at the same time minimize the maximum bias due to lack of fit among all potential departures of a given meaningful magnitude. A model-free interaction index surface is used to capture the interaction of two drugs. The nonparametric function of interaction index is estimated using the technique developed in thin plate splines. These methods are applicable to both in vivo and in vitro experiments. One example enclosed provides a guideline to use this program.


References:

1.	Fang HB, Ross DD, Sausville E, and Tan M (2008). Experimental design and interaction analysis of combination studies of drugs with log-linear dose response. Statistics in Medicine 27:3071-3083, 2008. 
2.	Fang HB, Tian GL, Li W, and Tan M (2009). Dose finding and sample size determination for evaluating combinations of drugs of linear and loglinear dose response curves. Journal of Biopharmaceutical Statistics 19, in press.
3.	Tan M, Fang HB,  Tian  GL, Houghton  PJ (2003).  Experimental design and sample size determination for testing synergy in drug combination studies based on uniform measures.  Statistics in Medicine 22: 2091-2100. 

***********************************************************************

Datasets: 
Drug1: dataset of Drug 1, column 1-- dose, column 2--response     
Drug2: dataset of Drug 2, column 1-- dose, column 2--response  
Combination: dataset of Combination response of drug 1 and 2; column1--dose of drug1, column 2--dose of drug2, column 3--response 

SYN.design program: Experimental design module
 
 Input: Drug1, Drug2, h, c0, r 
            data: single-dose experiments of two drugs
            Drug1: dataset of Drug 1, column 1-- dose, column 2--response     
            Drug2: dataset of Drug 2, column 1-- dose, column 2--response  
            h: dose range for experiments h; 
                for example, h=c(20, 80) means the dose range is from ED20 to ED80. 
            c0: the smallest meaningful difference to be detected; 
                  for example, c0(=15, default, the cell viability).
            r: number of repetitions plus one at each dose-mixture (r-1 replicates); 
               for example, r=6 means 5 replicates at each combination. (Note: r3).
plot.design: switching on/off the marginal model and combination scheme plots. The plots are displayed if Plot.Design = TRUE; they are not displayed if Plot.Design = FALSE.
	 method: specifying the fitted models for the drug data;
	         "loglog" for linear log models on both drugs,
	         "lili" for linear models on both drugs,
	         "lilog" for linear model on Drug1 and log model on Drug2,
	         "logli" for log model on Drug1 and linear model on Drug2.

Output: Mixture number (m): number of dose-mixtures for combination study
              Experiment number (m*r): totla sample size
              Discrepancy: the central L2-discrepancy of m mixtures
              Dose.mixtures: dose mixtures for combination experiments 
 
Default: data set "UniformDesign" is needed. 
               The number of combinations is less than 30.
               power: 80%; significance level: 5%

Upon completion of the combination experiments, the combination dataset is obtained.

SYN.test program: module to test for additive action 

Input: Drug1, Drug2, Combination, r
           Drug1: dataset of Drug 1, column 1-- dose, colummn 2--response     
           Drug2: dataset of Drug 2, column 1-- dose, colummn 2--response  
           Combination: dataset of Combination response of drug 1 and 2;
           column1--dose of drug1, column 2--dose of drug2, column 3--response 
           r: number of repetitions plus one at each dose-mixture (r-1 replicates)    
           method: specifying the fitted models for the drug data;
	        "loglog" for linear log models on both drugs,
	        "lili" for linear models on both drugs,
	        "lilog" for linear model on Drug1 and log model on Drug2
	        "logli" for log model on Drug1 and linear model on Drug2

Output: Degrees of freedom of the F test 
              F-statistic value
              P-value


If p-value 0.05, the result is that the two drugs are additive in this dose area. Otherwise (p<0.05),  use SYN.analysis to analyze the combination index.


SYN.analysis program:  module to obtain Interaction Index Surface 

Input: Drug1, Drug2, Combination, r, levels
           Drug1: dataset of Drug 1, column 1-- dose, column 2--response     
           Drug2: dataset of Drug 2, column 1-- dose, column 2--response  
           Combination: dataset of Combination response of drug 1 and 2;
           column1--dose of drug1, column 2--dose of drug2, column 3--response 
           r: number of repetitions plus one at each dose-mixture (r-1 replicates)
           method: specifying the fitted models for the drug data;
	         "loglog" for linear log models on both drugs,
	         "lili" for linear models on both drugs,
	         "logli" for log model on Drug1 and linear model on Drug2,
	         "lilog" for linear model on Drgu1 and log model on Drug2.
 levels: specifying the levels in the resulting contour plot. For the best effect of color display, the length of the vector assigned to levels should not exceed 8. By default, levels = c(0, 0.3, 0.6, 0.8, 1.3, 1.53, 1.76, 2).
contour.xlim: dose ranges to be displayed on the Drug1 axis. By default the program uses the minimum and maximum Drug1 doses found in the Combination data.
contour.ylim: dose ranges to be used on the Drug2 axis. By default the program uses the minimum and maximum Drug2 doses found in the Combination data.
plot.response: switching ON/OFF the plotting of the response surface. The response surface is plotted if plot.response = TRUE; no plotting is generated if plot.response = FALSE.
 plot.index: switching ON/OFF the plotting of the estimated interaction index surface. The index surface is plotted if plot.index = TRUE; no plotting is generated if plot.index = FALSE
Output: contour plot of the interaction index. The solid curve indicates exactly additive combinations (=1). The dotted curves indicate 95% confidence interval for the additive action.
index: the estimated indices on a  grids covers all the combination doses.
x1, x2: the coordinates of the grid points on which the variable index is defined.
**********************************************************************



TUTORIAL: How to use ComStat for the design and analysis of drug combination experiments. 

Introduction: ComStat not only analyzes combination data, but also helps you to design the combination experiment itself such that it provides optimally informative data. Let us say that you want to determine whether two drugs (drug 01 and drug 02) have synergistic, additive or antagonistic interactions in their effect when combined together. 

Step one: Once you have chosen two drugs to study in combination, first determine the relationship between dose and effect (eg., cytotoxicity, enzyme inhibition, etc.) for EACH drug used alone. Be sure to determine experimentally the doses of each drug that encompass the ED20 to ED80 dose ranges. Place these data in an Excel file exactly as given in the example files, below. Then use the SYN.design.R program to help you set up the combination experiments. Work through the example below to show you how to do this. 


Step Three. Analyzing the nature of drug interaction for non-additive combinations. To obtain a contour plot of values of Tau (see page 2 above) displayed for concentrations of drug 01 vs drug 02, and to display the dose-response surface for drug 01 vs drug 02.
