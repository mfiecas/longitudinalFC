#load necessary libraries
library(RcppArmadillo)
library(Rcpp)
library(MTS)
library(Matrix)
library(MASS)
library(matrixcalc)
library(statmod)

#bring in Rcpp functions
sourceCpp('C:/Users/Brian/Google Drive/PhD Work/LongitudinalfMRI/LongitInclude.cpp')
sourceCpp('C:/Users/Brian/Desktop/LongitInclude.cpp')


##########################   Simulate cross-sectional data ##############
###Set all parameters for simulated data to be analyzed
main=0.1      #baseline group difference in first ROI pair
N=20          #Number of subjects per group
P=3           #Number of ROIs

T=120         #length of time series
autocor=0.3   #autocorrelation in time series

N1=N          #Number of subjects in group 1
N2=N          #Number of subjects in group 2
Q=choose(P,2)         #Number of ROI pairs
phi=diag(autocor,P) 

#Define true variance structure
sigma2G1=0.01     #Diagonal elements of true error variance
kappaG1=0.001     #Off diagonal elements of diagonal blocks of true error variance

PSI0G1=matrix(kappaG1,Q,Q)-diag(kappaG1,Q)+diag(sigma2G1,Q)
PSI0G2=PSI0G1

mu1=rep(0,Q)    #Mean vector for simulated true FC for group 1
mu2=c(main,rep(0,Q-1)) #Mean vector for simulated true FC for group 2

sig1=mvrnorm(n=N1, mu=mu1, Sigma=PSI0G1) #simulate group 1 FC
sig2=mvrnorm(n=N2, mu=mu2, Sigma=PSI0G2) #simulate group 2 FC

DataList=list()       #initiate list for data

#simulate group 1 time series
for(i in 1:N1){
  sigmat = diag(1,P)
  sigmat[lower.tri(sigmat)] = sig1[i,]
  sigmat[upper.tri(sigmat)] <- t(sigmat)[upper.tri(sigmat)]
  while(!is.positive.definite(sigmat)){sigmat=sigmat+diag(0.05,P)}
  DataList[[i]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat, theta=matrix(0,P,P))$series
}
for(i in 1:N2){
  sigmat = diag(1,P)
  sigmat[lower.tri(sigmat)] = sig2[i,]
  sigmat[upper.tri(sigmat)] <- t(sigmat)[upper.tri(sigmat)]
  while(!is.positive.definite(sigmat)){sigmat=sigmat+diag(0.05,P)}
  DataList[[N1+i]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat, theta=matrix(0,P,P))$series
}


###############   Sample code to fit cross-sectional model   ###############################

#datalist is a list where each element is a TxP matrix of data for one subject
#N1 is the number of subjects in group 1
#N2 is the number of subjects in group 2
#Nperms is the number of permutations to run for the permutation test
#lag is the maximum number of lags to be used in calculations
#bw is the bandwidth to be used for the MB windowing function
#step2_tol is the stoping criteria for iteratively solve for beta and Psi
#Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
  #"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
#MaxIter is maximum iterations to perform solving for Psi and beta
#verbose is boolean to signify if permutation test status updates should be printed to the console
#SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"

NPerm=1000
CSFit=FCanalysis(datalist=DataList, N1=N1, N2=N2, Nperms=NPerm)

CSFit$T_global
permp((1-ecdf(CSFit$T_dist_global)(CSFit$T_global))*NPerm,NPerm,N1,N2)

#Calculate p values with correction from Phipson 2010
PLocal=rep(NA,Q)  #initialize vector of pvalues for the local tests
for(q in 1:Q){
  PLocal[q]=permp((1-ecdf(CSFit$T_dist_local[,q])(CSFit$T_local[q]))*NPerm,NPerm,N1,N2)
}

#apply false discovery rate correction, shows Q interaction effect p-values then Q main effect p-values
Padj=round(p.adjust(PLocal,method="BH"),4)


