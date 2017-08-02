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


#######################   Code to Simulate some data   ########################
#This code simulates up to 3 visits per subject

###Set all parameters for simulated data to be analyzed
main=0.1      #baseline group difference in first ROI pair
trend=0.1     #group difference in longitudinal trend in first ROI pair
N=20          #Number of subjects per group
P=3           #Number of ROIs
  
T=120         #length of time series
autocor=0.3   #autocorrelation in time series
v2prob=1      #probability each subject attends visit 2, set to 0 if you want cross-sectional data
v3prob=1      #probability each subject attends visit 3, set to 0 if you want cross-sectional data

N1=N          #Number of subjects in group 1
N2=N          #Number of subjects in group 2
Q=choose(P,2)         #Number of ROI pairs
phi=diag(autocor,P) 

#Define true variance structure
sigma2G1=0.01     #Diagonal elements of true error variance
kappaG1=0.001     #Off diagonal elements of diagonal blocks of true error variance
nuG1=0.005        #diagonal elements of true Psi1
omegaG1=0.001     #off diagonal elements of true Psi1
PSI0G1=matrix(kappaG1,Q,Q)-diag(kappaG1,Q)+diag(sigma2G1,Q)
PSI1G1=matrix(omegaG1,Q,Q)-diag(omegaG1,Q)+diag(nuG1,Q)
PSIG1=rbind(cbind(PSI0G1,PSI1G1,PSI1G1),
          cbind(PSI1G1,PSI0G1,PSI1G1),
          cbind(PSI1G1,PSI1G1,PSI0G1))
PSIG2=PSIG1

mu1=rep(0,Q*3)    #Mean vector for simulated true FC for group 1
mu2=c(main,rep(0,Q-1),main+trend,rep(0,Q-1),main+2*trend,rep(0,Q-1)) #Mean vector for simulated true FC for group 2

sig1=mvrnorm(n=N1, mu=mu1, Sigma=PSIG1) #simulate group 1 FC
sig2=mvrnorm(n=N2, mu=mu2, Sigma=PSIG2) #simulate group 2 FC
  
counter=1
time=matrix(,0,1)     #initiate time/visit matrix
visvec=rep(1,N1+N2)   #initiate visvec, vector with # of visits for each subject
DataList=list()       #initiate list for data

#simulate group 1 time series
for(i in 1:N1){
  sigmat11 = diag(1,P)
  sigmat11[lower.tri(sigmat11)] = sig1[i,1:Q]
  sigmat11[upper.tri(sigmat11)] <- t(sigmat11)[upper.tri(sigmat11)]
  while(!is.positive.definite(sigmat11)){sigmat11=sigmat11+diag(0.05,P)}
  DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat11, theta=matrix(0,P,P))$series
  time=rbind(time,c(0))
  counter=counter+1

  if(runif(1) < v2prob){
    sigmat12 = diag(1,P)
    sigmat12[lower.tri(sigmat12)] = sig1[i,(Q+1):(2*Q)]
    sigmat12[upper.tri(sigmat12)] <- t(sigmat12)[upper.tri(sigmat12)]
    while(!is.positive.definite(sigmat12)){sigmat12=sigmat12+diag(0.05,P)}
    DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat12, theta=matrix(0,P,P))$series
    time=rbind(time,c(1))
    visvec[i]=visvec[i]+1
    counter=counter+1
  }

  if(runif(1) < v3prob){
    sigmat13 = diag(1,P)
    sigmat13[lower.tri(sigmat13)] = sig1[i,(2*Q+1):(3*Q)]
    sigmat13[upper.tri(sigmat13)] <- t(sigmat13)[upper.tri(sigmat13)]
    while(!is.positive.definite(sigmat13)){sigmat13=sigmat13+diag(0.05,P)}
    DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat13, theta=matrix(0,P,P))$series
    time=rbind(time,c(2))
    visvec[i]=visvec[i]+1
    counter=counter+1
  }
}
  
#simulate group 2 time series
for(i in 1:N2){
  sigmat21 = diag(1,P)
  sigmat21[lower.tri(sigmat21)] = sig2[i,1:Q]
  sigmat21[upper.tri(sigmat21)] = t(sigmat21)[upper.tri(sigmat21)]
  while(!is.positive.definite(sigmat21)){sigmat21=sigmat21+diag(0.05,P)}
  DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat21, theta=matrix(0,P,P))$series
  time=rbind(time,c(0))
  counter=counter+1
  
  if(runif(1) < v2prob){
    sigmat22 = diag(1,P)
    sigmat22[lower.tri(sigmat22)] = sig2[i,(Q+1):(2*Q)]
    sigmat22[upper.tri(sigmat22)] = t(sigmat22)[upper.tri(sigmat22)]
    while(!is.positive.definite(sigmat22)){sigmat22=sigmat22+diag(0.05,P)}
    DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat22, theta=matrix(0,P,P))$series
    time=rbind(time,c(1))
    visvec[N1+i]=visvec[N1+i]+1
    counter=counter+1
  }
    
  if(runif(1) < v3prob){
    sigmat23 = diag(1,P)
    sigmat23[lower.tri(sigmat23)] = sig2[i,(2*Q+1):(3*Q)]
    sigmat23[upper.tri(sigmat23)] = t(sigmat23)[upper.tri(sigmat23)]
    while(!is.positive.definite(sigmat23)){sigmat23=sigmat23+diag(0.05,P)}
    DataList[[counter]]=VARMAsim(T, arlags=1, malags=0, phi=phi, sigma=sigmat23, theta=matrix(0,P,P))$series
    time=rbind(time,c(2))
    visvec[N1+i]=visvec[N1+i]+1
    counter=counter+1
  }
}



###############   Sample code to fit longitudinal model   ###############################
########## inputs for FCLongit function
#datalist is a list where each element is a TxP matrix of data for one subject
#N1 is the number of subjects in group 1
#N2 is the number of subjects in group 2
#Nperms is the number of permutations to run for the permutation test
#time is a vector where each entry is the time (e.g. age or visit number) variable for a visit
#visvec is vector with the number of visits for each subject
#lag is the maximum number of lags to be used in calculations
#bw is the bandwidth to be used for the MB windowing function
#step2_tol is the stoping criteria for iteratively solve for beta and Psi
#MaxIter is maximum iterations to perform solving for Psi and beta
#verbose is boolean to signify if permutation test status updates should be printed to the console
#Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
  #"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
#Psi1type denotes the structure to be used for Psi1, either "CS" for compound symmetry,
  #"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
#SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"

NPerm=1000
LongitFit=FCLongit(DataList, N1=N1, N2=N2, Nperms=NPerm, time=time, visvec=visvec)

#global hypothesis tests
LongitFit$TmainGlobal
permp((1-ecdf(LongitFit$TmainGlobal_dist)(LongitFit$TmainGlobal))*NPerm,NPerm,N1,N2)

LongitFit$TintGlobal
permp((1-ecdf(LongitFit$TintGlobal_dist)(LongitFit$TintGlobal))*NPerm,NPerm,N1,N2)


#Calculate p values with correction from Phipson 2010
PmainLocal=rep(NA,Q)  #initialize vector of pvalues for the local main effect tests
PintLocal=rep(NA,Q)   #initialize vector of pvalue for the local interaction tests
for(q in 1:Q){
  PmainLocal[q]=permp((1-ecdf(LongitFit$TmainLocal_dist[,q])(LongitFit$TmainLocal[q]))*NPerm,NPerm,N1,N2)
  PintLocal[q]=permp((1-ecdf(LongitFit$TintLocal_dist[,q])(LongitFit$TintLocal[q]))*NPerm,NPerm,N1,N2)
}

#apply false discovery rate correction, shows Q interaction effect p-values then Q main effect p-values
Padj=round(p.adjust(c(PintLocal,PmainLocal),method="BH"),4)




