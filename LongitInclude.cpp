#include <iostream>
#include <string>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


//////////////////////////////    Functions used in cross sectional and Longitudinal Models //////////

//calculate cross covariance for an individual visit
// [[Rcpp::export]]
cube ccv(const mat& xt, int lag){ 
  //xt is matrix for a given visit moving across ROIs and down time
  //lag is maximum number of lags to use in the calculations
  
  int P = xt.n_cols; //number of ROIs
  int T = xt.n_rows; //length of time series
  cube covMat(P, P, 2*lag+1, fill::zeros);
  vec xi(T, fill::zeros);
  double xibar;
  vec xj(T, fill::zeros);
  double xjbar;
  
  for(int i = 0; i < P; i++){
    for(int j = 0; j < P; j++){
      xi = xt.col(i);
      xj = xt.col(j);
      xibar=mean(xi);
      xjbar=mean(xj);
      for(int l = 0; l <= lag; l++){
        vec xil = xi.head(T-l);
        vec xjl = xj.tail(T-l);
        double sum = 0;
        for(int t=0; t<(T-l); t++){
          sum += (xil(t)-xibar)*(xjl(t)-xjbar);
        }
        covMat(i,j,lag+l) = covMat(j,i,lag-l) = sum/T;//cross covariance for ROI i and j at lag l
      }
    }
  }
  return covMat; //return cube with ROIs in rows and columns and lags in slices (lag 0 is middle slice)
}


//function to calculate the lag 0 cross correlation matrix for a given visit
// [[Rcpp::export]]
mat ccr(const cube& C){
  //C is a PxPxL cross covariance cube output by the ccv function
  
  int P = C.n_cols; //Number of ROIs
  int L = C.n_slices; //2 times the number of lags plus 1
  mat R(P,P,fill::zeros);
  mat lag0cov = C.slice((L-1)/2); //take the lag 0 slice from covariance cube
  for(int i = 0; i < P; i++){
    for(int j = 0; j < P; j++){
      R(i,j)=lag0cov(i,j)/std::sqrt(lag0cov(i,i)*lag0cov(j,j)); //convert to correlation
    }
  }
  return R; //return correlation matrix for lag 0
}

 
//calculate weights to be used in the estimation of delta
// [[Rcpp::export]]
arma::vec MB(double lag, double b){
  //lag is maximum number of lags to use in calculation
  //b is the bandwidth of the windowing function times sqrt(T) where T is the length of the time series
  
  vec weights(2*lag+1, fill::zeros);
  for(int i = -lag; i < lag+1; i++){
    if(std::abs(i/b)<=1){weights[i+lag] = 1-std::abs(i/b);}
    else{weights[i+lag] = 0;}
  }
  return weights; //return vector of weights of length 2 times lag + 1 (lag 0 in center)
}

  
//function to calculate delta which are used in calculation of Roy covariance cov(r_ij, r_lm)
// [[Rcpp::export]]
double delta(int i, int j, int l, int m, const vec& weights, const cube& C, int lag, int lagln){
    //i, j, l, and m are the ROIs for the desired delta value
    //weights is a vector of output weights from the MB function
    //C is cross covariance cube from ccv function
    //lag is maximum lag to be used in calculations
    //lagln is 2*lag+1
    
    double sum=0;
    for(int t=0; t<lagln; t++){
      sum += weights(t)*C(i-1,j-1,t)*C(l-1,m-1,t); //estimate theta
    }
    double result = sum/std::sqrt(C(i-1,i-1,lag)*C(j-1,j-1,lag)*C(l-1,l-1,lag)*C(m-1,m-1,lag)); //convert to correlation
    return result;
}


//function calculate Roy covariance cov(r_ij, r_lm)
// [[Rcpp::export]]
double roycov(int i, int j, int l, int m, const cube& C, const mat& R, double T, 
              const vec& weights, int lag, int lagln){
    //i,j,l,m are ROIs of the desired covariance
    //C is cross covariance cube from ccv function
    //R is lag 0 cross correlation from ccr function
    //T is the length of the time series
    //weights is a vector of weights from the MB function(or other windowing function)
    //lag is maximum lag to be used in calculations
    //lagln is 2*lag+1
    
    //formula for Roy covariance
    double result=(0.5*R(i-1,j-1)*R(l-1,m-1)*(delta(i,l,i,l,weights,C,lag,lagln)+
                                              delta(i,m,i,m,weights,C,lag,lagln)+
                                              delta(j,l,j,l,weights,C,lag,lagln)+
                                              delta(j,m,j,m,weights,C,lag,lagln)) -
      R(i-1,j-1)*(delta(i,l,i,m,weights,C,lag,lagln)+delta(j,l,j,m,weights,C,lag,lagln)) -
      R(l-1,m-1)*(delta(j,l,i,l,weights,C,lag,lagln)+delta(j,m,i,m,weights,C,lag,lagln)) +
      delta(i,l,j,m,weights,C,lag,lagln) + delta(j,l,i,m,weights,C,lag,lagln))/T;
    return result;
}


//construct within subject covariance matrix using Roy covariance function, unstructured
// [[Rcpp::export]]
mat RoyMat(const cube& C, const mat& R, double T, const vec& weights, int lag, string SigmaType){
  //C is cross covariance cube from ccv function
  //R is lag 0 cross correlation matrix from ccr function
  //T is the length of the time  series
  //weights is a vector of weights from the MB function(or other windowing function)
  //lag is maximum lag to be used in calculations
  //SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"
  
  int P=C.n_rows; //number of ROIs
  int Q=P*(P-1)/2; //number of ROI pairs
  int lagln=2*lag+1; 
  
  mat result(Q,Q,fill::zeros);
  
  if(SigmaType=="Unstructured"){
    vec resultvec(Q*Q, fill::zeros);
  
    //create vector of roy covariances
    int count=0;
    for(int i=0; i<(P-1); i++){
      for(int j=i+1; j<P; j++){
        for(int l=0; l<(P-1); l++){
          for(int m=l+1; m<P; m++){
            if((l > i) || (l == i && m >= j)){
            resultvec(count) = roycov(i+1,j+1,l+1,m+1,C,R,T,weights,lag,lagln);
            count += 1;
            }
          }
        }
      }
    }
  
    //populate within subject Roy variance matrix from the vector of Roy variances
    for(int i=0; i<Q; i++){
      for(int j=i; j<Q; j++){
        result(i,j)=resultvec[0];
        result(j,i)=resultvec[0];
        resultvec=resultvec.tail(resultvec.n_elem-1);
      }
    }
  }
  
  if(SigmaType=="Diagonal"){
    int count=0;
    for(int i=0; i<(P-1); i++){
      for(int j=i+1; j<P; j++){
        result(count,count) = roycov(i+1,j+1,i+1,j+1,C,R,T,weights,lag,lagln);
        count += 1;
      }
    }
  }
  
  if(SigmaType=="Zero"){
    result.fill(0);
  }
  
  return result;
}


//Apply previous functions to get within subject Roy variance for a single visit
// [[Rcpp::export]]
List Roy(const mat& x, int lag, double bw, string SigmaType){ 
  //x is matrix for a given subject moving across ROIs and down time
  //lag is maximum number of lags to use in the calculations
  //SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"
  //bw is the bandwidth to be used by the MB function
  
  int T=x.n_rows; //length of time series
  lag=min(lag, T);
  cube C=ccv(x,lag); //cross covariance cube
  mat R=ccr(C); //lag 0 correlation matrix
  
  double b = bw*sqrt(T); //bandwidth used in windowing function
  vec weights=pow(MB(lag,b),2); 
  
  List out(2);
  out("Sigma")=RoyMat(C,R,T,weights,lag,SigmaType);
  
  out("R")=R;
  return out;
}


//function to calculate test statistic for vector beta
// [[Rcpp::export]]
double getStat(const vec& beta1, const vec& beta2, const mat& beta1var, const mat& beta2var){
  //beta1 and beta1var are the estimates and variances of beta for group 1
  //beta2 and beta2var are the estimates and variances of beta for group 2
  
  mat Contrast=eye(beta1.n_elem, beta1.n_elem);
  mat VarComb = beta1var + beta2var; 
  vec betadiff = beta1 - beta2;
  mat T(1,1);
  //calculate test statistic
  T = trans(Contrast*betadiff) * inv(Contrast*VarComb*Contrast) * Contrast*betadiff; 
  return T(0,0);
}


//function to calculate test statistic for scaler beta
// [[Rcpp::export]]
double getStatUni(double beta1, double beta2, double beta1var, double beta2var){
  //beta1 and beta1var are the estimates and variances of beta for group 1
  //beta2 and beta2var are the estimates and variances of beta for group 2
  
  return pow(beta1 - beta2, 2) / (beta1var + beta2var);
}


//function to estimate Psi0
// [[Rcpp::export]]
mat get_Psi0(const mat& error, const cube& RoyVars, int Visits, int Q, string Psi0type){
  //error is a residual matrix where each column represents a visit
  //RoyVars is a cube with each slice being the Roy within visit covariance matrix for one visit
  //Visits is the total number of visits (or number of subjects in cross sectional model)
  //Q is the number of ROI pairs  
  //Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  
  mat Empir = error*error.t()/Visits; //empirical error covariance matrix given last beta
  mat AvgRoy = mean(RoyVars,2); //average Roy variance matrix across all N subjects
  
  mat Psi0(Q, Q);
  
  if(Psi0type == "Zero"){//Solve for Psi0 if it is assumed to be 0
    Psi0.fill(0);
  }
  
  if(Psi0type == "Unstructured"){//Solve for Psi0 if unstructured structure is assumed
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        Psi0(i,j) = Empir(i,j) - AvgRoy(i,j);
      }
    }
    for(int i = 0; i < Q; i++){
      if(Psi0(i,i) <= 0){Psi0(i,i) = 0;}
    }
  }
  
  if(Psi0type == "Diagonal"){//Solve for Psi0 if diagonal structure is assumed    
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        if(i == j){Psi0(i,j) = Empir(i,j) - AvgRoy(i,j);}else{Psi0(i,j) = 0;}
      }
    }
    for(int i = 0; i < Q; i++){
      if(Psi0(i,i) <= 0){Psi0(i,i) = 0;}
    }
  }
  
  if(Psi0type == "Scaled"){//Solve for Psi0 if scaled identity structure is assumed 
    vec diffsigma2(Q);
    double sigma2;
    for(int i = 0; i < Q; i++){
      diffsigma2(i) = Empir(i,i) - AvgRoy(i,i);
    }
    if(mean(diffsigma2) < 0){sigma2 = 0;}else{sigma2 = mean(diffsigma2);}
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        if(i == j){Psi0(i,j) = sigma2;}else{Psi0(i,j) = 0;}
      }
    }
  }
  
  if(Psi0type == "CS"){//Solve for Psi0 if compound symmetry structure is assumed 
    double kappa=0;
    if(Q != 1){
      vec diffkappa(Q*(Q-1)/2);
      int count = 0;
      for(int i = 0; i < Q; i++){
        for(int j = i+1; j < Q; j++){
          diffkappa(count) = Empir(i,j) - AvgRoy(i,j); 
          count++;
        }
      }
      kappa = mean(diffkappa);
    }
    vec diffsigma2(Q);
    double sigma2;
    for(int i = 0; i < Q; i++){
      diffsigma2(i) = Empir(i,i) - AvgRoy(i,i);
    }
    if(mean(diffsigma2) < 0){sigma2 = 0;}else{sigma2 = mean(diffsigma2);}
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        if(i == j){Psi0(i,j) = sigma2;}else{Psi0(i,j) = kappa;}
      }
    }
  }
  
  return Psi0;
}



/////////////////////////   Functions used in cross sectional model   //////////////////

//function to estimate beta assuming diagonal design matrix in cross sectional model
// [[Rcpp::export]]
List getBeta(int N, int Q, const cube& invTotVar, const mat& r){
  //N is number of subjects in the group
  //Q is number of ROI pairs
  //invTotVar is a cube where each slice is the inverse of Roy+Psi for one subject
  //r is an NxQ matrix of the observed correlations
  
  mat SumInvVar=sum(invTotVar,2); //sum up inverse covariance matrices
  mat betaVar=inv(SumInvVar); //invert to estimate the variance matrix of beta
  
  //estimate beta
  cube temp(Q,1,N,fill::zeros);
  for(int n = 0; n < N; n++){
    temp.slice(n)=invTotVar.slice(n)*r.row(n).t(); 
  }
  vec tempsum=sum(temp,2);
  vec beta=betaVar*tempsum;
  
  return List::create(
    _("betaVar")=betaVar,
    _("beta")=beta
  );
}


//function to iteratively update beta and Psi using method of moments for cross sectional model
// [[Rcpp::export]]
List getPsiBeta(int N, int Q, const mat& r, const cube& RoyVars, double step2_tol, 
                string Psi0type, int MaxIter){  
  //N is number of subjects in the group
  //Q is number of ROI pairs
  //r is matrix of observed correlations with a row for each subject
  //RoyVars is a cube where each slice is the estimated Roy within subject variance for one subject
  //step2_tol is used for the stopping criteria
  //Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //MaxIter is maximum iterations to perform
  
  bool done = FALSE;
  int iter = 0; // Current iteration
  int code = 0; // 0 - converged, 1 - otherwise
  
  //initiate estimates
  List betaout(2);
  vec beta(Q,fill::zeros);
  vec betaold(Q,fill::zeros);
  mat Psi0(Q,Q,fill::zeros);
  mat Psi0old(Q, Q, fill::zeros);
  cube TotVar(Q,Q,N,fill::zeros);
  cube invTotVar(Q,Q,N,fill::zeros);
  
  betaold=mean(r,0).t();
  mat error(Q, N, fill::zeros);
  for(int n = 0; n < N; n++){
    for(int q = 0; q < Q; q++){
      error(q,n) = r(n,q) - betaold(q);
    }
  }
  
  //iterate until convergence or max number of iterations hit
  while(!done){
    for(int n = 0; n < N; n++){
      for(int i = 0; i < Q; i++){
        error(i,n) = r(n,i) - beta(i);
      }
    }
    
    Psi0=get_Psi0(error, RoyVars, N, Q, Psi0type=Psi0type); //update Psi0
    for(int n = 0; n < N; n++){
      TotVar.slice(n)=RoyVars.slice(n)+Psi0;
      invTotVar.slice(n)=inv(TotVar.slice(n));
    }
    
    betaout=getBeta(N, Q, invTotVar, r); 
    beta=as<vec>(betaout("beta")); //update beta
    
    for(int n = 0; n < N; n++){
      for(int q = 0; q < Q; q++){
        error(q,n) = r(n,q) - beta(q);
      }
    }
    
    //stop if change in estimates is small
    if((norm(beta-betaold,2) < step2_tol) && 
       (abs(Psi0 - Psi0old).max() < step2_tol)){
      done = TRUE;
      code = 0;
    }
    iter++;
    //stop if max number of iterations is reached
    if(iter >= MaxIter){
      done = TRUE;
      code = 1;
    }
    betaold = beta;
    Psi0old = Psi0;
  }
  
  return List::create(
    _["beta"] = beta,
    _["betaVar"] = as<mat>(betaout("betaVar")),
    _["Psi0"] = Psi0,
    _["TotVar"] = TotVar,
    _["code"] = code
  ); 
}


//function to permute data for permutation test for cross sectional model
// [[Rcpp::export]]
List permute(int N, int Q, const mat& r, const cube& RoyVars, vec ids){  
  //N is number of subjects in both groups combined
  //Q is number of ROI pairs
  //r is observed correlations with a row for each subject
  //RoyVars is a cube of within subject Roy variances with a slice for each subject
  //ids is a vector containing (0, 1,... , n-1)
  
  vec ids_new=shuffle(ids); //permute indices
  mat r_new(N,Q);
  cube RoyVars_new(Q,Q,N);
  int id;
  for(int n = 0; n < N; n++){
    id=ids_new(n);
    r_new.row(n)=r.row(id); //reorder r
    RoyVars_new.slice(n)=RoyVars.slice(id); //reorder RoyVars
  } 
  
  return List::create(
    _["RoyVars_new"] = RoyVars_new,
    _["r_new"] = r_new
  ); 
}


//function to permute data and calculate test statistic for permuted data for cross sectional model
// [[Rcpp::export]]
double permT(int N1, int N2, int Q, const mat& r, const cube& RoyVars, double step2_tol, 
             vec ids, string Psi0type, bool verbose, int MaxIter){  
  //N1 is number of subjects in group 1
  //N2 is number of subjects in group 2
  //Q is number of ROI pairs
  //r is observed correlations with a row for each subject
  //RoyVars is a cube of within subject Roy variances with a slice for each subject
  //step2_tol is used for stoping criteria for solving for beta and Psi
  //ids is a vector containing (0, 1,... , n-1)  
  //Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //verbose is a boolean to denote if permutation progress updates should be printed
  //MaxIter is the maximum number of IWLS iterations to run
  
  int N=N1+N2; 
  List permdat=permute(N, Q, r, RoyVars, ids); //permute the data
  mat r_new=permdat("r_new");
  cube RoyVars_new=permdat("RoyVars_new");
  
  //recalculate beta1 and beta1var with new group assignment
  List PsiBeta1_new=getPsiBeta(N1, Q, r_new.head_rows(N1), RoyVars_new.head_slices(N1), 
                               step2_tol, Psi0type, MaxIter);
  vec beta1_new=PsiBeta1_new("beta");
  mat beta1Var_new=PsiBeta1_new("betaVar");
  
  //recalculate beta2 and beta2var with new group assignment
  List PsiBeta2_new=getPsiBeta(N2, Q, r_new.tail_rows(N2), RoyVars_new.tail_slices(N2), 
                               step2_tol, Psi0type, MaxIter);
  vec beta2_new=PsiBeta2_new("beta");
  mat beta2Var_new=PsiBeta2_new("betaVar");
  
  //recalculate test statistic with new group assignments
  double out=getStat(beta1_new, beta2_new, beta1Var_new, beta2Var_new);
    
  return out;
}


//master function for cross sectional model
// [[Rcpp::export]]
List FCanalysis(const List& datalist, int N1, int N2, int Nperms, int lag=50, 
              double bw=5, double step2_tol=1e-5, string Psi0type="CS", int MaxIter=20,
              bool verbose=false, string SigmaType="Unstructured") {
  //datalist is a list where each element is a TxP matrix of data for one subject
  //N1 is the number of subjects in group 1
  //N2 is the number of subjects in group 2
  //Nperms is the number of permutations to run for the permutation test
  //lag is the maximum number of lags to be used in calculations
  //bw is the bandwidth to be used for the MB windowing function
  //step2_tol is the stoping criteria for iteratively solve for beta and Psi
  //Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //MaxIter is maximum iterations to perform solving for Psi and beta
  //verbose is boolean to signify if permutation test status updates should be printed to the console
  //SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"
  
  ///Define dimension parameters and output list
  List outList;
  int N = N1+N2;
  int P=as<mat>(datalist(0)).n_cols;
  int Q=P*(P-1)/2;
  
  ///Estimate Roy Variances
  cube RoyVars(Q,Q,N);
  cube R(P,P,N);
  for(int n = 0; n < N; n++) {
    List Royn=Roy(datalist(n), lag, bw, SigmaType);
    RoyVars.slice(n) = as<mat>(Royn("Sigma"));
    R.slice(n) = as<mat>(Royn("R"));
  }
  cube RoyVarsG1(Q,Q,N1);
  RoyVarsG1=RoyVars.head_slices(N1);
  cube RoyVarsG2(Q,Q,N2);
  RoyVarsG2=RoyVars.tail_slices(N2);
  outList("SigmaG1")=RoyVarsG1;
  outList("SigmaG2")=RoyVarsG2;
  
  ////build r estimate matrix
  mat r(N,Q);
  for(int n = 0; n < N; n++){
    int count=0;
    for(int j = 0; j < P; j++){
      for(int k = j+1; k < P; k++){
        r(n,count)=R(j,k,n);
        count++;
      }
    }
  }
  mat rG1=r.head_rows(N1);
  mat rG2=r.tail_rows(N2);
  outList("rG1")=rG1;
  outList("rG2")=rG2;
  
  ///Estimate Psi and beta for both groups
  List PsiBetaG1=getPsiBeta(N1, Q, rG1, RoyVarsG1, step2_tol, Psi0type, MaxIter);
  vec betaG1=PsiBetaG1("beta");
  mat betaG1Var=PsiBetaG1("betaVar");
  mat Psi0G1=PsiBetaG1("Psi0");
  cube TotVarG1=PsiBetaG1("TotVar");
  double codeG1=PsiBetaG1("code");  
  if((codeG1==1) & (MaxIter>1)){cout << "WARNING: Group 1 estimates did not converge in " << MaxIter << " iterations" << endl;}
  
  List PsiBetaG2=getPsiBeta(N2, Q, rG2, RoyVarsG2, step2_tol, Psi0type, MaxIter);
  vec betaG2=PsiBetaG2("beta");
  mat betaG2Var=PsiBetaG2("betaVar");
  mat Psi0G2=PsiBetaG2("Psi0");
  cube TotVarG2=PsiBetaG2("TotVar");
  double codeG2=PsiBetaG2("code");
  if((codeG2==1) & (MaxIter>1)){cout << "WARNING: Group 2 estimates did not converge in " << MaxIter << " iterations" << endl;}
  
  outList("Psi0G1")=Psi0G1;
  outList("Psi0G2")=Psi0G2;
  outList("TotalVarianceG1")=TotVarG1;
  outList("TotalVarianceG2")=TotVarG2;
  outList("betaG1var")=betaG1Var;
  outList("betaG1")=betaG1;
  outList("betaG2var")=betaG2Var;
  outList("betaG2")=betaG2;
  
  //Calculate Test Statistic
  double T;
  T=getStat(betaG1, betaG2, betaG1Var, betaG2Var);
  outList("T")=T;
  
  //Run permutation test
  vec T_dist(Nperms, fill::zeros);
  vec ids(N,fill::zeros);
  for(int n = 0;n < N; n++){
    ids(n)=n;
  }
  for(int nperm = 0; nperm < Nperms; nperm++){
    T_dist(nperm)=permT(N1, N2, Q, r, RoyVars, step2_tol, ids, Psi0type, verbose, MaxIter);
    if(verbose==true){cout << "Permutation " << nperm+1 << " of " << Nperms << endl;}
  }
  outList("Perm_dist")=T_dist;
  
  return outList;
}





/////////////////////////////////   Functions used in Longitudinal Models   //////////////////////////

//Function to estimate Psi1 for longitudinal model
// [[Rcpp::export]]
mat get_Psi1(const vec& visvec, int Visits, int maxVis, int N, int Q, const mat& error, 
              string Psi1type){
  //visvec is a vector with the number of visits for each subject
  //Visits is the number of total visits, i.e. sum(visvec)
  //maxVis is the maximum number of visits by any one subject
  //N is the number of subjects
  //Q is the number of ROI pairs
  //error is a residual matrix where each column represents a visit
  //Psi1type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  
  //create a cube of Psi1 estimates from each possible combination or visits
  cube Psi1cube(Q, Q, maxVis*(maxVis-1)/2, fill::zeros);
  int countouter = 0;
  int countinner;
  mat v1mat;
  mat v2mat;
  vec v2ind(N, fill::zeros);
  int Nv2Tot = 0;
  int Nv2;
  mat AvgPsi1(Q, Q, fill::zeros);
  
  if(Psi1type != "Zero"){
    for(int v1 = 0; v1 < maxVis-1; v1++){
      for(int v2 = v1+1; v2 < maxVis; v2++){
        for(int n = 0; n < N; n++){
          if(visvec(n) > v2){v2ind(n) = 1;}else{v2ind(n) = 0;}
        }
        Nv2 = sum(v2ind); //number of subjects with a visits at v1 and v2
        Nv2Tot += Nv2;    
        v1mat.resize(Q, Nv2);
        v2mat.resize(Q, Nv2);
        countinner = 0;
        for(int n = 0; n < N; n++){
          if(v2ind(n) == 1){
            v1mat.col(countinner) = error.col(sum(visvec.head(n))+v1);
            v2mat.col(countinner) = error.col(sum(visvec.head(n))+v2);
            countinner++;
          }
        }
        Psi1cube.slice(countouter) = v1mat*v2mat.t();
        countouter++;
      }
    }
    AvgPsi1=sum(Psi1cube,2)/Nv2Tot;
  }
  
  mat Psi1(Q,Q,fill::zeros);  //Initialize Psi1
  
  if(Psi1type == "Zero"){
    Psi1.fill(0);
  }
  
  if(Psi1type == "CS"){ //Solve for Psi1 if compound symmetry structure for Psi1 is assumed
    double omega;               //Initialize omega, the off diagonal element of Psi1
    double nu;                  //Initialize nu, the diagonal element of Psi1
    if(Q == 1){omega = 0;} //if P=2 then Q=1 and Psi1 is a scaler so there are no off diagonal elements
    else{
      vec omegavec(Q*(Q-1)/2);
      int count = 0;
      for(int i = 0; i < Q; i++){
        for(int j = i+1; j < Q; j++){
          omegavec(count) = AvgPsi1(i,j); 
          count++;
        }
      }
      omega = mean(omegavec); //omega is mean of all off diagonal values of all slices of Psi1cube
    }
    vec nuvec(Q);
    for(int i = 0; i < Q; i++){
      nuvec(i) = AvgPsi1(i,i);
    }
    nu = mean(nuvec); //nu is mean of all diagonal values of all slices of Psi1cube
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        if(i == j){Psi1(i,j) = nu;}else{Psi1(i,j) = omega;}
      }
    }
  }
  
  if(Psi1type == "Scaled"){ //Solve for Psi1 if scaled identity structure for Psi1 is assumed
    double nu;                  //Initialize nu, the diagonal element of Psi1
    vec nuvec(Q);
    for(int i = 0; i < Q; i++){
      nuvec(i) = AvgPsi1(i,i);
    }
    nu = mean(nuvec); //nu is mean of all diagonal values of all slices of Psi1cube
    for(int i = 0; i < Q; i++){
      Psi1(i,i) = nu;
    }
  }
  
  if(Psi1type == "Unstructured"){//Solve for Psi1 if unstructured structure is assumed
    Psi1=AvgPsi1;
  }
  
  if(Psi1type == "Diagonal"){//Solve for Psi1 if diagonal structure is assumed
    Psi1=AvgPsi1;
    for(int i = 0; i < Q; i++){
      for(int j = 0; j < Q; j++){
        if(i != j){Psi1(i,j) = 0;}
      }
    }
  }
  
  return Psi1;
}



//function to iteratively update beta and Psi using method of moments for longitudinal model
// [[Rcpp::export]]
List getPsiBetaLongit(int N, int Q, const mat& r, const cube& RoyVars, double step2_tol, 
                      int MaxIter, const vec& visvec, const mat& X,
                      string Psi0type, string Psi1type, int Visits, int maxVis){  
  //N is number of subjects in the group
  //Q is number of ROI pairs
  //r is matrix of observed correlations with a row for each visit
  //RoyVars is a cube where each slice is the estimated Roy within visit variance for one visit
  //step2_tol is used for the stopping criteria
  //MaxIter is maximum iterations to perform
  //visvec is a vector with the number of visits for each subject
  //X is the design matrix
  //Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //Psi1type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //Visits is the total number of visits, i.e. sum visvec
  //maxVis is the maximum number of visits for any one subject
  
  bool done = FALSE;
  int iter = 0; // Current iteration
  int code = 0; // 0 - converged, 1 - otherwise
  
  //Initiate beta estimates at the OLS estimates
  vec Y = vectorise(r.t());
  vec betaold = inv(X.t() * X) * X.t() * Y;  
  mat fits = X * betaold;  
  mat error(Q, Visits, fill::zeros);
  for(int v = 0; v < Visits; v++){
    for(int q = 0; q < Q; q++){
      error(q,v) = r(v,q) - fits(v*Q+q);
    }
  }
  
  //Initiate other parameters
  vec beta(betaold.n_elem, fill::zeros);
  mat betaVar(betaold.n_elem, betaold.n_elem, fill::zeros);
  mat Psi0old(Q, Q, fill::zeros);
  mat Psi1old(Q, Q, fill::zeros);
  mat Psi0(Q, Q);
  mat Psi1(Q, Q);
  mat TotVar(Q*Visits, Q*Visits, fill::zeros);
  mat InvTotVar(Q*Visits, Q*Visits, fill::zeros);
  
  //iterate until convergence or max number of iterations hit
  while(!done){
    Psi0=get_Psi0(error, RoyVars, Visits, Q, Psi0type); //update Psi0
    Psi1=get_Psi1(visvec, Visits, maxVis, N, Q, error, Psi1type); //update Psi1
    int counter=0;
    for(int n = 0; n < N; n++){
      for(int v = 0; v < visvec(n); v++){
        TotVar.submat((counter+v)*Q, (counter+v)*Q, size(Q,Q))=RoyVars.slice(counter+v)+Psi0;
        for(int v2 = v+1; v2 < visvec(n); v2++){
          TotVar.submat((counter+v)*Q, (counter+v2)*Q, size(Q,Q)) = Psi1;
          TotVar.submat((counter+v2)*Q, (counter+v)*Q, size(Q,Q)) = Psi1;
        }
      }
      InvTotVar.submat(counter*Q, counter*Q, size(Q*visvec(n), Q*visvec(n))) = 
        inv(TotVar.submat(counter*Q, counter*Q, size(Q*visvec(n), Q*visvec(n))));
      counter += visvec(n);
    }
    
    betaVar = inv(X.t()*InvTotVar*X); //Estimate variance of beta
    beta = betaVar*X.t()*InvTotVar*Y; //Update beta
    fits = X * beta;  
    for(int v = 0; v < Visits; v++){
      for(int q = 0; q < Q; q++){
        error(q,v) = r(v,q) - fits(v*Q+q);
      }
    }
    
    //stop if change in estimates is small
    if((norm(beta-betaold,2) < step2_tol) && 
       (abs(Psi0 - Psi0old).max() < step2_tol) && 
       (abs(Psi1 - Psi1old).max() < step2_tol)){
      done = TRUE;
      code = 0;
    }
    iter++;
    //stop if max number of iterations is reached
    if(iter >= MaxIter){
      done = TRUE;
      code = 1;
    }
    betaold = beta;
    Psi0old = Psi0;
    Psi1old = Psi1;
  }

  return List::create(
    _["beta"] = beta,
    _["betaVar"] = betaVar,
    _["Psi0"] = Psi0,
    _["Psi1"] = Psi1,
    _["TotVar"] = TotVar,
    _["code"] = code,
    _["Residuals"] = error.t()
  ); 
}


//function to resample the data for bootstrap or permutation test
// [[Rcpp::export]]
List Resample(int N, int Q, const mat& r, const cube& RoyVars, const mat& Xall, 
              const vec& visvec, bool replace, const vec& ids){
  //N is the number of subjects 
  //Q is the number of ROI pairs
  //r is a matrix of observed correlation coefficients with a row for each visit
  //RoyVars is a cube of estimated within visit Roy variances with a slice for each visit
  //Xall is the design matrices from both groups stacked
  //visvec is a vector with the number of visits for each subject
  //replace is a boolean, if false then data is permuted, if true then data is bootstraped
  //ids is a vector with integer values 0 through N-1
  
  vec ids_new(N);
  if(replace == true){
    vec randoms(N);
    randoms.randu();
    ids_new = floor(randoms*N);
  }
  if(replace == false){
    ids_new = shuffle(ids); //permute indices
  }
  
  vec visvec_new(N);
  for(int n = 0; n < N; n++){
    visvec_new(n) = visvec(ids_new(n));
  }
  int SampVisits=sum(visvec_new);
  
  double counter=0;
  double lookup=0;
  mat r_new(SampVisits,Q);
  cube RoyVars_new(Q,Q,SampVisits);
  mat Xall_new(SampVisits*Q, Xall.n_cols);
  
  for(int n = 0; n < N; n++){
    lookup = sum(visvec.head(ids_new(n)+1)) - visvec(ids_new(n));
    for(int v = 0; v < visvec_new(n); v++){
      r_new.row(counter+v) = r.row(lookup+v);                 //reorder r
      RoyVars_new.slice(counter+v) = RoyVars.slice(lookup+v); //reorder RoyVars
      Xall_new.submat((counter+v)*Q, 0, size(Q,Xall.n_cols)) = Xall.submat((lookup+v)*Q, 0, size(Q,Xall.n_cols));
    }
    counter+=visvec_new(n);
  } 
  
  return List::create(
    _["RoyVars_new"] = RoyVars_new,
    _["r_new"] = r_new,
    _["Xall_new"] = Xall_new,
    _["visvec_new"] = visvec_new
  ); 
}


//function to calculate test statistic on permuted full model residuals for multivariate test
// [[Rcpp::export]]
List PermTestLongit(int N1, int N2, int Q, const mat& r_resid_new, const cube& RoyVars_new, 
                   const mat& Xall_new, const vec& visvec_new, double step2_tol, 
                   string Psi0type, string Psi1type, int MaxIter, 
                   const vec& betaG1, const vec& betaG2, int start, int end,
                   int Group1Visits_new, int Group2Visits_new, int Group1maxVis_new,
                   int Group2maxVis_new){  
  //N1 and N2 are the number of subjects in group 1 and group 2 respectively
  //Q is number of ROI pairs
  //r_resid_new are permuted full model residuals with a row for each visit
  //RoyVars_new is permuted a cube of within visit Roy variances with a slice for each visit
  //Xall_new is the permuted design matrices from the two groups stacked, output from Resample
  //visvec_new is a permuted vector with the number of visits for each subject
  //step2_tol is used for stoping criteria for solving for beta and Psi
  //Psi0type and Psi1type denote the structures to be used for Psi0 and Psi1, either "CS" for compound symmetry,
    //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
  //MaxIter is the maximum number of IWLS iterations to run
  //betaG1 and betaG2 are vectors of estimated coefficients from full model for groups 1 and 2 respectively
  //start and end are the first and last indices of beta involved in the test
  //Group1Visits_new and Group2Visits_new are the number of visits in the permuted data for groups 1 and 2 respectively
  //Group1maxVis_new and Group2maxVis_new are the maximum number of visits in the permuted data for groups 1 and 2 respectively
  
  int counter=0;
  vec betaG1null = betaG1;
  betaG1null.subvec(start, end).fill(0);
  vec betaG2null = betaG2;
  betaG2null.subvec(start, end).fill(0);
  vec rG1_null_fit = Xall_new.rows(0, Group1Visits_new*Q-1) * betaG1null;
  vec rG2_null_fit = Xall_new.rows(Group1Visits_new*Q, Xall_new.n_rows-1) * betaG2null;
  vec r_null_fit = join_cols(rG1_null_fit, rG2_null_fit);
  mat r_null(Group1Visits_new + Group2Visits_new, Q, fill::zeros);
  for(int v = 0; v < Group1Visits_new + Group2Visits_new; v++){
    for(int q = 0; q < Q; q++){
      r_null(v, q) = r_resid_new(v, q) + r_null_fit(counter);
      counter ++;
    }
  }
  
  List PsiBetaG1resid=getPsiBetaLongit(N1, Q, r_null.head_rows(Group1Visits_new), 
                                 RoyVars_new.head_slices(Group1Visits_new), step2_tol, 
                                 MaxIter, visvec_new.head(N1), 
                                 Xall_new.head_rows(Group1Visits_new*Q), 
                                 Psi0type, Psi1type, Group1Visits_new, Group1maxVis_new);
  vec betaG1resid=PsiBetaG1resid("beta");
  mat betaG1Varresid=PsiBetaG1resid("betaVar");
  
  List PsiBetaG2resid=getPsiBetaLongit(N2, Q, r_null.tail_rows(Group2Visits_new), 
                                 RoyVars_new.tail_slices(Group2Visits_new), step2_tol, 
                                 MaxIter, visvec_new.tail(N2), 
                                 Xall_new.tail_rows(Group2Visits_new*Q), 
                                 Psi0type, Psi1type, Group2Visits_new, Group2maxVis_new);
  vec betaG2resid=PsiBetaG2resid("beta");
  mat betaG2Varresid=PsiBetaG2resid("betaVar");
  
  double Tglobal = getStat(betaG1resid.subvec(start, end), betaG2resid.subvec(start, end), 
                        betaG1Varresid.submat(start,start,end,end), 
                        betaG2Varresid.submat(start,start,end,end));
  int NPar=end-start+1;
  vec Tlocal(NPar, fill::zeros);
  for(int nPar = 0; nPar < NPar; nPar++){
    Tlocal(nPar) = getStatUni(betaG1resid(start+nPar), betaG2resid(start+nPar), 
           betaG1Varresid(start+nPar, start+nPar), 
           betaG2Varresid(start+nPar, start+nPar));
  }
  
  return List::create(
    _["Tglobal"] = Tglobal,
    _["Tlocal"] = Tlocal
  ); 
}


//function to permute data and calculate global main effect and global interaction test statistics for permuted data for longitudinal model
// [[Rcpp::export]]
List LongitTests(int N1, int N2, int Q, const mat& r_resid, const cube& RoyVars, 
                 const mat& Xall, const vec& visvec, const vec& ids, double step2_tol, 
                 string Psi0type, string Psi1type, int MaxIter, 
                 const vec& betaG1, const vec& betaG2){  
  
  List permdat = Resample(N1+N2, Q, r_resid, RoyVars, Xall, visvec, false, ids);
  
  cube RoyVars_new = permdat("RoyVars_new");
  mat r_resid_new = permdat("r_new");
  mat Xall_new = permdat("Xall_new");
  vec visvec_new = permdat("visvec_new");
  
  int Group1Visits_new = sum(visvec_new.head(N1));
  int Group1maxVis_new = max(visvec_new.head(N1));
  int Group2Visits_new = sum(visvec_new.tail(N2));
  int Group2maxVis_new = max(visvec_new.tail(N2));
  
  List Tmain = PermTestLongit(N1, N2, Q, r_resid_new, RoyVars_new, Xall_new, visvec_new, 
                                  step2_tol, Psi0type, Psi1type, MaxIter, betaG1, 
                                  betaG2, 0, Q-1, Group1Visits_new, Group2Visits_new, 
                                  Group1maxVis_new, Group2maxVis_new);
  double TmainGlobal=Tmain("Tglobal");
  vec TmainLocal=Tmain("Tlocal");
  
  List Tint = PermTestLongit(N1, N2, Q, r_resid_new, RoyVars_new, Xall_new, visvec_new, 
                                  step2_tol, Psi0type, Psi1type, MaxIter, betaG1, 
                                  betaG2, Q, 2*Q-1, Group1Visits_new, Group2Visits_new, 
                                  Group1maxVis_new, Group2maxVis_new);
  double TintGlobal=Tint("Tglobal");
  vec TintLocal=Tint("Tlocal");
  
  return List::create(
    _["TmainGlobal"] = TmainGlobal,
    _["TmainLocal"] = TmainLocal,
    _["TintGlobal"] = TintGlobal,
    _["TintLocal"] = TintLocal
  ); 
}


//master function for longitudinal model
// [[Rcpp::export]]
List FCLongit(const List& datalist, int N1, int N2, int Nperms, const vec& time,
              const vec& visvec, int lag=50, double bw=5, double step2_tol=1e-5, 
              int MaxIter=20, bool verbose=false, string Psi0type="CS", string Psi1type="CS", 
              string SigmaType="Unstructured"){
//datalist is a list where each element is a TxP matrix of data for one subject
//N1 is the number of subjects in group 1
//N2 is the number of subjects in group 2
//Nperms is the number of permutations to run for the permutation test
//time is a vector where each entry is the time (e.g. age or visit number) variable for a visit
//visvec is vector with the number of visits for each subject
//lag is the maximum number of lags to be used in calculations
//bw is the bandwidth to be used for the MB windowing function
//step2_tol is the stoping criteria for iteratively solve for beta and Psi
//MaxIter is maximum iterations to perform solving for Psi and beta
//verbose is boolean to signify if permutation test status updates should be printed to the console
//Psi0type denotes the structure to be used for Psi0, either "CS" for compound symmetry,
  //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
//Psi1type denotes the structure to be used for Psi1, either "CS" for compound symmetry,
  //"Scaled" for Scaled Identity, "Zero", "Unstructured", or "Diagonal"
//SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"


  //Define dimension parameters and output list
  List outList;
  int P=as<mat>(datalist(0)).n_cols;
  int Q=P*(P-1)/2;
  int TotalVisits=sum(visvec);
  int Group1Visits=sum(visvec.head(N1));
  int Group2Visits=sum(visvec.tail(N2));
  int Group1maxVis=max(visvec.head(N1));
  int Group2maxVis=max(visvec.tail(N2));
  
  //Estimate Roy Variances
  cube RoyVars(Q,Q,TotalVisits);
  cube R(P,P,TotalVisits);
  for(int n = 0; n < TotalVisits; n++) {
    List Royn=Roy(datalist(n), lag, bw, SigmaType);
    RoyVars.slice(n) = as<mat>(Royn("Sigma"));
    R.slice(n) = as<mat>(Royn("R"));
  }
  cube RoyVarsG1=RoyVars.head_slices(Group1Visits);
  cube RoyVarsG2=RoyVars.tail_slices(Group2Visits);
  outList("SigmaG1")=RoyVarsG1;
  outList("SigmaG2")=RoyVarsG2;

  //build r matrix of observed correlations
  mat r(TotalVisits,Q);
  for(int n = 0; n < TotalVisits; n++){
    int count=0;
    for(int j = 0; j < P; j++){
      for(int k = j+1; k < P; k++){
        r(n,count)=R(j,k,n);
        count++;
      }
    }
  }
  mat rG1=r.head_rows(Group1Visits);
  mat rG2=r.tail_rows(Group2Visits);
  outList("rG1")=rG1;
  outList("rG2")=rG2;
  
  //build design matrix for each group
  mat X0(Q, Q);
  X0.eye();
  
  mat XG1(Group1Visits*Q , 2*Q, fill::zeros);
  XG1.submat(0,0,size(Group1Visits*Q, Q)) = repmat(X0, Group1Visits, 1);;
  for(int v = 0; v < Group1Visits; v++){
    XG1.submat(v*Q, Q, size(Q, Q)) = time(v) * eye(Q,Q);
  }
  
  mat XG2(Group2Visits*Q , 2*Q, fill::zeros);
  XG2.submat(0,0,size(Group2Visits*Q, Q)) = repmat(X0, Group2Visits, 1);;
  for(int v = 0; v < Group2Visits; v++){
    XG2.submat(v*Q, Q, size(Q, Q)) = time(v) * eye(Q,Q);
  }
  
  mat Xall = join_cols(XG1, XG2);
  
  ///Estimate Psi and beta for both groups
  List PsiBetaG1=getPsiBetaLongit(N1, Q, rG1, RoyVarsG1, step2_tol, MaxIter, visvec.head(N1), 
                                 XG1, Psi0type, Psi1type, Group1Visits, Group1maxVis);
  vec betaG1=PsiBetaG1("beta");
  mat betaG1Var=PsiBetaG1("betaVar");
  mat Psi0G1=PsiBetaG1("Psi0");
  mat Psi1G1=PsiBetaG1("Psi1");
  mat TotVarG1=PsiBetaG1("TotVar");
  double codeG1=PsiBetaG1("code");
  mat ResidsG1=PsiBetaG1("Residuals");
  if((codeG1==1) & (MaxIter>1)){cout << "WARNING: Group 1 estimates did not converge in " << MaxIter << " iterations" << endl;}
  
  List PsiBetaG2=getPsiBetaLongit(N2, Q, rG2, RoyVarsG2, step2_tol, MaxIter, visvec.tail(N2), 
                                 XG2, Psi0type, Psi1type, Group2Visits, Group2maxVis);
  vec betaG2=PsiBetaG2("beta");
  mat betaG2Var=PsiBetaG2("betaVar");
  mat Psi0G2=PsiBetaG2("Psi0");
  mat Psi1G2=PsiBetaG2("Psi1");
  mat TotVarG2=PsiBetaG2("TotVar");
  double codeG2=PsiBetaG2("code");
  mat ResidsG2=PsiBetaG2("Residuals");
  if((codeG2==1) & (MaxIter>1)){cout << "WARNING: Group 2 estimates did not converge in " << MaxIter << " iterations" << endl;}
  
  outList("Psi0G1")=Psi0G1;
  outList("Psi0G2")=Psi0G2;
  outList("Psi1G1")=Psi1G1;
  outList("Psi1G2")=Psi1G2;
  outList("TotalVarianceG1")=TotVarG1;
  outList("TotalVarianceG2")=TotVarG2;
  outList("betaG1var")=betaG1Var;
  outList("betaG1")=betaG1;
  outList("betaG2var")=betaG2Var;
  outList("betaG2")=betaG2;
  
  //Calculate Test Statistics
  double TmainGlobal;
  TmainGlobal = getStat(betaG1.head(Q), betaG2.head(Q), 
                  betaG1Var.submat(0,0,size(Q,Q)), betaG2Var.submat(0,0,size(Q,Q)));
  outList("TmainGlobal")=TmainGlobal;
  
  double TintGlobal = getStat(betaG1.subvec(Q,2*Q-1), betaG2.subvec(Q,2*Q-1), 
                                 betaG1Var.submat(Q,Q,2*Q-1,2*Q-1), 
                                 betaG2Var.submat(Q,Q,2*Q-1,2*Q-1));
  outList("TintGlobal")=TintGlobal;
  
  vec TmainLocal(Q, fill::zeros);
  vec TintLocal(Q, fill::zeros);
  for(int q = 0; q < Q; q++){
    TmainLocal(q) = getStatUni(betaG1(q), betaG2(q), betaG1Var(q, q), betaG2Var(q, q));
    TintLocal(q) = getStatUni(betaG1(Q+q), betaG2(Q+q), betaG1Var(Q+q, Q+q), betaG2Var(Q+q, Q+q));
  }
  outList("TmainLocal")=TmainLocal.t();
  outList("TintLocal")=TintLocal.t();
  
  //Run permutation test
  mat r_resid = join_cols(ResidsG1,ResidsG2);
  double mean_resid = mean(mean(r_resid,0));
  r_resid = r_resid - mean_resid;
  outList("Residuals")=r_resid;
  
  
  vec TmainGlobal_dist(Nperms, fill::zeros);
  mat TmainLocal_dist(Nperms,Q,fill::zeros);
  vec TintGlobal_dist(Nperms, fill::zeros);
  mat TintLocal_dist(Nperms,Q,fill::zeros);
  
  vec ids(N1+N2,fill::zeros);
  for(int n = 0;n < (N1+N2); n++){
    ids(n)=n;
  }
  List perm(4);
  
  for(int nperm = 0; nperm < Nperms; nperm++){
    perm=LongitTests(N1, N2, Q, r_resid, RoyVars, Xall, visvec, ids, step2_tol, Psi0type, 
                     Psi1type, MaxIter, betaG1, betaG2);
    TmainGlobal_dist(nperm)=perm("TmainGlobal");
    TmainLocal_dist.row(nperm)=as<vec>(perm("TmainLocal")).t();
    TintGlobal_dist(nperm)=perm("TintGlobal");
    TintLocal_dist.row(nperm)=as<vec>(perm("TintLocal")).t();
    
    if(verbose==true){cout << "Permutation " << nperm+1 << " of " << Nperms << endl;}
  }
  
  outList("TmainGlobal_dist")=TmainGlobal_dist;
  outList("TmainLocal_dist")=TmainLocal_dist;
  outList("TintGlobal_dist")=TintGlobal_dist;
  outList("TintLocal_dist")=TintLocal_dist;
    
  return outList;
}



