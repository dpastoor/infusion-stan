data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> nSubjects;
  // int nIIV;
  int<lower = 1> iObs[nObs];
  int<lower = 1> first[nSubjects];
  int<lower = 1> last[nSubjects];
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  int nIIV = 5;
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;
  int nCmt = 3;
  int nti[nSubjects];
  real biovar[nCmt];
  real tlag[nCmt];

  for(i in 1:nSubjects) nti[i] = last[i] - first[i] + 1;

  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0, upper = 500> CLHat;
  real<lower = 0, upper = 500> QHat;
  real<lower = 0, upper = 3500> V1Hat;
  real<lower = 0, upper = 3500> V2Hat;
  real<lower = 0, upper = 100> kaHat;
  real<lower = 0> sigma;
  
  # Inter-Individual variability
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0.01, upper = 2>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters{
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; // variable required for Matt's trick
  real<lower = 0> theta[nTheta];
  matrix<lower = 0>[nt, nCmt] x;
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  ## Matt's trick to use unit scale 
  thetaM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))'; 
  
  for(j in 1:nSubjects)
  {
    theta[1] = thetaM[j, 1]; # CL
    theta[2] = thetaM[j, 2]; # Q
    theta[3] = thetaM[j, 3]; # V1
    theta[4] = thetaM[j, 4]; # V2
    theta[5] = thetaM[j, 5]; # ka

    x[first[j]:last[j]] = PKModelTwoCpt(time[first[j]:last[j]], 
                                       amt[first[j]:last[j]],
                                       rate[first[j]:last[j]],
                                       ii[first[j]:last[j]],
                                       evid[first[j]:last[j]],
                                       cmt[first[j]:last[j]],
                                       addl[first[j]:last[j]],
                                       ss[first[j]:last[j]],
                                       theta, biovar, tlag);
                                       
    cHat[first[j]:last[j]] = col(x[first[j]:last[j]], 2) ./ theta[3]; ## divide by V1
  }

  cHatObs  = cHat[iObs];
}

model{
  ## Prior
  CLHat ~ lognormal(log(10), 0.25);
  QHat ~ lognormal(log(15), 0.5);
  V1Hat ~ lognormal(log(35), 0.25);
  V2Hat ~ lognormal(log(105), 0.5);
  kaHat ~ lognormal(log(2.5), 1);
  
  L ~ lkj_corr_cholesky(1);
  
  ## Inter-individual variability (see transformed parameters block
  ## for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);
  
  sigma ~ cauchy(0, 5);
  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
    real cObsCond[nObs];
    vector[nt] cHatPred;
    real cObsPred[nObs];
    matrix[nt, 3] xPred;
    matrix[nIIV, nSubjects] etaStdPred;
    matrix<lower=0>[nSubjects, nIIV] thetaPredM;
    corr_matrix[nIIV] rho;
    real<lower = 0> thetaPred[nTheta];

    rho = L * L';

    for(i in 1:nSubjects){
      for(j in 1:nIIV){ 
        etaStdPred[j, i] = normal_rng(0, 1);
      }
    }

    thetaPredM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

    for(j in 1:nSubjects){
      
      thetaPred[1] = thetaPredM[j,1]; # CL
      thetaPred[2] = thetaPredM[j,2]; # Q 
      thetaPred[3] = thetaPredM[j,3]; # V1
      thetaPred[4] = thetaPredM[j,4]; # V2
      thetaPred[5] = thetaPredM[j,5]; # ka 
    
      xPred[first[j]:last[j],] = PKModelTwoCpt(time[first[j]:last[j]],
                                              amt[first[j]:last[j]],
                                              rate[first[j]:last[j]],
                                              ii[first[j]:last[j]],
                                              evid[first[j]:last[j]],
                                              cmt[first[j]:last[j]],
                                              addl[first[j]:last[j]],
                                              ss[first[j]:last[j]],
                                              theta, biovar, tlag);

     cHatPred = xPred[ ,2] / thetaPred[3];
  }
  
  for(i in 1:nObs){
      cObsCond[i] = exp(normal_rng(log(cHatObs[i]), sigma)); # individual predictions
      cObsPred[i] = exp(normal_rng(log(cHatPred[iObs[i]]), sigma)); # population predictions
    }
}