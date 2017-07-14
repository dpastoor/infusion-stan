data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> nSubjects;
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
  real WT[nt];  // FIX ME - should one WT per subject, not event
}

transformed data{
  int nIIV = 2;
  vector[nObs] logCObs = log(cObs);
  int nTheta = 3;
  int nCmt = 2;
  int nti[nSubjects];
  real biovar[nCmt];
  real tlag[nCmt];
  
  real ka = 1.0;

  for(i in 1:nSubjects) nti[i] = last[i] - first[i] + 1;

  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0, upper = 500> CLHat;
  real<lower = 0, upper = 3500> VHat;
  real<lower = 0> sigma;
  real<lower = 0> sigma_prop;

  # Inter-Individual variability
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0.01, upper = 2>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters{
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM;  // Matt's trick
  real<lower = 0> theta[nTheta];
  matrix<lower = 0>[nt, nCmt] x;
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  vector<lower = 0>[nObs] sigma_prop_cHatObs_square;

  thetaHat[1] = CLHat;
  thetaHat[2] = VHat;

  ## Matt's trick to use unit scale 
  thetaM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))'; 

  for(j in 1:nSubjects)
  {
    theta[1] = thetaM[j, 1]; # CL
    theta[2] = thetaM[j, 2]; # V
    theta[3] = ka;

    x[first[j]:last[j]] = PKModelOneCpt(time[first[j]:last[j]], 
                                        amt[first[j]:last[j]],
                                        rate[first[j]:last[j]],
                                        ii[first[j]:last[j]],
                                        evid[first[j]:last[j]],
                                        cmt[first[j]:last[j]],
                                        addl[first[j]:last[j]],
                                        ss[first[j]:last[j]],
                                        theta, biovar, tlag);

    cHat[first[j]:last[j]] = col(x[first[j]:last[j]], 2) ./ theta[2]; ## divide by V
  }

  cHatObs  = cHat[iObs];

  // Construct the variance proportional to ChatObs, and square it.
  // Need to use a loop because there are no ^2 vectorized operation
  // (probably because its an ambiguity)
  for (i in 1:nObs)
    sigma_prop_cHatObs_square[i] = (cHatObs[i] * sigma_prop)^2;
}

model{
  ## Prior
  CLHat ~ lognormal(log(10), 0.25);
  VHat ~ lognormal(log(35), 0.25);

  L ~ lkj_corr_cholesky(1);

  ## Inter-individual variability (see transformed parameters block
  ## for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);

  sigma ~ cauchy(0, 5);
  // logCObs ~ normal(log(cHatObs), sigma);
  sigma_prop ~ cauchy(0, 2);
  logCObs ~ normal(log(cHatObs), sqrt(sigma_prop_cHatObs_square + sigma^2));
}