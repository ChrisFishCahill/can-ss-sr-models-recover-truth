/*
   State-space Ricker stock-recruit model
      based on Su and Peterson (2012)
            Cahill QFC MSU
 */
data {
  int<lower=0> k;
  int<lower=0> n_year;
  vector[n_year] E;
  vector[n_year-k] C; 
  // vector[k] ar_prior; 
  // vector[k] ln_sdp_prior; 
  // vector[k] ln_sdo_prior; 
}
parameters {
  real ar; 
  vector[k] ln_So; 
  real ln_sdp;
  real ln_sdo; 
  vector<lower=C>[n_year-k] R; 
}
transformed parameters {
 vector[n_year - k] mu; 
 vector[n_year] S;
 real b = ar/ exp(ln_So[1]); 
 real sdo = exp(ln_sdo); 
 real sdp = exp(ln_sdp); 
 real hmsy = 0.5*ar - 0.07*pow(ar, 2); 
 real smsy = hmsy / b; 
 
 // initialize 
 for(t in 1:n_year){S[t] = exp(ln_So[k]);}
 
 for(t in 1:(n_year-k)){
   mu[t] = ar + log(S[t]) - b*S[t];   
   S[t+k] = R[t] - C[t]; //R and C are t + k NOT t
 }
}
model {
  // priors all vague by default
  // ar ~ normal(ar_prior[1], ar_prior[2]); 
  // ln_sdp ~ normal(ln_sdp_prior[1], ln_sdp_prior[2]);
  // ln_sdo ~ normal(ln_sdo_prior[1], ln_sdo_prior[2]);
  // ln_So ~ normal(ln_So_prior[1], ln_So_prior[2]);   
  
  // likelihoods
  R ~ lognormal(mu, exp(ln_sdp));
  E ~ lognormal(log(S), exp(ln_sdo));
}


