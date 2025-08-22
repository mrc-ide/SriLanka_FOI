// IMAI MODEL  attempt with upper bound chi from 2 to 1.5, change lower bound to 0.5 and the priors n rhos and gamma

data {
  
  int nA; // N age groups
  int nL; // N locations
  int nP; //number of provinces
  int N;
  int nT;
   array[nL] int prov_N; //province index vector
  // int pop[nL, nT, nA]; // population data
  array[2, nA] int ageLims; // lower & upper bounds of age groups
  // row_vector[100] age; // age as sequence from 0 to 99
  array[nL, nT, nA] int cases_inc; //incidence
  array[nA] real age_band;
}

transformed data {
  
  array[nL] int sum_cases;
  array[nL, nA] int cases; //incidence
  
  for (l in 1:nL) for (a in 1:nA) cases[l,a] = to_int(sum(cases_inc[l,,a])/ nT);
  for (l in 1:nL) sum_cases[l] = sum(cases[l,]);
  
}

parameters {
  
  array[nL] real<lower=0,upper=0.25> lam_H; // historic average FOI for each location
  array[nP] real<lower=0, upper=1> rho; // reporting rate of 2nd infections (same for all locations)
  array[nP] real<lower=0, upper=1> gamma; // relative reporting rate of 1st infections (same for all locations)
  array[nL] real<lower=0.2,upper=1.5> chi1;
  array[nL] real<lower=0.2,upper=1.5> chi2;
  
}



transformed parameters {
  
  array[nL] vector<lower=0>[nA] Ecases; // Expected reported cases per age group per location
  array[nL] real<lower=0> Ecases_tot; // Expected total reported cases per location
  array[nL] vector<lower=0>[nA] Ecases_prob; // Expected probability of reported cases per age group and location

  array[nL] vector<lower=0, upper=1>[nA] inc1; // incidence of primary infections
  array[nL] vector<lower=0, upper=1>[nA] inc2; // incidence of secondary infections
  
  // need to for loop as no way to do elementwise exponential in rstan?
    
   for(l in 1:nL){
   // Compute expected reported cases
    for (a in 1:nA) {
      inc1[l,a] = exp(-4*lam_H[l]*ageLims[1, a]) - exp(-4*lam_H[l]*ageLims[2, a]);
      inc2[l,a] = 4*(exp(-3*lam_H[l]*ageLims[1, a]) - exp(-3*lam_H[l]*ageLims[2, a])) - 3*(exp(-4*lam_H[l]*ageLims[1, a]) - exp(-4*lam_H[l]*ageLims[2, a]));
    }
   
   
    // Expected reported cases -> should we didve by age band???
      for (a in 1:4) {
      Ecases[l, a] = N * rho[prov_N[l]] * chi1[l] * (inc2[l, a] + gamma[prov_N[l]] * inc1[l, a]) * (1/  age_band[a]);
      }
      
     Ecases[l, 5] = N * rho[prov_N[l]]* (inc2[l, 5] + gamma[prov_N[l]] * inc1[l, 5]) * (1/  age_band[5]);
    
    for (a in 6:nA) {
       Ecases[l, a] = N * rho[prov_N[l]] * chi2[l] * (inc2[l, a] + gamma[prov_N[l]] * inc1[l, a]) * (1/  age_band[a]);
      }
      
    
  
  Ecases_tot[l] = sum(Ecases[l,]);

 for(a in 1:nA){
  
  if (Ecases[l,a] == 0){
    Ecases[l,a] = 0.0001;
  } else{
    Ecases[l,a] = Ecases[l,a];
  }
  
  Ecases_prob[l,a] = Ecases[l,a] / Ecases_tot[l];
  
   }}
   }
  
  



model {
  
  //--- priors
  
   for (L in 1:nL) lam_H[L] ~ normal(0, 0.5);
  for (L in 1:nP) rho[L] ~ normal(0.5, 0.5);
  for (L in 1:nP) gamma[L] ~ normal(0.5, 0.5);
   for (L in 1:nL) chi1[L] ~ normal(1, 1);
   for (L in 1:nL) chi2[L] ~ normal(1, 1);


  //--- likelihood on combined cases across all years!
    for (l in 1:nL){
  target += poisson_lpmf(sum_cases[l] | Ecases_tot[l]); //for total number of cases
  target += multinomial_lpmf(cases[l,] | to_vector(Ecases_prob[l,])); //for cases by age-groups
    }
}



