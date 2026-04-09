//time-varying model

data {
  int nA; // N age groups
  int nT; // N time points
  int nL; // N locations
  int nP; // N province
  int N;
  array[nL] int prov_N; //province index vector
  // int cases[nL, nT, nA]; // reported case data
  array[nL, nT, nA] int cases_inc; //incidence
  // array [nL, nT, nA] int pop; // population data
  array[2, nA] int ageLims; // lower & upper bounds of age groups
  row_vector[100] age; // age as sequence from 0 to 99
  array [nL, nT]int sum_cases;
  array[nA] real age_band;
  int X; //reference age group with no age modifier in reportint
}

parameters {
  array[nL] real<lower=0.001,upper=0.25> lam_H; // historic average FOI for each location
  array[nL, nT] real<lower=0.001,upper=0.25> lam_t; // time varying FOI for each location and time
  array[nP] real<lower=0, upper=1> rho; // reporting rate of 2nd infections (same for all locations)
  array[nP] real<lower=0, upper=1> gamma; // relative reporting rate of 1st infections (same for all locations)
  array[nL] real<lower=0.45,upper=2> chi1;
  array[nL] real<lower=0.45,upper=2> chi2;
}

transformed parameters {
 array[nL] matrix<lower=0>[nT, nA] Ecases; // expected reported cases per year and age group
 array[nL] vector<lower=0>[nT] Ecases_sum; // expected reported cases per year
 array[nL]  matrix<lower=0>[nT, nA] prob_Ecases; // expected probability of reported cases per year and age group

{
  // Define location-specific arrays
 array[nL] row_vector[100] susc0; // proportion susceptible at time 0 for each location
 array[nL] row_vector[100] mono0; // proportion monotypic at time 0 for each location
 array[nL] row_vector[100] multi0; // proportion multitypic at time 0 for each location
 array[nL] matrix[nT, 100] susc; // proportion susceptible
 array[nL] matrix[nT, 100] mono; // proportion monotypic
 array[nL] matrix[nT, 100] multi; // proportion multitypic
 array[nL] matrix[nT, 100] inc1; // incidence of primary infections
 array[nL] matrix[nT, 100] inc2; // incidence of secondary infections  
 

  for (l in 1:nL) {
    // Compute initial proportions for each location
    susc0[l] = exp(-4 * lam_H[l] * age);
    mono0[l] = 4 * exp(-3 * lam_H[l] * age) .* (1 - exp(-lam_H[l] * age));
    multi0[l] = 1 - (susc0[l] + mono0[l]);
    
    // Compute proportions over time for each location
    for (t in 1:nT) {
      susc[l, t, 1] = 1;
      mono[l, t, 1] = 0;
      multi[l, t, 1] = 0;
    }

        susc[l,1,2:100] =  susc0[l,1:99] - 4 * lam_t[l,1] * susc0[l,1:99];
        mono[l,1,2:100] = mono0[l,1:99] + 4 * lam_t[l,1] * susc0[l,1:99] - 3 * lam_t[l,1] * mono0[l,1:99];
        multi[l,1,2:100] = multi0[l,1:99] + 3 * lam_t[l,1] * mono0[l,1:99];

        inc1[l,1,  : ] = 4 * lam_t[l,1] * susc0[l];
        inc2[l,1,  : ] = 3 * lam_t[l,1] * mono0[l];


      for (t in 2:nT) { for (age_index in 2:100) { 
        susc[l, t, age_index] = susc[l, t-1, age_index - 1] - 4 * lam_t[l, t] * susc[l, t-1, age_index - 1];
        mono[l, t, age_index] = mono[l, t-1, age_index - 1] + 4 * lam_t[l, t] * susc[l, t-1, age_index - 1] - 3 * lam_t[l, t] * mono[l, t-1, age_index - 1];
        multi[l, t, age_index] = multi[l, t-1, age_index - 1] + 3 * lam_t[l, t] * mono[l, t-1, age_index - 1];
      }

      inc1[l, t, :] = 4 * lam_t[l, t] * susc[l, t-1, :];
      inc2[l, t, :] = 3 * lam_t[l, t] * mono[l, t-1, :];
    } 



    // Expected reported cases
    for (t in 1:nT) {
      for (a in 1:(X-1)) {
      Ecases[l, t, a] = N * rho[prov_N[l]] * chi1[l] * (mean(inc2[l, t, ageLims[1, a]:ageLims[2, a]]) + gamma[prov_N[l]] * mean(inc1[l, t, ageLims[1, a]:ageLims[2, a]])); 
      }
      
    Ecases[l, t, X] = N * rho[prov_N[l]] *(mean(inc2[l, t, ageLims[1, X]:ageLims[2, X]]) + gamma[prov_N[l]] * mean(inc1[l, t, ageLims[1, X]:ageLims[2, X]])); 
    
    for (a in (X+1):nA) {
      Ecases[l, t, a] = N * rho[prov_N[l]] * chi2[l] * (mean(inc2[l, t, ageLims[1, a]:ageLims[2, a]]) + gamma[prov_N[l]] * mean(inc1[l, t, ageLims[1, a]:ageLims[2, a]])); 
      }
      
    }

  
    for (t in 1 : nT)  {
      for (a in 1 : nA) {

      if (Ecases[l,t, a] == 0) {
        Ecases[l,t, a] = 0.0001;

      } else {
        Ecases[l,t, a] = Ecases[l,t, a];
      }}

  Ecases_sum[l, t] = sum(Ecases[l, t, :]);
  prob_Ecases[l, t, :] = Ecases[l, t, :] / Ecases_sum[l, t];
  
    }}
  
}

}


model {
  
  // Priors
  for (L in 1:nL) lam_H[L] ~ normal(0, 0.1);
  for (L in 1:nL) for (T in 1:nT)  lam_t[L,T] ~ normal(0, 0.1);
  for (L in 1:nP) rho[L] ~ normal(0.8, 0.2);
  for (L in 1:nP) gamma[L] ~ normal(0.5, 0.2);
   for (L in 1:nL) chi1[L] ~ normal(1, 1);
   for (L in 1:nL) chi2[L] ~ normal(1, 1);

  // Likelihood
  for (l in 1 : nL) {
    for (t in 1 : nT) {
       target += poisson_lpmf(sum_cases[l,t] | Ecases_sum[l,t]);
       target += multinomial_lpmf(cases_inc[l, t, :] | to_vector(prob_Ecases[l, t, :]));
    }
  }
}



generated quantities{
  
  array[nL, nT] real log_lik;
  
   for (l in 1 : nL) {
    for (t in 1 : nT) {
       log_lik[l,t] = poisson_lpmf(sum_cases[l,t] | Ecases_sum[l,t]) + multinomial_lpmf(cases_inc[l, t, :] | to_vector(prob_Ecases[l, t, :]));
    }
  }
  
  
}



