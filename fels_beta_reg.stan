functions {
  real Btilde(int M, int j, real t){
    real a = (0.5/(M+1))*exp(beta_lpdf((t+1)/2 | j+1, M-j+1));
    return(a);
  }
  
  vector Balpha(int M, real t) {
  vector[M+1] bb = rep_vector(0.1, M+1); 
  for(j in 1:(M+1)){
    bb[j] = Btilde(M,j-1,t);
   } 
  return(bb);
  }
}

data{
  int<lower=0> m;                         // No of subjects                   
  int<lower=0> p;                         // No of predictors (p)
  int<lower=0> N;                         // Total no of observations including all the time points
  int T_vec[m];                           // irregular no of time points
  vector<lower=0, upper=1>[N] y;          // response variable (including missing values)
  int<lower=0, upper=1> y_obs[N];         // missingness indicator for Y
  matrix[N,p] x;                          // predictor matrix
  int<lower=0> M;                         // no of basis
  matrix[M+1,M+1] A; 
}

transformed data {
  int n_obs = sum(y_obs);  // number of observed cases
  int ns[n_obs];           // indices of observed cases
  int ny = 1; 
  for (n in 1:N) {
    if (y_obs[n]) {
      ns[ny] = n;
      ny += 1;
    }
  }
}

parameters{
 vector<lower=0>[M+1] phi;                // basis coefficients 
 vector[p] beta;                          // regression coefficients
 real<lower=0> psi;                       // parameter of beta distribution other than mean  
 vector[m] z;                             //random effects 
 real<lower=0> sigmaz;                    // sd of random effects  
}

transformed parameters{
  vector<lower=0, upper=1>[N] mu;     // Mean of beta regression
  real<lower=0> c[N];                 // 1st parameter of beta density
  real<lower=0> d[N];                 // 2nd parameter of beta density
  vector[p] alpha;                    // normalized regression coefficients    
  vector<lower=0>[N] psit;   
  vector[N] zt; 
  
  real betanorm = sqrt(sum(square(beta)));
  alpha = beta/betanorm;               
  
  for(i in 1:m){
    zt[(sum(T_vec[1:(i-1)])+1):(sum(T_vec[1:i]))] = rep_vector(z[i],T_vec[i]);
  }
  
  psit = rep_vector(psi,N);
  
  for (i in 1:N) {
  real si = (Balpha(M,x[i,]*alpha))' * (A*phi);  // Single index
  mu[i] = inv_logit(si + zt[i]);                 // model of mean of the beta distn
  }

  for (i in 1:N) {
    c[i] = mu[i] * psit[i];
    d[i] = (1-mu[i]) * psit[i];
  }

}

model {
  
 y[ns] ~ beta(c[ns], d[ns]);
 
 phi[1] ~ normal(0,5);
 for(i in 2:(M+1)){
 phi[i] ~ normal(0,5)T[0,];  // basis coefficients 
 }
 beta ~ multi_normal(rep_vector(0,p), diag_matrix(rep_vector(0.1,p)) ); // regression coefficients
 psi ~ gamma(2.5,0.5); 
 
 // construct random effects
 z ~ multi_normal(rep_vector(0,m), diag_matrix(rep_vector(sigmaz,m)));
}

