data {
  int N;
  real theta_hat[N];
  real sigma[N];
}


parameters {
  real theta[N];
  real mu;
  real<lower = 0> nu;
}


model {
  for(i in 1:N){
    theta_hat[i] ~ normal(theta[i], sigma[i]);
    theta[i] ~ normal(mu, nu);
  }
  
  mu ~ cauchy(0, 10);
  nu ~ cauchy(0, 10);
}

