data{
  int<lower = 0> N; 
  real pvar[N];
  real mm[N]; 
}

parameters{
  real<lower = 0> alpha;
  real<upper = 0> Beta; 
  real<lower = 1> phi; 
}

model{
  for(i in 1:N){
  pvar[i] ~ beta(phi*inv_logit(alpha + Beta*mm[i]), phi*(1-inv_logit(alpha+Beta*mm[i]))); 
  }
}