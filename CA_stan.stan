data{
  int<lower = 0> N; 
  real PSN[N];
  real CA[N]; 
}

parameters{
  real alpha;
  real Beta; 
  real Betas;
  real<lower = 1> phi; 
}

model{
  
  for(i in 1:N){
      CA[i] ~ beta(phi*inv_logit(alpha + Beta*PSN[i]), phi*(1-inv_logit(alpha+Beta*PSN[i]))); 
    }
  }
