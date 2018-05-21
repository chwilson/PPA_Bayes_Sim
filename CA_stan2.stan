data{
  int<lower = 0> N; 
  real PSN[N];
  real CA[N]; 
  real mm[N]; 
}

parameters{
  real alpha;
  real alpha2; 
  real Beta; 
  real beta2; 
  real<lower = 1> phi; 
  real beta3;
  real Betas; 
//  real<lower = 0, upper = 1> theta; 
  real alpha3;
  real beta4; 
  real beta5; 
}

model{
  
  
  for(i in 1:N){
   if (CA[i] < 0.05)
   1 ~ bernoulli(inv_logit(alpha2 + beta2*PSN[i] + beta3*mm[i])); 
  # else { if(CA[i] > 0.95)
  #        1 ~ bernoulli(inv_logit(alpha3 + beta4*PSN[i] + beta5*mm[i]));
   else{0 ~ bernoulli(inv_logit(alpha2 + beta2*PSN[i] + beta3*mm[i]));
             # 0 ~ bernoulli(inv_logit(alpha3 + beta4*PSN[i] + beta5*mm[i]));
    CA[i] ~ beta(phi*Phi(alpha + Beta*PSN[i] + Betas*mm[i]), phi*(1-Phi(alpha+Beta*PSN[i] + Betas*mm[i]))); 
      }
   }
}
//}