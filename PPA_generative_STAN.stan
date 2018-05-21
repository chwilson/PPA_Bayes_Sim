data{
  int<lower = 1> yrs;// # years 
  real<lower = 0, upper = 1> PSN[yrs]; 
  real c1;
  real plot_area; 
  real D; 
}

parameters{

real<lower = 0> muR; 
real<lower = 0> muGRO;
real<lower = 0, upper = 1> SUR; 
real<lower = 1, upper = 50> phi; 
  
}

transformed parameters{
  
 real muSUR = 1 - SUR; 
 
  
}

model{

vector[yrs] mu; // define locally and do not save 

muR ~ normal(1,0.5); 
muGRO ~ normal(0.5,0.25); 
SUR ~ beta(8,2); 
phi ~ normal(30,10); 
  
for(i in 1:yrs){
mu[i] = ((c1*(muR/plot_area))*(muGRO+D*muSUR-exp(-(i*muSUR))*(muGRO+D*muSUR + muGRO*(i*muSUR)))/(muSUR^2));
PSN[i] ~ beta(mu[i]*phi,(1-mu[i])*phi);
}
  
}

