data{

int<lower = 0> nYEARS; // nyears for L5 matrix
int<lower = 0> nPIX; // npixels for L5 matrix 
int<lower = 0> time[nYEARS]; // counter for yrs
matrix[nYEARS,nPIX] CA;//  canopy area data

// ICAR spatial adjacency information 
  int<lower=0> N_edges;
  int<lower=1, upper=nPIX> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nPIX> node2[N_edges];  

// passed-in priors
real recruit_mu;
real recruit_sigma; 
real gro_mu; 
real gro_sigma;  
real mort_mu;
real mort_sigma; 

real D; 
real c1; 
real plot_area; 
}


parameters{

real B; 
real<lower = 0, upper = 10> sigma_obs; 
real<lower = 0,upper = 0.2> muSUR; //it really wants mortality to be high! how do we avoid that? is the boundary ok?
real<lower = 0> muG; 

// spatial parameters 
 real<lower = 0> muR[nYEARS];
 matrix<lower = 0, upper = 450>[nYEARS,nPIX] R; 

 real<lower=0> tau_theta[nYEARS];   // precision of heterogeneous effects
 real<lower=0> tau_phi[nYEARS];     // precision of spatial effects
 //vector[nPIX - 1] phi_raw[nYEARS]; // raw spatial effects
 matrix[nYEARS,nPIX -1] phi_raw; 

}

transformed parameters {
real<lower = 0> sigma_theta[nYEARS];
real<lower = 0> sigma_phi[nYEARS];
//vector[nPIX] phi[nYEARS];
matrix[nYEARS,nPIX] phi; 

for(i in 1:nYEARS){
  sigma_theta[i] = inv(sqrt(tau_theta[i]));
  sigma_phi[i] = inv(sqrt(tau_phi[i]));
}

  for(i in 1:nYEARS){
  phi[i,1:(nPIX - 1)] = phi_raw[i];
  phi[i,nPIX] = -sum(phi_raw[i]);
  }
  
}


model {
// uses these as local variables but does not store 

matrix[nYEARS,nPIX] mu; 
matrix[nYEARS,nYEARS] nmu[nPIX]; // array of nXn matrices 


for(j in 1:nPIX) {
for(i in 1:nYEARS) {
 for(kk in 1:nYEARS) {

if(kk >= i){
nmu[j,i,kk] = R[i,j]*exp(-muSUR*((kk-i+1)))*c1*(D+muG*((kk-i+1)));
} 
else {
nmu[j,i,kk] = 0;
}
      }

	  }
	}

for(i in 1:nYEARS){
  for(j in 1:nPIX) {
  mu[i,j] = sum(col(nmu[j],i)/plot_area);// + b_year[i]; 
    }
  }



// Priors on growth and mortality 
muSUR ~ normal(mort_mu, mort_sigma);
muG ~ normal(gro_mu, gro_sigma); 


// global mean of recruitment by year 
muR ~ student_t(3,recruit_mu, recruit_sigma);

// spatial model for recruitment 
for(i in 1:nYEARS) {
 R[i] ~ normal(muR[i] + phi[i] * sigma_phi[i], sigma_theta[i]);
  // both a spatial and non-spatial random effect 
  target += -0.5 * dot_self(phi[i,node1] - phi[i,node2]);
}
 
  // priors on precisions for structured and unstructured REs
  // these will need some adjusting here 

  tau_theta ~ gamma(3.28, 1800);  // Carlin WinBUGS priors
  tau_phi ~ gamma(3.28*0.49*6.25, 1800);  

sigma_obs ~ normal(0,1); 
 
for(i in 1:nYEARS) {
for(j in 1:nPIX) {
CA[i,j] ~ normal(mu[i,j], sigma_obs) T[0,];

		}
	}
}