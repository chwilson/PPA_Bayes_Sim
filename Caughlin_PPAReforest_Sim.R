
library(mvtnorm)   # to draw multivariate normal outcomes
library(raster)    # to plot stuff
library(rasterVis) # to plot fancy stuff


########### Simulation of reforestation per PPA 
## In this version, recruitment follows a spatially autocorrelated surface
## as does initial canopy area. Seed output is function of canopy area,
## and seed input is the result of summing from all sources using dispersal 
## density kernel. 
## Growth and mortality are held constant and there is no process error 


nYears <- 10; 
seedlingMat_yr <- array(0, c(side^2,nYears));
seed_inputMat_yr <- array(0, c(side^2,nYears));
CA_yr <- array(0, c(30^2,nYears));
recruit_prob <- (0.007 + 0.002*cor.surface(30,0,0.1));

levelplot(recruit_prob) # Fixed recruitment probability on correlated surface 
 seedlingMat_yr[,1] <- seedlingMat/100; ## Initializes simulation  with recruits
CA_yr[,1] <- 0.1 + 0.05*cor.surface(30,0,0.05); ## Initializes simulation with CA 

M <- 0.05 
G <- 0.5


for(k in 2:nYears){ 
  
  # Compute output as a function of canopy area in time t - 1 
  out_corr <- matrix(CA_yr[,k-1]*10^5,nrow = side, ncol = side); 
  
  # setting up an array of matrices, one for each pixel for seed inputs
  seed_inputArray <- array(0, c(side^2,side,side));
  for(i in 1:(side^2)){ 
    seed_disp <- Alt2Dt2(distMat[i,],a=1.2,b=2.5);
    seed_dispMat <- matrix(seed_disp,nrow = side, ncol = side);
    seed_inputArray[i,,] <- drop(as.vector(out_corr)[i])*seed_dispMat;
    # loops through and generates 
  }
  
  seed_inputMat <- matrix(0, nrow=side, ncol=side);
  for (i in 1:side) {
    for (j in 1:side) {
      seed_inputMat[i,j] <- sum(seed_inputArray[,i,j]);
    }
  }
  #str(seed_inputMat);
  
  seedlingMat <- matrix(0,nrow = side,ncol = side);
  for(i in 1:nrow(seedlingMat)){
    for(j in 1:ncol(seedlingMat)){
      seedlingMat[i,j] <- rbinom(1,size =round(seed_inputMat[i,j]),prob=recruit_prob);
    }
  }
  
  seedlingMat_yr[,k] <- seedlingMat;
  seed_inputMat_yr[,k] <- seed_inputMat; 
  
  for(pp in 1:(side^2)){
    CA_yr[pp,k] <- sum_cohorts(seedlingMat_yr[pp,1:k]/10,k)$ca2[k];
    # This is cohort summing function developed earlier that represents the
    # discrete form of PPA 
    
  }
  
}

############### Functions used above 


# Dispersal density kernel
Alt2Dt2 <- function(x,a,b,Q=1){
  disp <- Q*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} # does not sum to 1 


# Cohort summing function for discrete PPA 

ag <- 0
M <- 0.2
D <- 1
c1 <- 1.79
G <- 0.3
plot_area <- 900

R <- c(5,8,100,150,20,10,10,8,7,15)
R <- muRs
time <- 12


sum_cohorts <- function(R,time){
  nmu <- matrix(0,length(R),time) # each row is a cohort path 
  
  ca <- matrix(0,length(R),time) # each row is a cohort path 
  # added at a different time
  ca2 <- rep(0, time) # each row is a cohort path 
  
  t_mat <- matrix(0,length(R),time) # checking the diagonalized
  # time matrix
  
  x <- seq(1,time,1)
  
  for(i in 1:length(R)){
    for(j in 1:length(x)){
      if(j >= i) {
        nmu[i,j] <-  exp(-M*(x[j-i+1]-ag))*c1*(D + G*(x[j-i+1]-ag))/plot_area
        ca[i,j] <-  R[i]*exp(-M*(x[j-i+1]-ag))*c1*(D + G*(x[j-i+1]-ag))/plot_area
        t_mat[i,j] <-  x[j-i+1]
        
      } else{
        nmu[i,j] <- 0
        ca[i,j] <- 0
      }
    }
    ca2[i] <-  sum(ca[,i])
  }
  return(list(ca=ca,t=t_mat, nmu = nmu, ca2 = ca2))
}


#function that simulates the autocorrelated 2D array with a given side,
# and with exponential decay given by lambda
# (the mean mu is constant over the array, it equals to global.mu)
cor.surface <- function(side, global.mu, lambda)
{
  D <- dist.matrix(side)
  # scaling the distance matrix by the exponential decay
  SIGMA <- exp(-lambda*D)
  mu <- rep(global.mu, times=side*side)
  # sampling from the multivariate normal distribution
  M <- matrix(nrow=side, ncol=side)
  M[] <- rmvnorm(1, mu, SIGMA)
  return(M)
}
