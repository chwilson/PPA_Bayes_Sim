
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
side <- 30; # number of pixels to a side 

seedlingMat_yr <- array(0, c(side,side,nYears));
seed_inputMat_yr <- array(0, c(side,side,nYears));
CA_yr <- array(0, c(side^2,nYears));

recruit_prob <- (0.007 + 0.002*cor.surface(30,0,0.1));
levelplot(recruit_prob) # Fixed recruitment probability on correlated surface 


CA_yr[,1] <- 0.2 + 0.05*cor.surface(30,0,0.05); ## Initializes simulation with CA 
seedlingMat <- matrix(0,side,side);
seedlingMat_yr[,,1] <- seedlingMat/100; ## Initializes simulation  with recruits

distMat <- dist.matrix(side);


M <- 0.05 
G <- 0.5

########## Running the simulations 

for(k in 2:nYears){ 

  # Compute output as a function of canopy area in time t - 1 
  out_corr <- matrix(CA_yr[,k-1]*10^5,nrow = side, ncol = side); 
  
  # setting up an array of matrices, one for each pixel for seed inputs

  seed_inputArray <- array(0, c(side^2,side,side));
  for(i in 1:(side^2)){ 
    # Stochastic dispersal process
    seed_disp <- Alt2Dt2(distMat[i,],a=1.2,b=2.5);
    seed_dispMat <- matrix(seed_disp,nrow = side, ncol = side);
    
    seed_inputArray[i,,] <- drop(as.vector(out_corr)[i])*seed_dispMat;
    # loops through and generates 
    
  }
  
  # convert to a matrix 
  seed_inputMat <- matrix(0, nrow=side, ncol=side);
  for (i in 1:side) {
    for (j in 1:side) {
      seed_inputMat[i,j] <- rpois(1,lambda=sum(seed_inputArray[,i,j]));
    }
  }
  #print(seed_inputMat)
  #str(seed_inputMat);
  
  seedlingMat <- matrix(0,nrow = side,ncol = side);
  for(i in 1:nrow(seedlingMat)){
    for(j in 1:ncol(seedlingMat)){
      # stochastic recruitment process
      seedlingMat[i,j] <- rbinom(1,size =round(seed_inputMat[i,j]),prob=recruit_prob[i,j]);
    }
  }
  
  seedlingMat_yr[,,k] <- seedlingMat;
  seed_inputMat_yr[,,k] <- seed_inputMat; 
  
  seedlingMat_yr2 <- array(seedlingMat_yr, c(side^2,nYears));
  
  for(pp in 1:(side^2)){
    
    # deterministic growth of cohorts according to PPA 
    CA_yr[pp,k] <- sum_cohorts(seedlingMat_yr2[pp,1:k]/10,k)$ca2[k];
    # This is cohort summing function developed earlier that represents the
    # discrete form of PPA 
    
  }
  
}

########### Plotting Output 

# Canopy cover 

CA.list <- list(my.rast(matrix(CA_yr[,1],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,2],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,3],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,4],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,5],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,6],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,7],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,8],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,9],side,side),max.ext=side),
               my.rast(matrix(CA_yr[,10],side,side),max.ext=side))
CC <- stack(CA.list)
names(CC) <- c("year 1","year 2", "year 3",
               "year 4","year 5", "year 6",
               "year 7","year 8", "year 9",
               "year 10")

levelplot(CC)



# Seed Input  

SI.list <- list(my.rast(seed_inputMat_yr[,,1],max.ext=side),
                my.rast(seed_inputMat_yr[,,2],max.ext=side),
                my.rast(seed_inputMat_yr[,,3],max.ext=side),
                my.rast(seed_inputMat_yr[,,4],max.ext=side),
my.rast(seed_inputMat_yr[,,5],max.ext=side),
my.rast(seed_inputMat_yr[,,6],max.ext=side),
my.rast(seed_inputMat_yr[,,7],max.ext=side),
my.rast(seed_inputMat_yr[,,8],max.ext=side),
my.rast(seed_inputMat_yr[,,9],max.ext=side),
my.rast(seed_inputMat_yr[,,10],max.ext=side))
      
SI <- stack(SI.list)
names(SI) <- c("year 1","year 2", "year 3",
               "year 4","year 5", "year 6",
               "year 7","year 8", "year 9",
               "year 10")

levelplot(SI)







############### Functions used above 


# Dispersal density kernel
Alt2Dt2 <- function(x,a,b,Q=1){
  disp <- Q*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} # does not sum to 1 

## A general observation: the more "long-tailed" the dispersal kernel,
## the less spatial structure for *seed input* is evident in the random landscape. This 
## kinda makes sense in that it spreads out the contributions of the 
## higher CA patches, effectively "mixing" or "shmooshing" seed rain 
## together. This implies that spatial structure in canopy structure 
## under reforestation arises either because a) seed dispersal is very 
## short distance, or b) recruitment/growth are very spatially structured. 

## To me, this suggests that the critical piece of field data to collect is,
## ironically, recruitment post seed addition! Dispersal, growth and mortality 
## can be inferred given reasonable ability to monitor forest structure via RS. 





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

### THe functions below are from Petr Keil http://www.petrkeil.com/?p=1910

# function that makes distance matrix for a side*side 2D array  
dist.matrix <- function(side)
{
  row.coords <- rep(1:side, times=side)
  col.coords <- rep(1:side, each=side)
  row.col <- data.frame(row.coords, col.coords)
  D <- dist(row.col, method="euclidean", diag=TRUE, upper=TRUE)
  D <- as.matrix(D)
  return(D)
}

#function that simulates the autocorrelated 2D array with a given side,
# and with exponential decay given by lambda
# (the mean mu is constant over the array, it equals to global.mu

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

# function that converts a matrix to raster and scales its sides to max.ext
my.rast <- function(mat, max.ext)
{
  rast <- raster(mat)
  rast@extent@xmax <- max.ext
  rast@extent@ymax <- max.ext
  return(rast)
}


PV <- rep(0,15*12)
PV[1] <- 27
r <- 0.1
for(i in 2:length(PV)){
  PV[i] <- PV[i-1] + (1/(1+(r/12))^i)*27
}
tail(PV)
plot(PV)


