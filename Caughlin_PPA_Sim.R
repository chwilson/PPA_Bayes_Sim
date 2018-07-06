library(mvtnorm)   # to draw multivariate normal outcomes
library(raster)    # to plot stuff
library(rasterVis) # to plot fancy stuff
library(ggplot2)   # more fancy plots

### Functions from Petr Keil 
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

# function that simulates the autocorrelated 2D array with a given side,
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


side= 30;  # sides of the raster
distMat <- dist.matrix(side);
distMat[1:10,1:10];

str(distMat)
# Fix output from arbitrary pixel, the distance from that pixel to any other pixel 
# is given by that row 
out1 <- 10^7;
seed_disp <- Alt2Dt2(distMat[1,],a=2,b=1.5);
sum(seed_disp); # should be < 1 
seed_dispMat <- matrix(seed_disp,nrow = side, ncol = side);
seed_inputMat <- drop(10^7)*seed_dispMat;

levelplot(seed_inputMat) # check 

# probability of seed to seedling
recruit_prob <- 0.0001;
seedlingMat <- matrix(0,nrow = side,ncol = side);
for(i in 1:nrow(seedlingMat)){
  for(j in 1:ncol(seedlingMat)){
    seedlingMat[i,j] <- rbinom(1,size =round(seed_inputMat[i,j]),prob=recruit_prob);
  }
}

levelplot(seedlingMat); # looks good

####### Summing over all pixels

### First choose the generative matrix for seed output 
# Correlated surface 
out_corr <- (80 + 20*cor.surface(30,0,0.1))*10^3
levelplot(out_corr)

# Random surface 
out1 <- runif(side*side,min=0,max=10^7);


# setting up an array of matrices, one for each pixel 
seed_inputArray <- array(0, c(side^2,side,side));
str(seed_inputArray)
for(i in 1:(side^2)){ 
  seed_disp <- Alt2Dt2(distMat[i,],a=2,b=1.5);
  seed_dispMat <- matrix(seed_disp,nrow = side, ncol = side);
 seed_inputArray[i,,] <- drop(as.vector(out_corr)[i])*seed_dispMat;
# loops through and generates 
 }

seed_inputMat <- matrix(0, nrow=side, ncol=side)
for (i in 1:side) {
  for (j in 1:side) {
    seed_inputMat[i,j] <- sum(seed_inputArray[,i,j])
  }
}
str(seed_inputMat)
p1 <- levelplot(seed_inputMat) # OK, very smooth surface 

recruit_prob <- (0.001 + 0.0004*cor.surface(30,0,0.1))*10;
p2 <- levelplot(recruit_prob)
seedlingMat <- matrix(0,nrow = side,ncol = side);
for(i in 1:nrow(seedlingMat)){
  for(j in 1:ncol(seedlingMat)){
    seedlingMat[i,j] <- rbinom(1,size =round(seed_inputMat[i,j]),prob=recruit_prob);
  }
}

p3 <- levelplot(seedlingMat); # OK, this is pretty reasonable for 
# each pixel being 30X30m 
library(gridExtra);
grid.arrange(p1,p2,p3);


############### Combining above with a yearly dimension to generate vectors of recruits for 
########### each pixel 

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



levelplot(matrix(seed_inputMat_yr[,2],side,side))


ca1 <- levelplot(matrix(CA_yr[,1],nrow=side,ncol=side));
ca2 <- levelplot(matrix(CA_yr[,2],nrow=side,ncol=side));
ca3 <- levelplot(matrix(CA_yr[,3],nrow=side,ncol=side));
ca4 <- levelplot(matrix(CA_yr[,4],nrow=side,ncol=side));
ca5 <- levelplot(matrix(CA_yr[,5],nrow=side,ncol=side));
ca6 <- levelplot(matrix(CA_yr[,6],nrow=side,ncol=side));
ca7 <- levelplot(matrix(CA_yr[,7],nrow=side,ncol=side));
ca8 <- levelplot(matrix(CA_yr[,8],nrow=side,ncol=side));
ca9 <- levelplot(matrix(CA_yr[,9],nrow=side,ncol=side));
ca10 <- levelplot(matrix(CA_yr[,10],nrow=side,ncol=side));

grid.arrange(ca1,ca2,ca3,ca10);

t <- seq(1,10,1)
 plot(t,CA_yr[runif(1,1,900),])
# Initial decline from year 1 to year 2, and then subsequent 
# increases. 
# So, need to compute edge sets for adjacency information prior
# to ICAR in Bayes. Also, combine CA data into a single matrix
# 

distMat[1,30:40]

rast <- my.rast(matrix(CA_yr[,3],side,side),30) # see function below
levelplot(rast) 
coordinates(rast)

df <- as.data.frame(coordinates(rast))

adj_matrix <- apply(df, 1, function(pt) 
  (pt["x"] == df$x &  abs(pt["y"] - df$y) == 1) |
  #  (abs(pt["x"] - df$x) == 1 &  pt["y"] == df$y) 4-neighborhood
  (abs(pt["x"] - df$x) <= 1 & abs(pt["y"] - df$y) <= 1)
)
diag(adj_matrix) <- 0

str(adj_matrix)
# using igraph to convert to edge sets 
library(igraph)
g  <- graph.adjacency(adj_matrix)
df <- get.data.frame(g)
str(df) # nice! Has reduced a matric with 407044 cells to a list of
# two vectors with 3990 observations in each 
N_edges <- length(df$from) # edges between adjacent cells
N <- 900 # number of pixels 
str(CA_yr)

simul_stanDat <- list(
N_edges <- length(df$from), # edges between adjacent cells
node1 <- df$from,
node2 <- df$to,
nPIX <- 900,
nYEARS <- 9, 
time <- seq(1,9,1),
CA <- t(CA_yr[,2:10]), 


D <- 1,
c1 <- 1.79,
plot_area <- 900,

recruit_mu <- 100,
recruit_sigma <- 50,
gro_mu <- 0.5,
gro_sigma <- 0.2,
mort_mu <- 0.3,
mort_sigma <- 0.12
)
library(rstan)  
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
getwd()
setwd("C:/Users/Nick/Desktop/Caughlin PPA/Simulation Study")
PPA_reforest_Bayes1 <- stan(file = "PPA_reforest_Bayes.stan",
                            data = simul_stanDat, chains = 4,
                            iter = 500)

print(PPA_reforest_Bayes1, pars = c("muR", "muSUR", "muG", "sigma_obs"))
pairs(PPA_reforest_Bayes1, pars = c("muR", "muSUR", "muG", "sigma_obs"))
save(PPA_reforest_Bayes1, file = "PPA_reforest_Bayes_test1.RData")



















## Summing over time from all pixels
out_mean <- 10^3;
out_sd <- 0.3*out_mean;
output_surface <- cor.surface(30, global.mu = 0, lambda = 0.1);
outSurface <- out_mean + out_sd*output_surface;
levelplot(outSurface);
out1 <- as.vector(outSurface);

library(cowplot)
p1 <- levelplot(seed_inputMat)
p2 <- levelplot(seedlingMat)
str(p1)
plot_grid(outSurface,p1,p2)
#recruit_surface <- cor.surface(30, global.mu = 0.00001, lambda = 0.5)
#levelplot(output_surface);

#str(recruit_surface)
#recruit_prob <- matrix(runif(side^2, min = 0.000005, max = 0.0001),
#                       nrow = side, ncol = side);
recruit_prob <- 0.00001;

# constant or spatially-varying 
levelplot(seedlingMat)
seedlingMatnew <- seedlingMat;
time <- 5; # call this years
for(t in 1:time){
  
  seed_inputArray <- array(0, c(side^2,side,side));
  for(i in 1:(side^2)){
    seed_disp <- Alt2Dt2(distMat[i,],a=10,b=1.5);
    seed_dispMat <- matrix(seed_disp,nrow = side, ncol = side);
    seed_inputArray[i,,] <- drop(out1[i])*seed_dispMat;
    # loops through and generates 
  }
  
  seed_inputMat <- matrix(0, nrow=side, ncol=side)
  for (i in 1:side) {
    for (j in 1:side) {
      seed_inputMat[i,j] <- sum(seed_inputArray[,i,j])
    }
  }
  
  seedlingMatnew <- seedlingMat + seedlingMatnew;
  seedlingMat <- matrix(0,nrow = side,ncol = side);
  for(i in 1:nrow(seedlingMat)){
    for(j in 1:ncol(seedlingMat)){
      seedlingMat[i,j] <- rbinom(1,size =round(seed_inputMat[i,j]),prob=recruit_prob);
    }
  }
}
# Sums the recruits 
levelplot(seedlingMatnew);








# Thinning out to study autocorrelation 
seedlingMatThin <- seedlingMatnew[5:25,5:25];
levelplot(seedlingMatThin)

distMatThin <- distMat[5:25,5:25];

resids <- as.vector((seedlingMatThin) - mean((seedlingMatThin)));
distances <- as.vector(distMatThin);
plot(distances,resids)

ggplot(data = data.frame(x = distances, y=resids),aes(x=x,y=y)) + geom_point() +
  geom_smooth(se=F)

# seedling growth and survival 

# canopy area 


## OK, now I need to have output from all cells dispersing into 





# function that converts a matrix to raster and scales its sides to max.ext
my.rast <- function(mat, max.ext)
{
  rast <- raster(mat)
  rast@extent@xmax <- max.ext
  rast@extent@ymax <- max.ext
  return(rast)
}


# Dispersal distance kernel
Alt2Dt <- function(x,a,b,s=1){
  disp <- s*2*pi*x*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} # sums to 1

# Dispersal density kernel
Alt2Dt2 <- function(x,a,b,Q=1){
  disp <- Q*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} # does not sum to 1 



Alt2Dt2_NDD <- function(x,a,b,k){
  disp <- (1/(1+((1/0.01)-1)*exp(-k*x)))*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} #
integrate(Alt2Dt2_NDD, 0, 20, a = 2, b=2, k= 2)



Alt2Dt_NDD <- function(x,a,b,k){
  disp <- (2*pi*x)*(1/(1+((1/0.01)-1)*exp(-k*x)))*((b-1)/(pi*a^2))*(1 +(x^2/a^2))^-b;
  return(disp);
} 

a <- 2
b <- 2
curve((1/(1+((1/0.01)-1)*exp(-4*x))),0,5)
curve(Alt2Dt(x,a=2,b=2),0,5)

## Panel 1 Density and Survivorship 
Density_Survivor <- ggplot(data.frame(x=c(0,10)), aes(x)) + 
  stat_function(fun = Alt2Dt2, args = list(a=2,b=1.1,Q = 10^2)) +
  stat_function(fun = function(x){(1/(1+((1/0.01)-1)*exp(-3*x)))}) +
  ylab("Density (Seeds/Area or Survival)") + xlab("radial distance")


## Panel 2 comparing JC effect from just NDD to the area correction
Area_JC <- ggplot(data.frame(x = c(0,10)),aes(x)) + 
  stat_function(fun = Alt2Dt2_NDD, args = list(a=2,b=1.1,k=3)) +
  stat_function(fun = Alt2Dt, args = list(a=2,b=1.1,s=0.1),color = "blue") +
  ylab("Location Expectation (blue) versus NDD J-C effect") +
  xlab("radial distance")

library(cowplot)
plot_grid(Density_Survivor, Area_JC)

### These kernels are related by 2*pi*x in 2D
library(ggplot2)
#ggplot(data.frame(x = c(0,100)),aes(x)) + stat_function(fun = Clark2Dt, args = list(u = 4))
ggplot(data.frame(x = c(0,10)),aes(x)) + stat_function(fun = Alt2Dt2, args = list(a = 1, b = 2.5)) +
  stat_function(fun = Alt2Dt, args = list(a = 1, b = 2.5), color = "red") 
# stat_function(fun = Alt2Dt2_NDD, args = list(a=2,b=1.5,k=0.5), color = "blue")



integrate(Alt2Dt,0,20,a=2,b=2)


