
load(file = "dat.Rdata")
str(dat)

plot(dat$CA/676, dat$PsnCA)
max(dat$CA)
CA_std <- (dat$CA)/max(dat$CA);

m1 <- lm(dat$PsnCA ~ CA_std)
summary(m1)

length(which(dat$CA<0.05))

LiDAR_PSN <- data.frame(CA = dat$CA, PSN = dat$PsnCA)
library(ggplot2)
ggplot(LiDAR_PSN, aes(x=CA,y=PSN)) + geom_point() +
  geom_abline(aes(intercept = coef(m1)[1], slope = coef(m1)[2]), color = "red")



# Coefficients: 0.507[0.008] intercept, 0.301[0.016] slope, 0.12 SE; 
vcov(m1) # since covariance is minimal, can generate predictions with independent draws

CA_levels <- seq(0,1,0.05);
pred_PSN <- matrix(0,nrow = 10^3, ncol = length(CA_levels));
for(i in 1:ncol(pred_PSN)){
pred_PSN[,i] <- coef(m1)[1] + CA_levels[i]*coef(m1)[2] + rnorm(10^3,0,0.1227)
}
apply(pred_PSN,2,quantile,probs = c(0.1,0.5,0.9))


pheno_dat <- read.csv("SQUAREdatT98.csv")

ggplot(pheno_dat, aes(x=mm, y = pvar, group = sqid, color = yearmat)) + geom_point() +
  stat_function(fun = logit_link, args = list(a=1.4,b=-1.9))


pheno_m1 <- lm(pvar ~ mm, pheno_dat)
summary(pheno_m1)


# Examining possibility for logistic modeling of phenology 
logit_link <- function(x,a,b){
  logit <- 1/(1 + exp(-(a + b*x)));
  return(logit); 
}
# Reasonable values for a,b 
ggplot(data = data.frame(x = c(0,1.5)),aes(x)) + stat_function(fun = logit_link, args = list(a = 1.4, b = -1.89))

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Scrubbing NAs 
pheno_dat_noNA <- na.omit(pheno_dat)
# Bringing pvar values off the boundaries for computational reasons
pheno_dat_noNA$pvar[which(pheno_dat_noNA$pvar == 1)] <- pheno_dat_noNA$pvar[which(pheno_dat_noNA$pvar == 1)] - 0.01
pheno_dat_noNA$pvar[which(pheno_dat_noNA$pvar == 0)] <- pheno_dat_noNA$pvar[which(pheno_dat_noNA$pvar == 0)] + 0.01

# Reinstalling Rstan
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
# Checking toolchain 
fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) , '
	return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;
' )
fx( 2L, 5 ) #

getwd()
# Running model in Stan
pheno_datStan <- list(mm = pheno_dat_noNA$mm, pvar = pheno_dat_noNA$pvar, N = length(pheno_dat_noNA$pvar));
pheno_mod <- stan(file = "pheno_stan.stan", data = pheno_datStan, chains = 4, iter = 500);
print(pheno_mod)
pairs(pheno_mod)

library(installr)
updateR()

ggplot(pheno_dat, aes(x=mm, y = pvar, group = sqid, color = yearmat)) + geom_point() +
  stat_function(fun = logit_link, args = list(a=1.4,b=-1.9))



# Structural model in Stan

LiDAR_PSN$CA2 <- LiDAR_PSN$CA/676;
# Brining values off boundaries 
LiDAR_PSN$CA2[which(LiDAR_PSN$CA2 == 1)] <- LiDAR_PSN$CA2[which(LiDAR_PSN$CA2 == 1)] - 0.005
LiDAR_PSN$CA2[which(LiDAR_PSN$CA2 == 0)] <- LiDAR_PSN$CA2[which(LiDAR_PSN$CA2 == 0)] + 0.005


LiDAR_PSN$PSN[which(LiDAR_PSN$PSN == 0)] <- LiDAR_PSN$PSN[which(LiDAR_PSN$PSN == 0)] + 0.005
LiDAR_PSN$PSN[which(LiDAR_PSN$PSN == 1)] <- LiDAR_PSN$PSN[which(LiDAR_PSN$PSN == 1)] - 0.005


# Data
str(LiDAR_PSN)
LiDAR_PSN_stan <- list(CA = LiDAR_PSN$CA2, PSN = LiDAR_PSN$PSN, mm =  N = length(LiDAR_PSN$PSN));

# Fitting model

CA_model_stan <- stan(file = "CA_stan.stan", data = LiDAR_PSN_stan, chains = 4, iter = 500); 
print(CA_model_stan)



CA_model_stan2 <- stan(file = "CA_stan2.stan", data = LiDAR_PSN_stan, chains = 4, iter = 500); 
print(CA_model_stan2)
pairs(CA_model_stan2)

# Plotting model fit 
ggplot(LiDAR_PSN, aes(x=PSN, y=CA2)) + geom_point() +
  stat_function(fun = pnorm_link, args = list(a = -2.06, b = 3.38))

pnorm_link <- function(x,a,b){
  pnormey <- pnorm(a+b*x,0,1);
  return(pnormey);
}

logit_link(0.7,3.24,-7.11)


hurdle_function <- function(x,a,b,a2,b2){
  hurdle_expect <- (1-logit_link(x,a2,b2))*pnorm(a+b*x,0,1);
  return(hurdle_expect);
}

ggplot(LiDAR_PSN, aes(x=PSN, y=CA2)) + geom_point() +
  stat_function(fun = hurdle_function, args = list(a = -2.06, b = 3.38, a2 = 3.24, b2 = -7.11))


hist(dat$PSN)

## 
str(dat)
PSN_1 <- dat$PSN[,1];
plot(PSN_1)

CA <- matrix(0,nrow(dat$PSN),ncol(dat$PSN)); 

for(i in 1:nrow(CA)){
  for(j in 1:ncol(CA)){
   CA[i,j] <- pnorm_link(dat$PSN[i,j] + (logit_link(0,1.4,-1.89)-logit_link((dat$mm[i,j]-0.4),a=1.4,b=-1.89)),-2.07,3.39)
  }
}

str(CA)
library(reshape2)
CA_df <- as.data.frame(CA);
CA_df_melt <- melt(CA_df); 
CA_df_melt$yr <- rep(seq(1,12,1),414);

library(dplyr)
CA_df_melt2 <- CA_df_melt %>% filter(yr>1);

ggplot(CA_df_melt2,aes(x=yr,y=value,group=variable)) + geom_line()
str(CA_df_melt)



# pheno_correct
PSN_pheno_corr <- logit_link(0,1.4,-1.89) - logit_link((dat$mm[,1]-0.4),a=1.4,b=-1.89) 
cbind(PSN_pheno_corr, dat$mm[,1])

PSN_1_pheno <- PSN_1 + PSN_pheno_corr
# CA 

CA <- logit_link(PSN_1_pheno,-4.09,5.88)
plot(CA)





### 

getwd()
setwd("C:/Users/Chris Wilson/Downloads")
load(file="dmean2.76_mm.rdata")
str(dmean)

library(ggplot2)
ggplot(dmean, aes(x=Psn,y=CA/676)) + geom_point()

# Randomly subsample 10K from this
dmean_sub <- dmean[sample(nrow(dmean),10^4,replace=F),];
str(dmean_sub)

full_stanDat <- list(N = length(dmean_sub$Psn),PSN=dmean_sub$Psn,CA=dmean_sub$CA/676,mm=dmean_sub$mm);
range(full_stanDat$CA)

full_stanDat$CA[(which(full_stanDat$CA==1))] <- full_stanDat$CA[(which(full_stanDat$CA==1))] - 0.005;
# Now randomly subsample 10K from this list

full_CAmod <- stan(file="CA_stan2.stan", data = full_stanDat, chains = 4, iter = 500);
str(full_stanDat)



print(full_CAmod)
pairs(full_CAmod)


hurdle_function <- function(x,a,b,a2,b2){
  hurdle_expect <- (1-logit_link(x,a2,b2))*pnorm(a+b*x,0,1);
  return(hurdle_expect);
}

# This one includes the hurdle at 0, and the absorbing boundary at 1 :) 
hurdle_function2 <- function(x,a,b,a2,b2,a3,b3){
  hurdle_expect <- (1-logit_link(x,a2,b2))*(1-logit_link(x,a3,b3))*pnorm(a+b*x,0,1) + logit_link(x,a3,b3);
  return(hurdle_expect);
}

ggplot(dmean_sub, aes(x=Psn, y=CA/676)) + geom_point() +
  stat_function(fun = hurdle_function, args = list(a = -2.75, b = 3.57, a2 = 2.92, b2 = -5.85),
                color = "red",size=2)

ggplot(dmean_sub, aes(x=Psn, y=CA/676)) + geom_point() +
  stat_function(fun = pnorm_link, args = list(a = -2.75, b = 3.57),
                color = "red",size=2)

library(installr)
updateR()
