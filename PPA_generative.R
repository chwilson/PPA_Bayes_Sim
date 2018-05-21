### 

############# TO DO on 05/21/18: 1) Get Stan running again!; 2) Re-run CA and pheno models to get 
### reasonable values for phi; 3) Re-run prior predictive simulations using the phi values; 4) fit
### full model in Stan and sample posterior predictive for all parameters; 5) construct histograms of
### rank statistics for ALL parameters; 5) compile into a single PDF; 6) Merge everything into the Rmd 
### document. 




library(reshape2);library(ggplot2); library(dplyr)


PPA_traj <- function(muR,muSUR,muGRO,plot_area = 30,c1 = 1, D = 1.7, years = 15){
  mu <- rep(0,years); 
  for(i in 1:years){
  mu[i] = ((c1*(muR/plot_area))*(muGRO+D*muSUR-exp(-(i*muSUR))*(muGRO+D*muSUR + muGRO*(i*muSUR)))/(muSUR^2));
  }
  return(mu);
}

logit_link <- function(x,a,b){
  logit <- 1/(1 + exp(-(a + b*x)));
  return(logit); 
}

muR <- 2;
SUR <- 0.8;
muSUR <- 1 - SUR;
muGRO <- 0.5; 

CA_true <- PPA_traj(muR,muSUR,muGRO);


# Replicated datasets 
?rbeta
muR <- rnorm(10^3,1,0.5);
SUR <- rbeta(10^3, 8, 2);
muSUR <- 1-SUR;
muGRO_1 <- rnorm(10^3,0.5,0.25);
muGRO_2 <- runif(10^3,0.1,1.1);
#muGRO <- 0.5*muGRO_1 + 0.5*muGRO_2; # convolved prior
muGRO <- muGRO_1; 
phi_CA <- rnorm(10^3, )
phi_pheno <- rnorm(10^3, )
phi_annual <- rnorm(10^3, 30, 10); 

param_mat <- cbind(muR,muSUR,muGRO,phi);


param_mat2 <- ifelse(param_mat >= 0, param_mat, NA);
row.has.na <- apply(param_mat2, 1, function(x){any(is.na(x))});
naRows <- sum(row.has.na);
paramMat.filtered <- param_mat2[!row.has.na,]; 


CA_trueMAT <- matrix(0,10^3,15); 
for(i in 1:(10^3-naRows)){
  CA_trueMAT[i,] <- PPA_traj(paramMat.filtered[i,1],paramMat.filtered[i,2],paramMat.filtered[i,3]);
}


CA_trueDF <- as.data.frame(t(CA_trueMAT));
CA_trueDFmelt <- melt(CA_trueDF);
CA_trueDFmelt$t <- rep(seq(1,15,1),1000);

# Showing all trajectories
ggplot(CA_trueDFmelt, aes(x = t, y = value, group = variable)) + geom_line()

# examining quantiles 

CA_quantiles <- CA_trueDFmelt %>% group_by(t) %>% 
  do(data.frame(t(quantile(.$value, probs = c(0.2,0.5,0.8,0.99)))))

str(CA_quantiles)
ggplot(CA_quantiles, aes(x = t)) + geom_line(aes(y = X50.),color = "red") +
  geom_segment(aes(xend = t, y = X20., yend = X80.),linetype = "dashed")


# OK good 
# Now adding measurement error of three types: 1) Interannual stochasticity;
# 2) linear phenology, and 3) PSN --> CA relationship 

# Start with just interannual stochasticity 

# First need to transform into a constrained 0-1 scale, so CA > 1 --> PSN = max 


# Reasonable values for a,b 
ggplot(data = data.frame(x = c(0,1)),aes(x)) + stat_function(fun = logit_link, args = list(a = 1.4, b = -1.9))



# Starting with a "true" CA matrix, then adding CA --> PSN, then PSN phenology, then interannual stochasticity

CA_trueMAT <- matrix(0,10^3,15); 
CA_obsMAT <- matrix(0,10^3,15); 
for(i in 1:(10^3-naRows)){
  CA_trueMAT[i,] <- ifelse(PPA_traj(paramMat.filtered[i,1],paramMat.filtered[i,2],paramMat.filtered[i,3])<1,
                           PPA_traj(paramMat.filtered[i,1],paramMat.filtered[i,2],paramMat.filtered[i,3]),1);
  for(j in 1:15){
   PSN_trueMAT[i,] <- logit_link(CA_trueMAT[i,], a = -5, b = 7);
   PSN_phenoMAT[i,] <- rbeta(PSN_trueMAT[i,] - logit_link(mm[i,], a = 1.4, b = -1.9); 
 # CA_obsMAT[i,j] <- logit_link(x=rnorm(1,mean=CA_trueMAT[i,j], sd = CA_trueMAT[i,j]*0.2),a=-3,b=8); 
   CA_obsMAT[i,j] <- rbeta(1, shape1 = paramMat.filtered[i,4]*CA_trueMAT[i,j], shape2 = paramMat.filtered[i,4]*(1-CA_trueMAT[i,j])); 
  }
}

CA_obsDF <- as.data.frame(t(CA_obsMAT));
CA_obsDFmelt <- melt(CA_obsDF);
CA_obsDFmelt$t <- rep(seq(1,15,1),1000);
# Showing all trajectories
ggplot(CA_obsDFmelt, aes(x = t, y = value, group = variable)) + geom_line()

# Subset trajectories 
CA_obsDFmelt2 <- CA_obsDFmelt;
str(CA_obsDFmelt2)
vars <- round(runif(20,1,1000));
vars_select <- paste("V",vars, sep = "");
reduce_DF <- CA_obsDFmelt2 %>% filter(variable %in% vars_select);

ggplot(reduce_DF, aes(x = t, y = value, group = variable, color = variable)) + geom_line() +
  xlab("Year") + ylab("Observed PSN") + theme_bw();









# Stan data

PPAGen_data <- list(c1 = 1, D = 1.7, plot_area = 30, PSN = CA_obsMAT[1,], yrs = as.integer(15));
library(rstan);
options(mc.cores = parallel::detectCores());
rstan_options(auto_write = TRUE);
PPAGen_stan <- stan(file = "PPA_generative_STAN.stan", data = PPAGen_data, chains = 4, iter = 500);


rank_stat_muR <- rep(0,nrow(paramMat.filtered));
rank_stat_muSUR <- rep(0,nrow(paramMat.filtered));
rank_stat_muGRO <- rep(0,nrow(paramMat.filtered));
rank_stat_phi <- rep(0,nrow(paramMat.filtered));


# Iterating over the draws from the joint prior 

for(i in 1:length(rank_stat_muR)){
PPAGen_data <- list(c1 = 1, D = 1.7, plot_area = 30, PSN = CA_obsMAT[i,], yrs = as.integer(15));
PPAGen_stan2 <- stan(fit = PPAGen_stan, data = PPAGen_data, chains = 4, iter = 1000, control =
                       list(adapt_delta = 0.95));
#print(PPAGen_stan2)
post_muR <- extract(PPAGen_stan2, pars = "muR")[[1]];
rank_stat_muR[i] <- length(which(post_muR < paramMat.filtered[i,1]));

post_muSUR <- extract(PPAGen_stan2, pars = "muSUR")[[1]];
rank_stat_muSUR[i] <- length(which(post_muSUR < paramMat.filtered[i,2]));

post_muGRO <- extract(PPAGen_stan2, pars = "muGRO")[[1]];
rank_stat_muGRO[i] <- length(which(post_muGRO < paramMat.filtered[i,3]));


post_phi <- extract(PPAGen_stan2, pars = "phi")[[1]];
rank_stat_phi[i] <- length(which(post_phi < paramMat.filtered[i,4])); 

print(i); 
print(c(paste("here are the ranks"),rank_stat_muR[i], rank_stat_muSUR[i],
      rank_stat_muGRO[i], rank_stat_phi[i])); 
}


rank_statistics <- cbind(rank_stat_muR, rank_stat_muSUR, rank_stat_muGRO, rank_stat_phi);
save(rank_statistics, file = "rankStats1.Rdata")

rank_df <- as.data.frame(rank_statistics);

hist1 <- ggplot(rank_df, aes(x=rank_stat_muR)) + geom_histogram(binwidth = 50) + scale_x_continuous(limits = c(1,2000)); 
hist2 <- ggplot(rank_df, aes(x=rank_stat_muGRO)) + geom_histogram(binwidth = 50) + scale_x_continuous(limits = c(1,2000)); 
hist3 <- ggplot(rank_df, aes(x=rank_stat_muSUR)) + geom_histogram(binwidth = 50) + scale_x_continuous(limits = c(1,2000)); 
hist4 <- ggplot(rank_df, aes(x=rank_stat_phi)) + geom_histogram(binwidth = 50) + scale_x_continuous(limits = c(1,2000)); 

library(cowplot)
plot_hist4 <- plot_grid(hist1,hist2,hist3,hist4)


### Next steps: incorporate phenology model in generative approach, along with a non-linear relationship between
### canopy area (structure) and photosynthetic fraction (greenness). 



CA_trueMAT <- matrix(0,10^3,15); 
CA_obsMAT <- matrix(0,10^3,15); 
CA_trueGreen <- matrix(0,10^3,15); 


# Parameters for CA --> PSN:
int_green <- rnorm(10^3-naRows,0.507,0.008);
slope_green <- rnorm(10^3-naRows,0.302,0.016); 
err_green <- rnorm(10^3-naRows,0,0.123); 

# Parameters for phenology 


for(i in 1:(10^3-naRows)){
  CA_trueMAT[i,] <- ifelse(PPA_traj(paramMat.filtered[i,1],paramMat.filtered[i,2],paramMat.filtered[i,3])<1,
                           PPA_traj(paramMat.filtered[i,1],paramMat.filtered[i,2],paramMat.filtered[i,3]),1);
  for(j in 1:15){
    # CA_obsMAT[i,j] <- logit_link(x=rnorm(1,mean=CA_trueMAT[i,j], sd = CA_trueMAT[i,j]*0.2),a=-3,b=8); 
    CA_trueGreen[i,j] <- int_green[i] + slope_green[i]*CA_trueMAT[i,j] + err_green[i];
    CA_obsMAT[i,j] <- rbeta(1, shape1 = paramMat.filtered[i,4]*CA_trueGreen[i,j], shape2 = paramMat.filtered[i,4]*(1-CA_trueGreen[i,j])); 
  }
}







