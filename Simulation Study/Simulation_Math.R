


### Spatial field simulation and development

# Dispersal Kernel 
PDFclark2dt<-function(x,u) { p <- 1;
disp <- 0.5*pi*x*(p/(pi*u))/(1+(x^2/u))^(p+1);
return(disp);} 

# Alternate kernel 

Kernel2 <- function(x,a,b){
  disp <- (a/(2*pi*gamma(2/b)))*exp(-(x^b/a^b));
  return(disp);
}

integrate(PDFclark2dt,lower = 0, upper = 300, u = 50)

integrate(Kernel2,lower = 0, upper = 300, a = 1, b = 2)


# Set the kernel to have u = 10, and then have a 15X15 grid 
# Each pixel will then have a "seed rain" computed from summing
# over all other pixels contributions 

# Q(j), seed source in pixel j, Seed_input(i) <- sum(Q(j)*k(i,j))
# This computation could be done by double-looping or something 

# We then sum discrete cohorts over time

# Simulation 1: uniform dmeographic parameters
# Simulation 2: spatially-varying demographic parameters 

# Show raster heat maps of the underlying parameters, and the resulting
# trends in reforestation 













