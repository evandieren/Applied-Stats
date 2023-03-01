library(readr)
library(ggplot2)
snow_particles <- read_csv("1_snow_particles.csv")
n_bins <- length(snow_particles$startpoint)
snow_particles$endpoint[n_bins] <- snow_particles$endpoint[n_bins-1]+0.25

breakpoints <- c(snow_particles$startpoint[1],snow_particles$endpoint)

n_particules <- snow_particles$particles.detected[1]
n_snow_bin <- round(0.01*snow_particles$particles.detected*snow_particles$retained....)

width <- rep(0,length(n_snow_bin))
for (i in 1:length(n_snow_bin)){
  width[i] = snow_particles$endpoint[i]-snow_particles$startpoint[i]
}
bp <- barplot(height = n_snow_bin,
              width = width,
              xlim = c(breakpoints[1], breakpoints[length(breakpoints)]), 
              ylim = c(0, max(n_snow_bin) * 1.1), 
              xlab = "Snowflakes diameter", 
              ylab = "Number of snowflakes",
              main = "Snowflake distribution",
              cex.axis = 0.75,cex.names = 0.7,las=2)
axis(side=1, at =bp, labels=snow_particles$startpoint)

# Will not work due to the size of the bins 
# -> find something else
# barplot(n_snow_bin, 
#         names.arg = breakpoints[-length(breakpoints)], 
#         width = 0.06,
#         xlim = c(breakpoints[1], breakpoints[length(breakpoints)]), 
#         ylim = c(0, max(n_snow_bin) * 1.1), 
#         xlab = "Snowflakes diameter", 
#         ylab = "Number of snowflakes",
#         main = "Snowflake distribution",
#         cex.axis = 0.75,cex.names = 0.7,las=2)

# Seems OK at first sight
# Uniform distribution on each bin with # snow particles
set.seed(42)
X = rep(0,n_particules)
count <- 0
for (i in 1:length(n_snow_bin)){
  if (n_snow_bin[i] != 0){
    start <- count+1
    end <- start+n_snow_bin[i]-1
    X[start:end] <- runif(n_snow_bin[i],snow_particles$startpoint[i],snow_particles$endpoint[i])
    count <- count + n_snow_bin[i]
  }else{
    print("done")
    print(i)
    break
  }
}
# check uniform dist
hist(X,breaks = 100)

################
# EM - algorithm
################

dmixnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm(x,meanlog=mu1,sdlog=sigma1) + tau*dlnorm(x,meanlog=mu2,sdlog=sigma2)
  return(y)
}

mu_1_hat = -2; # mean from the biggest bin - bump #1
mu_2_hat = 0.1; # mean from the biggest bin - bump #2
sigma_1_hat = 0.3;
sigma_2_hat = 0.5;
tau_hat = 0.6;
N <- n_particules

p <- dlnorm(X,meanlog=mu_2_hat,sdlog=sigma_2_hat)*tau_hat/dmixnorm(X,mu_1_hat, mu_2_hat, sigma_1_hat, sigma_2_hat, tau_hat)

l_old <- Inf
l <- sum((1-tau_hat)*dlnorm(X,meanlog=mu_1_hat,sdlog=sigma_1_hat) + tau_hat*dlnorm(X,meanlog=mu_2_hat,sdlog=sigma_2_hat))

while (abs(l-l_old)>0.01){
  tau_hat <- sum(p)/N
  mu_1_hat <- sum((1-p)*log(X))/sum(1-p)
  mu_2_hat <- sum(p*log(X))/sum(p)
  sigma_1_hat <- sqrt(sum((1-p)*(log(X)-mu_1_hat)^2)/sum(1-p))
  sigma_2_hat <- sqrt(sum(p*(log(X)-mu_2_hat)^2)/sum(p))
  
  p <- dlnorm(X,meanlog=mu_2_hat,sdlog=sigma_2_hat)*tau_hat/dmixnorm(X, mu_1_hat, mu_2_hat, sigma_1_hat, sigma_2_hat, tau_hat)
  l_old <- l
  l <- sum((1-tau_hat)*dlnorm(X,meanlog=mu_1_hat,sdlog=sigma_1_hat) + tau_hat*dlnorm(X,meanlog=mu_2_hat,sdlog=sigma_2_hat))
}

init_vals <- c(mu_1_hat,mu_2_hat,sigma_1_hat,sigma_2_hat,tau_hat)

##############
# Optimization
##############

# Write a function given binned data

pmixnorm <- function(x,mu1,mu2,sigma1,sigma2,tau){
  y <- (1-tau)*plnorm(x,meanlog=mu1,sdlog=sigma1) + tau*plnorm(x,meanlog=mu2,sdlog=sigma2)
  return(y)
}
neg_loglikelihood <- function(par){
  # par = mu1,mu2,sigma1,sigma2,tau
  mu1<- par[1]
  mu2 <- par[2]
  sigma1 <- par[3]
  sigma2 <- par[4]
  tau <- par[5]
  y <- 0
  for (i in 1:n_bins) { # data = snow_particles
    y <- y + n_snow_bin[i]*log(pmixnorm(snow_particles$endpoint[i],mu1,mu2,sigma1,sigma2,tau)-pmixnorm(snow_particles$startpoint[i],mu1,mu2,sigma1,sigma2,tau))
  }
  return(-y)
}
out <- optim(init_vals,neg_loglikelihood)
expected_val_1 <- exp(out$par[1]+out$par[3]^2/2)
print(expected_val_1)
expected_val_2 <- exp(out$par[2]+out$par[4]^2/2)
print(expected_val_2)
std_var_1 <- sqrt((exp(out$par[3]^2)-1)*exp(2*out$par[1]+out$par[3]^2))
print(std_var_1)
std_var_2 <- sqrt((exp(out$par[4]^2)-1)*exp(2*out$par[2]+out$par[4]^2))
print(std_var_2)
tau_opt <- out$par[5]
print(tau_opt)

##########
# Plotting
##########

vals_hist <- unlist(hist(X,breaks = 100)["density"])
breaks_hist <- unlist(hist(X,breaks = 100)["breaks"])
vals <- sapply(breaks_hist[1:length(breaks_hist)-1],dmixnorm,
               mu1=out$par[1], mu2=out$par[2], 
               sigma1=out$par[3], sigma2=out$par[4], 
               tau=out$par[5])
# PLOT EM VALS pdf
plot(breaks_hist[1:length(breaks_hist)-1],vals_hist)
lines(x=breaks_hist[1:length(breaks_hist)-1],y=vals,type="l")
abline(v=mu1_opt)
abline(v=mu2_opt)


##################################
# LAST-STEP | Parametric bootstrap
##################################

#prob <- snow_particles$retained..../100
#X_class <- sample(x=1:length(n_snow_bin),size=N,replace = T,prob = prob)
#table(X_class)

# H_0 : diameters follow this bi-log normal distribution

# 1) Utiliser X jitter
# 2) Resample uniformly from X_jitter - X^b_jit
# 3) Rerun the EM algorithm and the optimization step on the X^b_jit 
# donc le f^b est la nouvelle phi binned avec le X^b_jit comparé à l'autre partie
# 4) On récupère nos estimateurs lambda_b hat

# 5) on calcule T*_b comme le suprémum et on le stock dans un array de taille B
# - Repeate 1 to 5 B times
# compute p-val and check whether it is high enough
