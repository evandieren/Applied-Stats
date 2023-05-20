library(readr)

snow_particles <- read_csv("1_snow_particles.csv")
n_bins <- length(snow_particles$startpoint)
snow_particles$endpoint[n_bins] <- snow_particles$endpoint[n_bins-1]+0.25

breakpoints <- c(snow_particles$startpoint[1],snow_particles$endpoint)

n_particules <- snow_particles$particles.detected[1] # 700k
snow_particles$n_snow_bin = round(0.01*snow_particles$particles.detected*snow_particles$retained....)

width <- rep(0,n_bins)
for (i in 1:n_bins){
  width[i] = snow_particles$endpoint[i]-snow_particles$startpoint[i]
}

#png(file="./plots/snow_distribution.png",width=1200, height=700)
bp <- barplot(height = snow_particles$n_snow_bin/n_particules,
              width = width,
              xlim = c(breakpoints[1], breakpoints[length(breakpoints)]), 
              #ylim = c(0, max(snow_particles$n_snow_bin) * 1.1), 
              xlab = "Snowflakes diameter", 
              ylab = "Snowflakes distribution",
              main = "Snowflake distribution",
              cex.axis = 1,cex.names = 0.8)
axis(side=1, at =bp,label=snow_particles$startpoint,cex.axis = 1)
#dev.off()

# Jittering

set.seed(42)
X = rep(0,n_particules)
count <- 0
for (i in 1:n_bins){
  if (snow_particles$n_snow_bin[i] != 0){
    start <- count+1
    end <- start+snow_particles$n_snow_bin[i]-1
    X[start:end] <- runif(snow_particles$n_snow_bin[i],snow_particles$startpoint[i],snow_particles$endpoint[i])
    count <- count + snow_particles$n_snow_bin[i]
  }else{
    print("done")
    print(i)
    break
  }
}

# check uniform dist for jittering
hist(X,breaks = 100)

################
# EM - algorithm
################

dmixnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm(x,meanlog=mu1,sdlog=sigma1) + tau*dlnorm(x,meanlog=mu2,sdlog=sigma2)
  return(y)
}

em_alg <- function(data){
  mu_1_hat = -2; # mean from the biggest bin - bump #1
  mu_2_hat = 0.1; # mean from the biggest bin - bump #2
  sigma_1_hat = 0.3;
  sigma_2_hat = 0.5;
  tau_hat = 0.6;
  N <- n_particules
  p <- dlnorm(data,meanlog=mu_2_hat,sdlog=sigma_2_hat)*tau_hat/dmixnorm(data,mu_1_hat, mu_2_hat, sigma_1_hat, sigma_2_hat, tau_hat)
  l_old <- Inf
  l <- sum(log((1-tau_hat)*dlnorm(data,meanlog=mu_1_hat,sdlog=sigma_1_hat) + tau_hat*dlnorm(data,meanlog=mu_2_hat,sdlog=sigma_2_hat)))
  while (abs(l-l_old)>0.01){
    tau_hat <- sum(p)/N
    mu_1_hat <- sum((1-p)*log(data))/sum(1-p)
    mu_2_hat <- sum(p*log(data))/sum(p)
    sigma_1_hat <- sqrt(sum((1-p)*(log(data)-mu_1_hat)^2)/sum(1-p))
    sigma_2_hat <- sqrt(sum(p*(log(data)-mu_2_hat)^2)/sum(p))
    
    p <- dlnorm(data,meanlog=mu_2_hat,sdlog=sigma_2_hat)*tau_hat/dmixnorm(data, mu_1_hat, mu_2_hat, sigma_1_hat, sigma_2_hat, tau_hat)
    l_old <- l
    l <- sum(log((1-tau_hat)*dlnorm(data,meanlog=mu_1_hat,sdlog=sigma_1_hat) + tau_hat*dlnorm(data,meanlog=mu_2_hat,sdlog=sigma_2_hat)))
  }
  return(c(mu_1_hat,mu_2_hat,sigma_1_hat,sigma_2_hat,tau_hat))
}

init_vals <- em_alg(X)
init_vals
##############
# Optimization
##############

# Functions to get neg_likelihood

# cdf of bilog_likelihood
pmixnorm <- function(x,mu1,mu2,sigma1,sigma2,tau){
  y <- (1-tau)*plnorm(x,meanlog=mu1,sdlog=sigma1) + tau*plnorm(x,meanlog=mu2,sdlog=sigma2)
  return(y)
}

neg_loglikelihood <- function(par,data){
  # par = mu1,mu2,sigma1,sigma2,tau
  mu1<- par[1]
  mu2 <- par[2]
  sigma1 <- par[3]
  sigma2 <- par[4]
  tau <- par[5]
  y <- 0
  for (i in 1:n_bins) { # data = snow_particles
    y <- y + data$n_snow_bin[i]*log(pmixnorm(data$endpoint[i],mu1,mu2,sigma1,sigma2,tau)-pmixnorm(data$startpoint[i],mu1,mu2,sigma1,sigma2,tau))
  }
  return(-y)
}

# optimizing from EM points
theta_hat <- optim(init_vals,neg_loglikelihood,data = snow_particles)

theta_hat$par

##########
# Plotting
##########

vals_hist <- unlist(hist(X,breaks = 100)["density"])
breaks_hist <- unlist(hist(X,breaks = 100)["breaks"])
vals <- sapply(breaks_hist[1:length(breaks_hist)-1],dmixnorm,
               mu1=theta_hat$par[1], mu2=theta_hat$par[2], 
               sigma1=theta_hat$par[3], sigma2=theta_hat$par[4], 
               tau=theta_hat$par[5])

vals_em <- sapply(breaks_hist[1:length(breaks_hist)-1],dmixnorm,
                  mu1=init_vals[1], mu2=init_vals[2], 
                  sigma1=init_vals[3], sigma2=init_vals[4], 
                  tau=init_vals[5])
# PLOT EM VALS pdf
#png("./plots/EM_opt_dist.png",width=1200, height=700)
plot(breaks_hist[1:length(breaks_hist)-1],vals_hist,xlab="Snowflake diameters",ylab="Distribution")
lines(x=breaks_hist[1:length(breaks_hist)-1],y=vals,type="l")
lines(x=breaks_hist[1:length(breaks_hist)-1],y=vals_em,type="l",col=2)
legend("topright", legend = c("Optimal","EM-based"), col=1:2,lty=1) # optional legend
#dev.off()
##################################
# LAST-STEP | Parametric bootstrap
##################################

# H_0 : diameters follow this bi-log normal distribution

# 1) Utiliser X jitter and resample uniformly from X_jitter - X^b_jit
# 2) refaire le dataset pour avoir les nouveaux binnings -> avoir les start et end pareil mais juste changer les # elems
# 3) Rerun the EM algorithm and the optimization step on the X^b_jit 
# donc le f^b est la nouvelle phi binned avec le X^b_jit comparé à l'autre partie
# 4) On récupère nos estimateurs lambda_b hat

# 5) on calcule T*_b comme le suprémum et on le stock dans un array de taille B
# - Repeat 1 to 5 B times
# compute p-val and check whether it is high enough

# Helper functions
count_bins <- function(X, data){
  counts <- rep(0,n_bins)
  for (i in 1:n_bins){
    start <- data$startpoint[i]
    end <- data$endpoint[i]
    counts[i] <- sum(start <= X & X <= end)
  }
  return(counts)
}

B <- 20
data <- snow_particles[c("startpoint","endpoint","n_snow_bin")]
T_list <- numeric(B)
x_eval <- seq(from=0,to=snow_particles$endpoint[n_bins],length.out=1000) # creating breakpoints
Fn <- ecdf(X)
Fn_values <- Fn(x_eval) # values at left breakpoint
F_hat_values <- sapply(x_eval,pmixnorm,
                       mu1=theta_hat$par[1], mu2=theta_hat$par[2], 
                       sigma1=theta_hat$par[3], sigma2=theta_hat$par[4], 
                       tau=theta_hat$par[5])
abs_err <- abs(Fn_values-F_hat_values)
abs_err_max <- abs_err[which.max(abs_err)]
T_ <- abs_err_max
T_
for (b in 1:B) {
  # 1) Resampling from jittering
  X_b <- sample(X,n_particules,replace = T)
    
  # 2) Updating the binning dataset from new X_b
  counts <- count_bins(X_b,data)
  data$n_snow_bin = counts
    
  # 3) Re-running the em algorithm 
  init_val_X_b <- em_alg(X_b)
  print(b)
  
  # 4) Getting the best param estimators by optimizing the new neg_loglikelihood
  opt_X_b <- optim(init_val_X_b,neg_loglikelihood,data=data)
  # 5) Computing T*
  Fn <- ecdf(X_b)
  Fn_values <- Fn(x_eval) # values at left breakpoint
  F_hat_values <- sapply(x_eval,pmixnorm,
                    mu1=opt_X_b$par[1], mu2=opt_X_b$par[2], 
                    sigma1=opt_X_b$par[3], sigma2=opt_X_b$par[4], 
                    tau=opt_X_b$par[5])
  abs_err <- abs(Fn_values-F_hat_values)
  abs_err_max <- abs_err[which.max(abs_err)]
  T_list[b] <- abs_err_max
}

p_val = 1/(B+1)*(1+sum(T_<=T_list))
p_val

