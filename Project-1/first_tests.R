library(readr)
library(ggplot2)
snow_particles <- read_csv("1_snow_particles.csv")
View(snow_particles)

breakpoints <- c(snow_particles$startpoint[1],snow_particles$endpoint)
breakpoints[length(breakpoints)] = breakpoints[length(breakpoints)-1]+0.25

n_particules <- snow_particles$particles.detected[1]
n_snow_bin <- round(0.01*snow_particles$particles.detected*snow_particles$retained....)


# Will not work due to the size of the bins 
# -> find something else
barplot(n_snow_bin, 
        names.arg = breakpoints[-length(breakpoints)], 
        width = 0.06,
        xlim = c(breakpoints[1], breakpoints[length(breakpoints)]), 
        ylim = c(0, max(n_snow_bin) * 1.1), 
        xlab = "Snowflakes diameter", 
        ylab = "Number of snowflakes",
        main = "Snowflake distribution",
        cex.axis = 0.75,cex.names = 0.7,las=2)

# Seems OK at first sight
# Uniform distrubtion on each bin with # snow particles
set.seed(42)
X = rep(0,n_particules)
count <- 0
for (i in 1:length(n_snow_bin)){
  if (n_snow_bin[i] != 0){
    start <- count+1
    end <- start+n_snow_bin[i]-1
    X[start:end] <- runif(n_snow_bin[i],snow_particles$startpoint[i],snow_particles$endpoint[i])
    count <- count + n_snow_bin[i]
  }
}
# check uniform dist
hist(X,breaks = 100,xlim = c(0.06,0.09))
