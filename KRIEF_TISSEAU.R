#rm(list=ls())

# Question 1

simulate_brownian=function(gaussian_vector){
  n <- length(gaussian_vector)
  cumulative_sum <- c(0, cumsum(gaussian_vector))
  W <- (1/sqrt(n)) * cumulative_sum
  return(W)
}

simulate_trajectories=function(n_trajectories, n, graph = FALSE,  plot_cone = TRUE){
  #gaussian_matrix <- replicate(n_trajectories, rnorm(n, 0,1))
  gaussian_matrix <- matrix(rnorm(n * n_trajectories, 0, 1), nrow = n)
  brownian_trajectories <- apply(gaussian_matrix, 2, simulate_brownian)
  if (graph == TRUE){
    colors <- rainbow(n_trajectories)
    time_points <- seq(0, 1, length.out = n + 1)
    
    matplot(seq(0, 1, length.out = n + 1), brownian_trajectories, type = "l", xlab = "Time", ylab = "Brownian Motion Value",
            main = "Brownian Motion Trajectories", col = colors, ylim = c(-3, 3))
    upper_cone <- 1.96 * sqrt(time_points)
    lower_cone <- -upper_cone
    lines(time_points, upper_cone, col = "black", lty = 2, lwd = 2)
    lines(time_points, lower_cone, col = "black", lty = 2, lwd = 2)
    legend("topright", legend = c("±1.96 * sqrt(t)"), col = "black", lty = 2, lwd = 2)
    
  }
}

# Example of Brownian trajectories simulation

n <- 10000
n_trajectories <- 100
t=proc.time()
simulate_trajectories(n_trajectories, n, TRUE)
proc.time()-t

simulate_risky_asset=function(brownians, r, sigma){
  n <- length(brownians)
  index <- 2:n
  S <- exp((r-(sigma^2)/2)*index/n + sigma*brownians[index])
  return(S)
}

simulate_mutiple_risky_assets <- function(r, sigma, brownian_matrix) {
  S <- apply(brownian_matrix, 2, function(column) {
    monte_carlo_init <- simulate_risky_asset(column, r, sigma)
  })
  return(S)
}

positive_part <- function(vector) {
  vector[vector < 0] <- 0
  return(vector)
}

classic_MC_estimator=function(S, K) {
  prices <- apply(S, 2, function(column) {
    monte_carlo_init <- mean(column) - K
    positive_part(monte_carlo_init)
  })
  return(c(mean(prices),var(prices)/length(prices))) 
}

# Example of implementation to calculate G

# r <- 0
# sigma <- 1/4
# K <- 1
# n_trajectories <- 10000
# n <- 100
# t=proc.time()
# 
# brownian_matrix <- simulate_trajectories(n_trajectories, n)
# S <- simulate_mutiple_risky_assets(r, sigma, brownian_matrix)
# 
# G_estimate <- classic_MC_estimator(S, K)
# print(G_estimate)
# proc.time()-t

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Question 2 

antithetic_estimator=function(brownian_matrix, K, r, sigma){
  results_monte_carlo <- apply(brownian_matrix, 2, function(column) {
    monte_carlo_init_one <- mean(simulate_risky_asset(column, r, sigma)) - K
    monte_carlo_init_two <- mean(simulate_risky_asset(-column, r, sigma)) - K
    (positive_part(monte_carlo_init_one) + positive_part(monte_carlo_init_two))/2
  })
  return(c(mean(results_monte_carlo),var(results_monte_carlo)/length(results_monte_carlo))) #/length(results_monte_carlo)
}

# Example of implementation

# r <- 0
# sigma <- 1/4
# K <- 1
# n_trajectories <- 10000
# n <- 100
# t=proc.time()
# 
# brownian_matrix <- simulate_trajectories(n_trajectories, n)
# 
# G_estimate_bis <- antithetic_estimator(brownian_matrix, K, r, sigma)
# print(G_estimate_bis)
# proc.time()-t

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Question 3

control_variable=function(prices, K){
  n <- length(prices)
  return((prod(prices))^(1/n) - K)
}

simulate_control_variables <- function(K, S) {
  control_variables <- apply(S, 2, function(column) {
    monte_carlo_init <- control_variable(column, K)
    positive_part(monte_carlo_init)
  })
  return(remove_Inf_elements(control_variables))
}

remove_Inf_elements=function(prices){
  index_Inf <- which(prices == Inf)
  if (length(index_Inf) > 0) {
    prices <- prices[-index_Inf]
  }
  return(prices)
}

simulate_Gs <- function(K, S) {
  Gs <- apply(S, 2, function(column) {
    monte_carlo_init <- mean(column) - K
    positive_part(monte_carlo_init)
  })
  return(Gs)
}

estimate_bstar=function(Gs, control_variables){
  numerator <- sum((Gs - mean(Gs))*(control_variables - mean(control_variables)))
  denominator <- sum((control_variables - mean(control_variables))^2)
  b <- -numerator/denominator
  return(b)
}

expectation_control_variable=function(n){
  alpha <- sqrt(((n+1)*(2*n + 1))/6)
  sigma <- 1/4
  return(exp(-(sigma^2)*(n+1)/(4*n) + ((alpha^2)*sigma^2)/(2*n^2))*pnorm(-sigma*(n+1)/(4*alpha) + (sigma*alpha)/n) - (1 - pnorm(sigma*(n+1)/(4*alpha))))
}

control_variable_estimator=function(S, K, n_zero, n){
  control_variables <- simulate_control_variables(K, S)
  Gs <- simulate_Gs(K, S)[1:length(control_variables)]
  #print(cor(Gs, control_variables))
  b <- estimate_bstar(Gs[1:n_zero], control_variables[1:n_zero])
  estimate <- Gs + b*(control_variables - expectation_control_variable(n)) #mean(control_variables)
  return(c(mean(estimate),var(estimate)/length(estimate))) #/length(estimate)
}

# Example of implementation

# r <- 0
# sigma <- 1/4
# K <- 1
# n_trajectories <- 10000
# n <- 100
# n_zero <- 200
# t=proc.time()
# 
# brownian_matrix <- simulate_trajectories(n_trajectories, n)
# S <- simulate_mutiple_risky_assets(r, sigma, brownian_matrix)
# 
# G_estimate_ter <- control_variable_estimator(S, K, n_zero, n)
# print(G_estimate_ter)
# proc.time()-t

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Recap 

r <- 0
sigma <- 1/4
K <- 1
n_trajectories <- 10000
n <- 1000
n_zero <- 200
t=proc.time()

brownian_matrix <- simulate_trajectories(n_trajectories, n)
S <- simulate_mutiple_risky_assets(r, sigma, brownian_matrix)

#Estimations
G_estimate <- classic_MC_estimator(S, K)
G_estimate_bis <- antithetic_estimator(brownian_matrix, K, r, sigma)
G_estimate_ter <- control_variable_estimator(S, K, n_zero, n)

#Results
print(G_estimate)
print(G_estimate_bis)
print(G_estimate_ter)

#Gains
((G_estimate_bis[2] - G_estimate[2])/G_estimate[2])*100
((G_estimate_ter[2] - G_estimate[2])/G_estimate[2])*100

proc.time()-t

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Question 4

# Function to calculate z(y) and S(y)
zS_vectors <- function(y, n, sigma = 1/4, K = 1) {
  d_t = 1/n
  z_vector <- numeric(n)
  z_vector[1] <- (sigma*sqrt(d_t)*(y+K))/y
  
  S_vector <- numeric(n)
  S_vector[1] <- exp((-(d_t*sigma^2))/2 + sigma*sqrt(d_t)*z_vector[1])
  
  for (j in 2:n) {
    z_vector[j] <- z_vector[j-1] - (sigma*sqrt(d_t)*S_vector[j-1])/(n*y)
    S_vector[j] <- exp((-(d_t*sigma^2))/2 + sigma*sqrt(d_t)*sum(z_vector[1:j]))
  }
  return(c(mean(S_vector), z_vector))
}

#Bisection method
bisection_method <- function(n, sigma, K, tol = 1e-6, max_iter = 100) {
  a <- 0.01  # Lower bound
  b <- 0.5 # Upper bound 
  
  # Bisection algorithm
  for (iter in 1:max_iter) {
    c <- (a + b) / 2
    S_c <- zS_vectors(c, n, sigma, K)[1] - K - c
    
    # Check if the solution is found
    if (abs(S_c) < tol) {
      cat("Converged after", iter, "iterations.\n")
      return(c)
    }
    
    # Update the interval
    if (S_c * (zS_vectors(a, n, sigma, K)[1] - K - a) < 0) {
      b <- c
    } else {
      a <- c
    }
  }
  
  # Return the best guess if the maximum number of iterations is reached
  warning("Maximum number of iterations reached.")
  return(c)
}

custom_function <- function(col1, col2, mu, K = 1) {
  monte_carlo_init <- mean(col1) - K
  result <- positive_part(monte_carlo_init) * exp(-(mu %*% col2)-((mu%*%mu)/2))
  return(result)
}

IS_MC_estimator=function(n_trajectories, n, mu, K = 1, sigma = 1/4) {
  gaussian_matrix_IS <- matrix(rnorm(n * n_trajectories, 0, 1), nrow = n)
  brownian_matrix_IS <- apply(gaussian_matrix_IS + mu, 2, simulate_brownian)

  S <- simulate_mutiple_risky_assets(0, sigma, brownian_matrix_IS)
  
  IS <- numeric(n_trajectories)
  for (i in 1:n_trajectories){
    IS[i] = custom_function(S[, i], gaussian_matrix_IS[, i], mu)
  }
  
  return(c(mean(IS),var(IS)/length(IS)))
}

# Example of implementation

#Choice of the interval used in bisection method
K <- 1
sigma <- 1/4
n <- 100
n_trajectories <- 10000
t=proc.time()
y_values <- seq(-0.5, 0.5, length.out = 100)
z_values <- sapply(y_values, function(y) {zS_vectors(y, n = 100, sigma = 0.25, K = 1)[1] - y - 1})
plot(y_values, z_values, type = "l", col = "blue", lwd = 2, xlab = "y", ylab = "zS_vectors(y)", main = "Plot of zS_vectors(y)")

abline(h = 0, col = "gray", lty = 2)
# Answer : we take the interval [0.01, 0.5]. a=0.01, b=0.5

# Importance sampling

# - Finding the solution of the equation and computing mu
sol <- bisection_method(n, sigma, K = 1)
mu <- zS_vectors(sol, n, sigma, K)[2:(n+1)]
#cat("Solution:", mu, "\n")

# - Estimating G using the importance sampling method

G_estimate_IS <- IS_MC_estimator(n_trajectories, n, mu, K)
print(G_estimate_IS)
proc.time()-t

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#Plots

r <- 0
sigma <- 1/4
K <- 1
n <- 100
n_zero <- 200
trajectory_steps <- seq(500, 10000, by=500)

# Initialize vectors to store the variances
variances_classic <- numeric(length(trajectory_steps))
variances_antithetic <- numeric(length(trajectory_steps))
variances_control <- numeric(length(trajectory_steps))
variances_IS <- numeric(length(trajectory_steps))

for (i in seq_along(trajectory_steps)) {
  n_trajectories <- trajectory_steps[i]
  
  brownian_matrix <- simulate_trajectories(n_trajectories, n)
  S <- simulate_mutiple_risky_assets(r, sigma, brownian_matrix)
  
  # Classic Monte Carlo estimator
  G_estimate <- classic_MC_estimator(S, K)
  variances_classic[i] <- G_estimate[2] 
  
  # Antithetic estimator
  G_estimate_bis <- antithetic_estimator(brownian_matrix, K, r, sigma)
  variances_antithetic[i] <- G_estimate_bis[2] 
  
  # Control variable estimator
  G_estimate_ter <- control_variable_estimator(S, K, n_zero, n)
  variances_control[i] <- G_estimate_ter[2] 
  
  # Importance sampling estimator
  G_estimate_IS <- IS_MC_estimator(n_trajectories, n, mu, K)
  variances_IS[i] <- G_estimate_IS[2] 
}

# Plotting the convergence
plot(trajectory_steps, variances_classic, type="l", col="red", xlab="Number of Trajectories", ylab="Variance of Estimator", main="Convergence of Estimators", ylim=c(min(c(variances_classic, variances_antithetic, variances_control, variances_IS)), 
                                                                                                                                                                     max(c(variances_classic, variances_antithetic, variances_control, variances_IS))))
lines(trajectory_steps, variances_antithetic, col="blue")
lines(trajectory_steps, variances_control, col="green")
lines(trajectory_steps, variances_IS, col="purple")

legend("topright", legend=c("Classic MC", "Antithetic", "Control Variable", "Importance Sampling"), col=c("red", "blue", "green", "purple"), lty=1)






