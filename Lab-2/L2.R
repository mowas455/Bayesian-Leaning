##### Q1 #####
library(mvtnorm)
data = read.table("TempLinkoping.txt", header = TRUE)

which.min(data$temp)

plot(x = data$time, y = data$temp, xlab="Time", ylab = "Temperature(Celcius)",
     main = "Time V. Temperature(Linkoping)", col = "red")

#a)
plot_data = function(data, v0,data_var_0, omega0, mu0,nDraws = 100){
  set.seed(2020-04-19)
  X = data.frame(c(1, data['time'], data['time']**2))
  n  = dim(X)[1]
  X_T = data.frame(t(X))
  XTX = as.matrix(X_T)%*% as.matrix(X)
  target = as.matrix(x = data['temp'], ncol=1)
  YTY = t(target) %*% target
  # Beta
  omega_o = omega0 * diag(3)
  omega_o_inv = solve(omega_o)
  omega_n = XTX + omega_o
  mu_o = matrix(data = mu0, ncol = 1)
  # variance
  vo = v0
  v_n = vo + n
  data_var_o = data_var_0
  nDraws = nDraws
  data_var = vo*data_var_o/rchisq(nDraws, df = vo)

  # Beta|sigma^2
  beta_prior = matrix(nrow = nDraws, ncol=3)
  for(draw in 1:nDraws){
    beta_prior[draw, ] = rmvnorm(1, mean = mu_o, sigma = data_var[draw]* omega_o_inv)
  }
  ### Prior Beta
  beta = as.matrix(beta_prior, ncol = 3)
  # Y_hat
  error = rnorm(nDraws, 0, sd = sqrt(data_var))
  y_hat = as.matrix(X) %*% t(beta) + error

  # Plot
  plot(y = data$temp, x=1:365,type = "l", pch = 1,  col = "red", ylim = c(-20,30),
       lwd = 3,
       xlab = "Time(Days)",ylab = "Temperature", main = "Time V. Temperature(Linkoping)")
  for(i in 1:nDraws){
    lines(x = 1:365, y = y_hat[,i], type="l",
          col = rgb(0.01,0.1,0.03,0.4), lwd = 0.5)
  }
  lines(x = 1:365, y = rowMeans(y_hat), type="l", col ="blue", lwd = 3)
  # legend("topleft",
  #        legend = c("Actual Temp", "Expected Temp", "Sample"), col = c("red", "blue", rgb(0.01,0.1,0.03,0.4)),
  #        lty = c(1,1,1))

}

# fig:2
par(mfrow=c(1,1))
plot_data(data = data, v0 = 4,data_var_0 = 0.28, omega0 = 0.01, mu0 = c(-10,100,-100), nDraws = 100)
# fig:3
plot_data(data = data, v0 = 10,data_var_0 = 0.1, omega0 = 1, mu0 = c(-15,120,-109), nDraws = 100)

## (b)
## Plot Posterior
plot_posterior_data = function(data, v0,data_var_0, omega0, mu0,nDraws = 1000){
  X = data.frame(c(1, data['time'], data['time']**2))
  n  = dim(X)[1] # No of observations
  D = dim(X)[2] # No of covariates
  X_T = data.frame(t(X))
  XTX = as.matrix(X_T)%*% as.matrix(X)
  target = as.matrix(x = data['temp'], ncol=1)
  YTY = t(target) %*% target
  # Beta
  omega_o = omega0 * diag(3)
  omega_o_inv = solve(omega_o)
  omega_n = XTX + omega_o
  mu_o = matrix(data = mu0, ncol = 1)
  # variance
  vo = v0
  v_n = vo + n
  data_var_o = data_var_0
  nDraws = nDraws
  data_var = vo*data_var_o/rchisq(nDraws, df = vo)

  # Beta|sigma^2
  beta_prior = matrix(nrow = nDraws, ncol=3)
  for(draw in 1:nDraws){
    beta_prior[draw, ] = rmvnorm(1, mean =  c(-10,100,-100), sigma = data_var[draw]* omega_o_inv)
  }
  ### Prior Beta
  beta = as.matrix(beta_prior, ncol = 3)
  # Y_hat
  error = rnorm(nDraws, 0, sd = sqrt(data_var))
  y_hat = as.matrix(X) %*% t(beta) + error

  # Calculating Parameters
  XTX_beta_hat =  as.matrix(X_T)%*%target
  var1 = solve(XTX + omega_o)
  var2 = XTX_beta_hat+(omega_o%*% mu_o)
  mu_n = var1%*%var2

  omega_n = XTX + omega_o
  omega_ni = solve(omega_n)

  v_n = vo + n

  sigma_n2 = ((vo*data_var_o) + t(target)%*%target + (t(mu_o)%*%omega_o%*%mu_o) - (t(mu_n)%*%omega_n%*%mu_n))/v_n
  sigma_n2 = sigma_n2[1]

  # sigma_2 posterior
  sigma2_post = matrix((v_n*sigma_n2)/rchisq(nDraws, df = v_n), ncol = 1)
  beta_dist = matrix(nrow = nDraws, ncol = D)

  # Beta Posterior
  for(i in 1:nDraws){
    beta_dist[i,] = rmvnorm(1, mean = mu_n, sigma = sigma2_post[1]*omega_ni)
  }

  dist_data = as.matrix(cbind(beta_dist, sigma2_post))

  par(mfrow = c(2,2))
  names = c("beta_0", "beta_1", "beta_2", "sigma_sq")
  for(i in 1:4){
    hist(x = dist_data[,i], xlab = names[i], breaks = 50,
         main = paste0("Histogram of parameter: ", names[i]), freq=F)
  }

  return(dist_data)
}
posterior_params = plot_posterior_data(data = data, v0 = 4, data_var_0 = 1, omega0 = 0.01, mu0 = c(-10,100,-100), nDraws = 1000)
## (b)(ii)
X = data.frame(c(1, data['time'], data['time']**2))
XT = t(X)

posterior_beta_draws = posterior_params[,1:3] # Draws of Betas
#Calculate f(time)
f_time  = posterior_beta_draws %*% XT
# Calculate median for each time-stamp's draw
median = apply(f_time, 2, median)
# Calculate 95% credible region for each time-stamp
prediction_intervals = apply(f_time, 2 , quantile, probs = c(0.025,0.975))
par(mfrow = c(1,1))
plot(y = data$temp, x=1:365 ,type = "b", pch = 1,  col = "red", ylim = c(-20,30),
     lwd = 3,
     xlab = "Time",ylab = "Temperature", main = "Time V. Temperature(Linkoping)")
for(i in 1:2){
  lines(x = 1:365, y = prediction_intervals[i,], type="l",
        col = rgb(0.1,0.2,0.03,0.4), lwd = 2)
}
lines(x = 1:365, y = median, col="blue", lwd=2)
legend("topleft",
       legend = c("Actual Temp", "Median Temp", "Credible Region"), col = c("red", "blue", rgb(0.01,0.1,0.03,0.4)),
       lty = c(1,1,1))

#(c)
# Drawing the samples of hottest days
d = density(apply(f_time, 1, which.max))
hist(apply(f_time, 1, which.max), freq = F, ylim = c(0,0.3),
     xlab="Day",col = "grey",
     main = "Histogram of Hottest day in each draw of f(time)") # Close to July 17
lines(x = d$x, y  = d$y, col="red")

# Analytically calculating the hottest days from f(time)
b1 = posterior_beta_draws[,2]
b2 = posterior_beta_draws[,3]
mean((-b1)/(2*b2))

which.min(data$time<=mean((-b1)/(2*b2)))

##### Q 2 #####
## (a)
library(mvtnorm) 
set.seed(2020-04-19)
WomenWorkData =read.table(file = "C:/Users/Mowniesh/OneDrive/Desktop/STIMA/4-Bayesian Learning/Lab-2/WomenWork.dat", header = 1)
summary(WomenWorkData)
covs= c(2:9) # Select which covariates/features to include
#standardize = TRUE # If TRUE, covariates/features are standardized to mean 0 and variance 1
tao = 10 # scaling factor for the prior of beta

Nobs = dim(WomenWorkData)[1]
y = WomenWorkData$Work # Target Vecotr

X = as.matrix(WomenWorkData[,covs]) # Co-variate Matrix
Xnames=colnames(X)

Npar = dim(X)[2] # No. of Parameters

# Setting up Prior
mu = as.matrix(rep(0,Npar))  # Prior mean vector
Sigma = 100*diag(Npar) # Prior covariance matrix

# Functions that returns the log posterior for the logistic regression.

logPostLogistic <- function(betas, y,X, mu,Sigma){
  set.seed(2020-04-19)
  linPred = X%*%betas
  logLik = sum( linPred*y - log(1+exp(linPred) ) )
  logPrior = dmvnorm(betas, mean = mu, Sigma, log = TRUE)
  return(logLik + logPrior)
}

# Initialize Betas
initVal = matrix(0,Npar,1)

OptimRes <- optim(initVal,logPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

#Printing the results to the screen
names(OptimRes$par) <- Xnames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimRes$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames # Naming the coefficient by covariates
print('The posterior mode is:')
print(OptimRes$par) ##### Numerical Values of Beta @ Mode
print('The approximate posterior standard deviation is:')
approxPostStd <- sqrt(diag(-solve(OptimRes$hessian)))
print(approxPostStd)
# Hessian Matrix
hessian = OptimRes$hessian
#solve(-hessian) == -solve(hessian)
# Neg Inv Hessian
posterior_var_mat = -solve(hessian)
# Beta draws
betas = rmvnorm(1000, mean = OptimRes$par, sigma = posterior_var_mat)
# 95% Equal Tail Posterior Probability for 'NSmalChild' variable
quantile(betas[,7], probs = c(0.025,0.975))
# Test Beta Vals
glm_model = glm(Work ~ 0 + ., data = WomenWorkData, family = binomial)

round(glm_model$coefficients - OptimRes$par, digits = 2)
#### b)

# Test Data
x <- matrix(c(1, 13, 8,11,(11/10)^2, 37,3,6), ncol = 1)

sim_draws<-function(X,beta){
  # For non-linear mapping
  sigmoid <- function(x){exp(x)/(1+exp(x))}
  XTB = beta%*%X
  p_y = apply(XTB, FUN = sigmoid, MARGIN = 2)
  return(p_y)
}

p_y = sim_draws(X=x, beta = betas)
d_py = density(p_y)
hist(p_y, freq = F, breaks = 100,
     xlab = "Probability",
     main = "Pr(y=1|x)")
lines(x=d_py$x, y= d_py$y, col="red")

### c)

# Expected P(y=1|x)
theta = mean(p_y)
hist(rbinom(1000, 8, p=theta), freq=F,
     main="", xlab="No of working women(Out of 8)")
