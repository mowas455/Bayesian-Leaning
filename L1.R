library(ggplot2)

##### 1
## (a)
set.seed(2020-04-10)
s = 8; n=24;f=n-s
alpha = 3; beta=3
### Theoretical mu
e_theta = (alpha+s)/(alpha+beta+n)
### theoretical standard deviation
sqrt((e_theta*(1- e_theta))/(alpha+beta+n+1))
# Posterior Draws
bernoulli.posterior<- function(nDraws, s, f, a, b){
  shape1 = a+s
  shape2 = b+f

  dist = rbeta(n = nDraws, shape1 = shape1,shape2 = shape2)
  posterior_mean = mean(dist)
  posterior_variance = var(dist)
  return(c(posterior_mean,posterior_variance))
}

draws = seq(10,10000, by = 10)
length(draws)
result_matrix = matrix(0,2,length(draws))

result_matrix = sapply(draws, FUN = bernoulli.posterior, s=s,f=f,a=alpha,b=beta)

plot(result_matrix[1,], col=2, x = draws, xlab="# Draws", ylab = "Posterior Mean")

plot(sqrt(result_matrix[2,]), col=3, x = draws, xlab="# Draws", ylab = "Posterior St. Dev.")
## (b)
set.seed(2020-04-10)
simulation<- function(nDraws, s, f, a, b){
  shape1 = a+s
  shape2 = b+f

  dist = rbeta(n = nDraws, shape1 = shape1,shape2 = shape2)
  return(dist)
}
nDraws = 10000
result_sim = simulation(nDraws = nDraws, s=s,f=f,a=alpha,b=beta)

length(result_sim[result_sim>0.4])/nDraws # Analytical calculation

1 - pbeta(0.4, shape1 = alpha+s, shape2 = beta+f) # Actual Probability

## (c)
set.seed(2020-04-10)
log.odds<-function(nDraws, s, f, a, b){
  shape1 = a+s
  shape2 = b+f

  dist = rbeta(n = nDraws, shape1 = shape1,shape2 = shape2)
  odds = dist/(1-dist)
  log_odds = sapply(odds, FUN = log)
  return(log_odds)
}

res_log_odds = log.odds(nDraws = 10000, s = s,f = f, a=alpha, b= beta)

# plot
hist(res_log_odds, breaks=100, ylim=c(0,1), freq=F, main="Posterior Distribution of Log Odds", xlab="", ylab="Posterior Distribution(Logit)")
lines(density(res_log_odds), col = "blue")

##### 2 ####
data = c(38,20,49,58,31,70,18,56,25,78)
## (a)
log.norn<-function(y,mu,sd){
  v = sd**2
  return((1/(y*sqrt(2*pi*v)))*exp(-0.5/v*((log(y) - mu)**2)))
}

log_data = log(data)
n = length(data)
mu = 3.8
tao_sq = sum((log_data - mu)**2)/n
# Simulating from Posterior
# Lecture-3, pg. 5
set.seed(2020-04-10)
sigma_sq = (n*tao_sq)/rchisq(10000,n)
hist(sigma_sq, breaks = 1000, xlim = c(0,2))


mean(sigma_sq) # Simulated Mean
(10*tao_sq)/8 ## Theoretical mu (BDA pg. 577)
var(sigma_sq) # Simulated Variance
(2*(10*tao_sq*10*tao_sq))/(64*6) # Theoretical Variance (BDA pg. 577)

## (b)
#sqrt(2)
gini_dist = (pnorm(sqrt(sigma_sq/2),mean = 0, sd = 1)*2)-1
hist(gini_dist, freq = F,breaks=100, xlab="Gini Index", ylab= "Posterior Gini Distribution(Data)",
     main = "")
abline(v=0.26,col="red",lwd=3)
lines(density(gini_dist), col="blue", lwd=3)

## (c)
# 90% Credible Interval
ci_90 = quantile(gini_dist, probs = c(0.05,0.95))


densities = density(gini_dist, kernel = "gaussian")



dat <- with(densities, data.frame(x, y))
ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_line()+
  geom_area(mapping = aes(x = ifelse(x>ci_90['5%'] & x< ci_90['95%'] , x, 0)), fill = "red")+
  ylim(c(0,10)) +
  labs(x = "Gini Coefficient", y = "Posterior Density",
       title = "90% Credible Interval",subtitle = "Gini Coefficient")+
  geom_vline(xintercept = ci_90['5%'], color = "black", size=1, linetype="dotted")+
  geom_vline(xintercept = ci_90['95%'], color = "black", size=1,linetype="dotted")

## Highest Posterior Density Interval
# Reference
#https://stats.stackexchange.com/questions/381520/how-can-i-estimate-the-highest-posterior-density-interval-from-a-set-of-x-y-valu

hdi = function(x,y,coverage){
  l_x = length(x)
  best = 0
  for(ai in 1:(l_x-1))
  {
    for(bi in (ai +1): l_x){
      mass = sum(diff(x[ai:bi]) * y[(ai+1):bi])
      if (mass >= coverage && mass/(x[bi] - x[ai]) > best)
      {
        best = mass / (x[bi] - x[ai])
        ai.best = ai
        bi.best = bi
      }
    }
  }
  c(x[ai.best], x[bi.best])
}

hdci_90 = hdi(x = densities$x, y = densities$y, coverage = 0.9)
# Plot HPDI
ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_line()+
  geom_area(mapping = aes(x = ifelse(x>hdci_90[1] & x< hdci_90[2] , x, 0)), fill = "red")+
  ylim(c(0,10)) +
  labs(x = "Gini Coefficient", y = "Posterior Density",
       title = "90% Highest Posterior Density Interval",subtitle = "Gini Coefficient") +
  geom_vline(xintercept = hdci_90[1], color = "black", size=1, linetype="dotted") +
  geom_vline(xintercept = hdci_90[2], color = "black", size=1, linetype="dotted")

##### 3

## (a)
set.seed(2020-04-10)
data_radian = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

mu = 2.39 # Given
n = length(data_radian)
nDraws = 10000
prior_kappa = rexp(nDraws, rate = 1) # Draw Prior Kappa values

const2 = sum(cos(data_radian-mu)) # y_i - mu

# Calculate Posterior
posterior = exp(prior_kappa* (const2 - 1))/ (2*pi*besselI(x =prior_kappa, nu = 0))^n
post = posterior/sum(posterior)# Normalize

plot(x = prior_kappa, y = post, col="red", xlab = "kappa", ylab = "Posterior Density")

## (b)
prior_kappa[which.max(post)]

