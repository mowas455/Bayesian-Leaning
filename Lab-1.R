
#a.Draw Random number from the Posterior distribution
calc_mean = 0.366
calc_sd   = 0.0865

set.seed(12345)
Posterior_value<-function(n){
  theta_samples <- rbeta(n,11,19)
  return(theta_samples)
}

Graph <- function(s){
  no_breaks<-min(length(s)/2,30)
  hist(s,breaks=no_breaks,freq=FALSE,xlab="Samples",main="Histogram of Posterior Samples")
  lines(density(s))
  abline(v=mean(s),col="blue",lwd=3)
  abline(v=calc_mean,col="red",lwd=4)
}


P1<-Posterior_value(100)
P1_mean<-mean(P1)
P1_sd<-sd(P1)
Graph(P1)

P2<-Posterior_value(1000)
P2_mean<-mean(P2)
P2_sd<-sd(P2)
Graph(P2)

P3<-Posterior_value(10000)
P3_mean<-mean(P3)
P3_sd<-sd(P3)
Graph(P3)

#b.Use simulation from the Exact value
actual_Posterior_pb<- 1- pbeta(0.4,11,19)

sample_Posterior_pb<-sum(P3>0.4)/10000

cat("The difference between the actual probability and sample probability is",abs(actual_Posterior_pb-sample_Posterior_pb))


#c.Compute the log odds

log_values<-sapply(P3,function(y) log(y/(1-y)))
log_values

hist(log_values,breaks=30,freq=FALSE,xlab="Samples",main="Histogram of Posterior for Log Odd values")
lines(density(log_values))
abline(v=mean(log_values),col="blue",lwd=3)   

#==============================================================================================================
#2.
#a.Posterior variance




NormalNonInfoPrior <- function(NDraws,Data){
  Datamean <- 3.8
  s2 <- (sum(log(Data)-Datamean)^2)/length(Data)
  n <- length(Data)
  PostDraws <- matrix(0,NDraws,2)
  PostDraws[,1] <- ((n)*s2)/rchisq(NDraws,n)
  #PostDraws[,1] <- rnorm(NDraws,mean=Datamean,sd=sqrt(PostDraws[,2]/n))
  
  return(PostDraws)
}



Ndraws <- 10000
Data<-log(x)
PostDraws <- NormalNonInfoPrior(Ndraws,Data) # Generating draws from the joint posterior of mu and sigma^2
hist(PostDraws[,1],breaks = 30)
lines(density(PostDraws[,1]))
# Plotting the histogram of mu-draws
PostDraws

x<-c(38,20,49,58,31,70,18,56,25,78)
N<-10000
n<-length(x)
data<-log(x)



#mean value
set.seed(12345)
mean_sample <- mean(data)    # sample mean  
mu<-3.8                      # actual mean
cat("The mean of the sample is",mean_sample,"\n")


# Standard Deviation
sd_sample <- sum((data-mu)^2)/n
tau<-sd_sample


# Sample Variance
samples_chisq <- rchisq(N,n)
PostDraws <- tau * n / samples_chisq
#samples_mean <- mean(samples_invchisq_scaled)
#cat("The mean of inverse chi-square-values",samples_mean,"\n")


X <- seq(0.1, 2, length = N)
theo_var2 <- LaplacesDemon::dinvchisq (X, df = n, scale = tau)
hist(theo_var2,freq = F)
hist(PostDraws,freq = F)
lines(density(PostDraws))
abline(v=mean(PostDraws), col="red")

d <- density(samples_invchisq_scaled)
hist(x = samples_invchisq_scaled,y=density(samples_invchisq_scaled))
lines(samples_invchisq_scaled,main="Density of posterior values",xlab="sigma")
abline(v=samples_mean, col="red")

samples_invchisq_scaled

#======================
library(ggplot2)
G <- 2 * pnorm(sqrt(PostDraws/2)) - 1
ggplot(as.data.frame(G)) +
  geom_histogram(aes(x = G, y=..density..), bins = 40, fill = "#ffffffff", colour = "black", size = 0.2) +
  geom_density(aes(x = G, y=..density..), colour = "#DD141D", size = 0.5) +
  labs(subtitle = "Gini Index",
       y = "Density",
       x = "G") +
  theme_bw()


d<-density(G)
?density


q5 <- quantile(G,.05)
q95 <- quantile(G,.95)
G<-as.data.frame(G)
p <- ggplot(data = G) + theme_bw() + 
  geom_density(aes(x=G, y = ..density..), color = 'black')
x.dens <- density()
df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
p + geom_area(data = subset(df.dens, x >= q5 & x <= q95), 
              aes(x=x,y=y), fill = 'blue') +
  geom_vline(xintercept = medx)



#=============================
require(graphics)
K<-seq(0.01,10,0.01)
y<-c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)


prior_fn<-function(x,lamda=1)
{
  e<-dexp(x,rate = lamda)
  return(e)
  
}

likelihood_fn<-function(k,data,mu=2.39)
  {
  nume<-exp(k*sum(cos(data-mu)))
  deno<-(2*pi*besselI(k,0))**length(data)
  von_ll<-nume/deno
  return(von_ll)
  }

posterior_fn<-function(k,data,mu=2.39)
{
  nume<-exp(k*sum(cos(data-mu))-k)
  deno<-(besselI(k,0))**(length(data))
  p<-nume/deno
  return(p)
  
}

prior_fn(K)
likelihood_fn(K,y)
posterior_fn(K,y)
hist(prior_fn(K),freq = F)
lines(likelihood_fn(K,y))



d<- data.frame(k=K,v1=prior_fn(K),v2=likelihood_fn(K,y),v3=posterior_fn(K,y))
library(ggplot2);library(reshape2)
data<- melt(d)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()


ggplot(d)+
  geom_line(aes(x=k,y=v1),size=0.5)+
  geom_line(aes(x=k,y=v2),size=0.5)+
  geom_line(aes(x=k,y=v3),size=0.5)



prior_fn(K)
hist(prior_fn(K))
hist(likelihood_fn(K,y))



```{r}
library(ggplot2)
d<- data.frame(k=K,v1=prior_fn(K),v2=likelihood_fn(K,y),v3=posterior_fn(K,y))
p1<-ggplot(d)+ geom_line(aes(x=k,y=v1),size=0.5)
p2<-ggplot(d)+
  geom_line(aes(x=k,y=v2),size=0.5)
p3<-ggplot(d)+
  geom_line(aes(x=k,y=v3),size=0.5)
p1
p2
p3
```
