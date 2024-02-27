library(rjags)
library(coda)
library(tidyverse)
library(runjags)
library(gtools)

y=c(67,34,193,250,141)
n=c(282,225,290,261,141)
x=c(0,62.5,125,250,500)
iterations <- 10000
burnin <- 1000
chains <- 2

model_code <- "model
{
#likelihood
for (i in 1:5){
y[i] ~ dbinom(p[i],n[i])
logit(p[i]) = alpha + beta*x[i]
}
#priors
alpha ~ dunif(-100,100)
beta ~ dunif(-100,100)
}"

cat(model_code, file = "model.txt")
model.fit <- jags.model(file="model.txt",
                        data=list(n=n,y=y,x=x), n.chains = chains,
                        inits=list(list(alpha=-2,beta=0.1),
                                   list(alpha=4,beta=-0.1)))
model.samples <- coda.samples(model.fit, c("alpha", "beta"), n.iter=iterations)
summary(window(model.samples, start = burnin))
plot(model.samples, trace=TRUE, density = TRUE)
gelman.diag(model.samples,confidence=0.95)
combined.samples = combine.mcmc(model.samples)
HPDinterval(combined.samples,prob=0.95)
## Observed vs. fitted malformation probability plot
observed=data.frame(dose_obs=c(0,62.5,125,250,500),
                    p_obs=c(67/282,34/225,193/290,250/261,141/141))
fitted=data.frame(dose_fitted=0:500,p_fitted=exp(-1.79177+0.01832*dose_fitted)/
                    (1+exp(-1.79177+0.01832*dose_fitted)))
plot = ggplot()+geom_point(data=observed,aes(x=dose_obs,y=p_obs),size=3)+
  geom_line(data=fitted,aes(x=dose_fitted,y=p_fitted),lwd=1,linetype="longdash")+
  theme_bw()+xlab("DYME dose")+ylab("Malformation probability")+
  scale_x_continuous(expand = c(0, 0), limits = c(-5,505)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1.02))+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2),
         linetype=guide_legend(keywidth = 3, keyheight = 2))
plot
q = 0.05
alphas = c(model.samples[[1]][,1], model.samples[[2]][,1])
betas = c(model.samples[[1]][,2], model.samples[[2]][,2])
BMD = c()
# There's 20 000 samples after combining both chains
for(i in 1:20000) {
  P0 = exp(alphas[i])/(1+exp(alphas[i]))
  BMD = c(BMD, (logit(q*(1-P0)+P0)-alphas[i])/betas[i])
}
hist(BMD)
median(BMD) # posterior median = 17.11037
max(x) # tau can range from 0 to max(x) = 500
model_code <- "model
{
for (i in 1:5){
y[i] ~ dbinom(p[i],n[i])
logit(p[i]) = alpha + beta*(x[i]-tau)*(x[i]>tau)
}
alpha ~ dunif(-100,100)
beta ~ dunif(-100,100)
tau ~ dunif(0, 500)
}"

cat(model_code, file = "threshold_model.txt")
threshold.fit <- jags.model(file="threshold_model.txt",
                            data=list(n=n,y=y,x=x),
                            n.chains = chains,
                            inits=list(list(alpha=-2,beta=0.1,tau=100),
                                       list(alpha=4,beta=-0.1,tau=200)))
update(threshold.fit, burnin)
model.samples <- coda.samples(threshold.fit, c("alpha", "beta", "tau"), n.iter=iterations)
summary(window(model.samples, start = burnin))
par(mar = c(5, 5, 2, 2))
plot(model.samples, trace=TRUE, density = TRUE)
gelman.diag(model.samples,confidence=0.95)
combined.samples = combine.mcmc(model.samples)
HPDinterval(combined.samples,prob=0.95)
## Observed vs. fitted malformation probability plot
observed=data.frame(dose_obs=c(0,62.5,125,250,500),
                    p_obs=c(67/282,34/225,193/290,250/261,141/141))
dose_fitted=0:500
fitted=data.frame(dose_fitted=dose_fitted,p_fitted=exp(-1.2930+0.0272*((dose_fitted)-61.0357)*
                                                         ((dose_fitted)>61.0357))/
                    (1+exp(-1.2930+0.0272*((dose_fitted)-61.0357)*((dose_fitted)>61.0357))))
plot = ggplot()+geom_point(data=observed,aes(x=dose_obs,y=p_obs),size=3)+
  geom_line(data=fitted,aes(x=dose_fitted,y=p_fitted),lwd=1,linetype="longdash")+
  theme_bw()+xlab("DYME dose")+ylab("Malformation probability")+
  scale_x_continuous(expand = c(0, 0), limits = c(-5,505)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1.02))+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2),
         linetype=guide_legend(keywidth = 3, keyheight = 2))
plot
