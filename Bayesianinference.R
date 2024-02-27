# Let's define the two groups and their sample birth weights
nonsmoker <- c(7.5, 6.2, 6.9, 7.4, 9.2, 8.3, 7.6)
smoker <- c(6.2, 6.8, 5.7, 4.9, 6.2, 7.1, 5.9, 5.4)

# Question 1
# Marginal posterior of the mean is a scaled, shifted t-distribution

library(metRology)

mean_nonsmoker <- mean(nonsmoker)
mean_smoker <- mean(smoker)

n_nonsmoker <- length(nonsmoker)
n_smoker <- length(smoker)

sd_nonsmoker <- sqrt(var(nonsmoker))
sd_smoker <- sqrt(var(smoker))


mu <- seq(5, 10, 0.01)
likelihood_nonsmoker <- dt.scaled(mu, df = n_nonsmoker - 1, mean = mean_nonsmoker, sd = sd_nonsmoker/sqrt(n_nonsmoker))
likelihood_smoker <- dt.scaled(mu, df = n_smoker - 1, mean = mean_smoker, sd = sd_smoker/sqrt(n_smoker))


posterior_nonsmoker <- dt.scaled(mu, df = n_nonsmoker - 1, mean = mean_nonsmoker, sd = sd_nonsmoker/sqrt(n_nonsmoker))
posterior_smoker <- dt.scaled(mu, df = n_smoker - 1, mean = mean_smoker, sd = sd_smoker/sqrt(n_smoker))


plot(mu, likelihood_nonsmoker, type = 'l', col = 'blue')
lines(mu, posterior_nonsmoker, type = 'l', col = 'red')
# As the posterior and likelihood are the same, the two lines overlap


# The two posteriors in the same plot
plot(mu, posterior_smoker, type = 'l', col = 'blue', ylab = 'Density', 
     main = 'Posterior distribution of mean birth weight by group')
lines(mu, posterior_nonsmoker, type = 'l', col = 'red')
legend("topright", c('Smoker', 'Nonsmoker'), col = c("Blue", "Red"), lty = c(1,1))

v1 = n_nonsmoker - 1
v2 = n_smoker - 1

# Question 2


# The posterior mean is equal to the mean of the observed data (= posterior median and mode bc it's symmetrical and unimodal)
mean_nonsmoker
mean_smoker


# The posterior variance is v/(v-2) the variance of the observed data divided by the sample size due to the properties of the t-dist 
postsd1 <- sqrt((v1 / (v1 - 2)) * var(nonsmoker)/n_nonsmoker)
postsd2 <- sqrt((v2 / (v2 - 2)) * var(smoker)/n_smoker)

#Credible intervals
nonsmoker_lower <- mean_nonsmoker - qt(df = v1,0.975)*postsd1
nonsmoker_upper <- mean_nonsmoker + qt(df = v1,0.975)*postsd1

smoker_lower <- mean_smoker - qt(df = v2,0.975)*postsd2
smoker_upper <- mean_smoker + qt(df = v2,0.975)*postsd2

nonsmoker_int <- cbind(nonsmoker_lower, nonsmoker_upper)
nonsmoker_int

smoker_int <- cbind(smoker_lower, smoker_upper)
smoker_int


# Question 3

#Sampling
samplenonsmoker <- rt.scaled(100000, df = v1,mean_nonsmoker, sd_nonsmoker/sqrt(n_nonsmoker))
samplesmoker <- rt.scaled(100000, df = v2, mean_smoker, sd_smoker/sqrt(n_smoker))
samplediff <- samplenonsmoker - samplesmoker

ddiff = density(samplediff)
credint <- quantile(samplediff, c(0.025, 0.975))
plot(ddiff, main = 'Empirical density of difference in means', xlim = c(-1, 4))
abline(v = credint[1], col = "red")
abline(v = credint[2], col = "red")
legend("topleft", c('95% Credible interval'), col = c("Red"), lty = c(1))

mean(samplediff)
sqrt(var(samplediff))

#Credible interval for difference

diff_int <- cbind(credint[1], credint[2])
diff_int

#0 is not within the credible interval so we cannot conclude there is an association


# Question 4

library(readr)
library(coda)
library(runjags)
library(MCMCvis)
library(ggmcmc)
library(basicMCMCplots)
library(rjags)


Nchains <- 2

model.data <- list('nonsmoker' = nonsmoker, 'smoker' = smoker, 'N_nonsmoker' = n_nonsmoker,
                   'N_smoker' = n_smoker)

model.inits <- model.inits <- list(mu1 = 0, mu2 = 0, tau1 = 1, tau2 = 1)

cat("model
  {
    for (i in 1:N_nonsmoker){
      nonsmoker[i] ~ dnorm(mu1, tau1)
     
    }
    for (i in 1:N_smoker){
    smoker[i] ~ dnorm(mu2, tau2)
    }
    mu1 ~ dnorm(0, 1.0E-6)
    mu2 ~ dnorm(0, 1.0E-6)
    tau1 ~ dgamma(1.0E-3, 1.0E-3)
    tau2 ~ dgamma(1.0E-3, 1.0E-3)
    diff <- mu1 - mu2
  }", file="meandiff.txt")

jags <- jags.model('meandiff.txt',
                   data = model.data,
                   inits = model.inits,
                   n.chains = Nchains)

update(jags,2000) 
diff.sim <- coda.samples(jags,
                          c('diff'),
                          n.iter=8000, 
                          thin=1) 

print(diff.sim,digits=3)
plot(diff.sim) #Trace doesn't appear to contain a pattern, low autocorrelation which is good
diff.mcmc <- as.mcmc.list(diff.sim)


# Question 5
#Convergence tests

acfplot(diff.mcmc)
autocorr.plot(diff.mcmc) #Low autocorrelation suggests good convergence

gelman.diag(diff.mcmc)
gelman.plot(diff.mcmc,ask=FALSE) #Gelman-Rubin statistic is 1 -> good convergence

effectiveSize(diff.mcmc) #Effective size very high, = real size, good convergence again


# Question 6

plot(diff.mcmc)
summary(diff.mcmc)
HPDinterval(diff.mcmc) #0 still isn't in the interval, interval corresponds to the analytical one


# Question 7
#Can't initialize mu1 to 0 because that would mean dividing by 0 when calculating reldiff

model.inits2 <- list(mu1 = 0.1, mu2 = 0.1, tau1 = 1, tau2 = 1)

cat("model
  {
    for (i in 1:N_nonsmoker){
      nonsmoker[i] ~ dnorm(mu1, tau1)
     
    }
    for (i in 1:N_smoker){
    smoker[i] ~ dnorm(mu2, tau2)
    }
    mu1 ~ dnorm(0, 1.0E-6)
    mu2 ~ dnorm(0, 1.0E-6)
    tau1 ~ dgamma(1.0E-3, 1.0E-3)
    tau2 ~ dgamma(1.0E-3, 1.0E-3)
    reldiff <- (mu1 - mu2)/mu1
  }", file="meanreldiff.txt")

jags2 <- jags.model('meanreldiff.txt',
                   data = model.data,
                   inits = model.inits2,
                   n.chains = Nchains)

update(jags2,2000) 
reldiff.sim <- coda.samples(jags2,
                         c('reldiff'),
                         n.iter=8000, 
                         thin=1) 

print(reldiff.sim,digits=3)
plot(reldiff.sim)
reldiff.mcmc <- as.mcmc.list(reldiff.sim)

summary(reldiff.mcmc)
HPDinterval(reldiff.mcmc) 
#Mean relative difference according to MCMC is ~0.19-0.20
