
### -------------------------------------------------------------------------------------------------------- ###
#        This Script is part of the manuscript:
#
#
#        Blochmannia endosymbionts reduce brood rearing success in a carpenter ant (Camponotus sp.)
#
#        Anika Preuss, Peter Czuppon, Ulrich R. Ernst, Juergen Gadau
#
#
#        It was used to
#                  analyse the survival and development of brood data set
#                  generate figure S3
#
### -------------------------------------------------------------------------------------------------------- ###



###################################################
########################  Clear workspace
###################################################
rm(list = ls())
#### NEED TO CHANGE THE DIRECTORY !!!
setwd("/home/pete/Documents/Projects/01_Mini_projects/Blochmannia_Gadau/statsNew/")

###################################################
########################  Load libraries
###################################################
library(lme4)           # GLMMs
library(blme)           # Bayesian GLMMs
library(readxl)         # reading xls files
library(ggplot2)        # plotting
library(DHARMa)         # model fit diagnostics
library(lattice)        # visualizing random effects
library(sjPlot)         # output table of mixed effects model

###################################################
########################  Load dataset
###################################################
#d <- read_excel("Data_fostering_reversed.xlsx", sheet = "Whole Dataset")
d <- read.csv("Data_fostering_reversed.csv")
# replicate should be considered as a factor!
d$Replicate <- as.factor(d$Replicate)

# response data in the form [survived, dead] larvae
Survival_larvae <- cbind(d$Survived_larvae, d$Not_Survived_larvae)
Survival_pupae_corrected <- cbind(d$Survived_pupae_corrected, d$Not_Survived_pupae_corrected)

###################################################
###################################################
######################## STATISTICAL MODEL --  Eggs -> Larvae
###################################################
###################################################

########################
#### Model fitting 
########################
## Model: Y = beta0 (colony_ants, colony_eggs) + beta1 * Treatment_ants + beta2 * Treatment_eggs + beta3 * replicate + residual
# Note that treatment variable = 1 for no treatment and 0 for treatment!
# random intercept = mixed effects model: intercept depends on colony_ants and colony_eggs
# replicate is added as a fixed effect because of low number of replicates (random effects make sense for at least 3 replicates, because a variance needs to be estimated)
# binomial response variable -> link function = logit; logistic regression
# nAGQ = 1 gives Laplace approximation of the likelihood function
model <- glmer(Survival_larvae ~ Treatment_ants + Treatment_eggs + Replicate + (1|Colony_ants)  + (1|Colony_eggs), 
               data=d, family = binomial, nAGQ = 1, na.action = na.omit)

# because we have data separation (all larvae survived in one treatment), we run an additional estimation method to test the model fit
# method: penalized likelihood (Chung et al. (2013))
model_bayes <- bglmer(Survival_larvae ~ Treatment_ants + Treatment_eggs + Replicate + (1|Colony_ants) + (1|Colony_eggs), 
                     data=d, family = "binomial",
                     na.action = na.omit)

summary(model)
summary(model_bayes)
# -> both model fits yield the similar results for the (relevant) fixed effects

# In addition, we compute the confidence intervals (takes time!)
# bootstrap because Wald not good for probability estimates close to the boundary and profile did not converge
confint(model, method="boot", nsim=250)
confint(model_bayes, method="boot", nsim=250)
# -> comparing confidence intervals also yields no big differences for the fixed parameters
# -> we will continue with the model estimated by the glmer-method

###################################################
######################## DIAGNOSTICS
###################################################
########################
#### goodness of model fit -- assessment by simulation
########################
# simulation
sim <- simulate(model,nsim=100)

# data arrangement
dobs <- aggregate(d$Survived_larvae,by = list(AntTr = d$Treatment_ants, EggTr = d$Treatment_eggs),mean)
sim_p <- lapply(sim, aggregate,
              by = list(AntTr = d$Treatment_ants,EggTr = d$Treatment_eggs),
              mean)
sim_df <- do.call(rbind, sim_p)

# Relabeling the axes
sim_df$EggTr <- factor(sim_df$EggTr, levels = c("t","u"),
                    labels = c("treated eggs","untreated eggs"))
dobs$EggTr <- factor(dobs$EggTr, levels = c("t","u"),
                     labels = c("treated eggs","untreated eggs"))
sim_df$AntTr <- factor(sim_df$AntTr, levels = c("u","t"),
                     labels = c("untreated","treated"))
dobs$AntTr <- factor(dobs$AntTr, levels = c("u","t"),
                     labels = c("untreated","treated"))

# Plot
ggplot() +
  geom_jitter(aes(y = V1, x = AntTr), data = sim_df, 
              width = 0.2, height = 0, shape = 1, alpha = 1/2, grid=FALSE) +
  geom_point(aes(y = x, x = AntTr), data = dobs, 
             size = 4, color = "blue") +
  facet_wrap(~EggTr) +
  theme_classic() +
  labs(x = "Ant Treatment", y = "Predicted survived larvae")

# -> looks good! The statistical model seems to capture the experimental observation.

########################
#### evaluation of random effects
########################
dotplot(ranef(model))$Colony_eggs
dotplot(ranef(model))$Colony_ants

# -> random effects are spread, which indicates that they are relevant

# are they normally distributed?
qqnorm(ranef(model)$Colony_eggs[,1])
qqline(ranef(model)$Colony_eggs[,1])

qqnorm(ranef(model)$Colony_ants[,1])
qqline(ranef(model)$Colony_ants[,1])

# -> both qq plots look acceptable

########################
#### diagnostics with the DHARMa package
########################
### generate simulations
# unconditional random effects
sim_uncondRE <- simulateResiduals(fittedModel = model, n = 1000, plot = F)

### run multiple tests:
# qq-plot and scaled residuals
plot(sim_uncondRE)
# -> both qqplot and scaled residuals looks good!

# over/under-dispersion?
testDispersion(sim_uncondRE)
# -> looks good

# zero-inflation?
testZeroInflation(sim_uncondRE)
# -> looks good

###### -> GLMM fit seems to be OK! -> we report the results from the glmer-estimation

###### produce sjPlot table
tab_model(model)

###### compute confidence intervals of the predictions
## define data treatment vector
nd <- data.frame(Treatment_ants = c("t","t","t","t","u","u","u","u"),Treatment_eggs = c("t","t","u","u","t","t","u","u"),Replicate = c("1","2","1","2","1","2","1","2"))

# need to average over the two replicates, which is done in this function
confFunc <- function(m){
  # ignore the random effects -> only take fixed effects into account.
  v <- predict(m,newdata=nd,re.form=~0,type="response")
  vN <-seq(1,4)
  for (i in seq(1,4)){
    vN[i] <- (v[2*i-1]+v[2*i])/2
  }
  return(vN)
}

# bootstrap models and apply the function above, which produces predictions for 1000 bootstrapped samples
# note that random effects are estimated anew for each sample
a <- bootMer(model,confFunc,nsim=1000,parallel="multicore",ncpus=8)

# compute confidence intervals on the response scale
predCI <- apply(a$t,MARGIN=2,FUN=quantile,probs=c(0.025,0.975))

###################################################
###################################################
######################## STATISTICAL MODEL --  Larvae -> Pupae
###################################################
###################################################

########################
#### Model fitting
########################
## Model: Y = beta0 (colony_ants, colony_eggs) + beta1 * Treatment_ants + beta2 * Treatment_eggs + beta3 * replicate + residual
# Note that treatment variable = 1 for no treatment and 0 for treatment!
# random intercept = mixed effects model: intercept depends on colony_ants and colony_eggs
# replicate is added as a fixed effect because of low number of replicates (random effects make sense for at least 3 replicates, because a variance needs to be estimated)
# binomial response variable -> link function = logit; logistic regression
# nAGQ = 1 gives Laplace approximation of the likelihood function
model <- glmer(Survival_pupae_corrected ~ Treatment_ants + Treatment_eggs + Replicate + (1|Colony_ants)  + (1|Colony_eggs), 
               data=d, family = binomial, nAGQ = 1, na.action = na.omit)

summary(model)

### we see that the random effects colony_ants and colony_eggs have zero variance, which causes the warning "isSingular"
# to check for parameter stability, we run a Bayesian maximum likelihood
# method: penalized likelihood
model_bayes <- bglmer(Survival_pupae_corrected ~ Treatment_ants + Treatment_eggs + Replicate + (1|Colony_ants)  + (1|Colony_eggs), 
                      data=d, family = binomial, na.action = na.omit)


### An alternative is to simplify the model by removing the zero-variance random effects 
# suggested by Barr et al. (2013): Random effects structure for confirmatory hypothesis testing: Keep it maximal 
## Model: Y = beta0 + beta1 * Treatment_ants + beta2 * Treatment_eggs + beta3 * replicate + residual
# Note that treatment variable = 1 for no treatment and 0 for treatment!
# no variance detected in the mixed model in colony_ants and colony_eggs -> maybe overparameterization
model_simplified <- glm(Survival_pupae_corrected ~ Treatment_ants + Treatment_eggs + Replicate, 
                          data=d, family = binomial, na.action = na.omit)

summary(model)
summary(model_bayes)
summary(model_simplified)
# -> similar results for the fixed effects estimates for all three statistical model estimations

###################################################
######################## DIAGNOSTICS
###################################################
########################
#### goodness of model fit -- assessment by simulation
########################
# simulation
sim = simulate(model_simplified, nsim = 1000, weights = as.integer(weights(model_simplified)))
sim_bayes = simulate(model_bayes, nsim = 1000, weights = as.integer(weights(model_bayes)))

# data arrangement
dobs = aggregate(d$Survived_pupae_corrected,by = list(AntTr = d$Treatment_ants, EggTr = d$Treatment_eggs),mean)

# simplified model
sim_p <- lapply(sim, aggregate,
                by = list(AntTr = d$Treatment_ants,EggTr = d$Treatment_eggs),
                mean)
sim_df <- do.call(rbind, sim_p)

# Bayes model
sim_p2 <- lapply(sim_bayes, aggregate,
                 by = list(AntTr = d$Treatment_ants,EggTr = d$Treatment_eggs),
                 mean)
sim_df2 <- do.call(rbind, sim_p2)

# re-ordering and labeling x-axis
treat <- character()
for (i in 1:length(sim_df[,1])){
  if (sim_df$AntTr[i]=="t"){
    if (sim_df$EggTr[i]=="t"){
      treat <- c(treat,"tt")
    }
    else{
      treat <- c(treat,"tu")
    }
  }
  else{
    if (sim_df$EggTr[i]=="t"){
      treat <- c(treat,"ut")
    }
    else{
      treat <- c(treat,"uu")
    }
  }
}

treat2 <- character()
for (i in 1:length(sim_df2[,1])){
  if (sim_df2$AntTr[i]=="t"){
    if (sim_df2$EggTr[i]=="t"){
      treat2 <- c(treat2,"tt")
    }
    else{
      treat2 <- c(treat2,"tu")
    }
  }
  else{
    if (sim_df2$EggTr[i]=="t"){
      treat2 <- c(treat2,"ut")
    }
    else{
      treat2 <- c(treat2,"uu")
    }
  }
}

# Adding full treatment condition to the data sets
sim_df$Treatment <- treat
sim_df2$Treatment <- treat2
dobs$Treatment <- c("tt","ut","tu","uu")

# Relabeling the axes
sim_df$Treatment <- factor(sim_df$Treatment, levels = c("ut","tt","uu","tu"))
sim_df2$Treatment <- factor(sim_df2$Treatment, levels = c("ut","tt","uu","tu"))

# x axis label for publication
Label_x_axis <- c("untreated workers \n eggs from treated colonies", "treated workers \n  eggs from treated colonies", 
                  "untreated workers \n eggs from untreated colonies", "treated workers \n eggs from untreated colonies") # rename x-Labels

# Plot: simplified model
ggplot() +
  geom_jitter(aes(y = V1, x = Treatment), data = sim_df, 
              width = 0.2, height = 0, shape = 1, alpha = 1/2) +
  geom_point(aes(y = x, x = Treatment), data = dobs, 
             size = 4, color = "red") +
  scale_x_discrete(labels = Label_x_axis) + # rename x-labels
  theme_classic() +
  labs(x = "Treatment", y = "Predicted survived pupae")

# -> untreated eggs look OK!
# -> not so good fit of treated eggs!

# Plot: Bayes model
ggplot() +
  geom_jitter(aes(y = V1, x = Treatment), data = sim_df2, 
              width = 0.2, height = 0, shape = 1, alpha = 1/2) +
  geom_point(aes(y = x, x = Treatment), data = dobs, 
             size = 4, color = "red") +
  scale_x_discrete(labels = Label_x_axis) + # rename x-labels
  theme_classic() +
  labs(x = "Treatment", y = "Predicted survived pupae")

# -> untreated eggs look OK (better than for the simplified model, but still not good)!
# -> treated eggs seem to be problematic!

# -> in general blme-estimate looks better than the simplified model
# -> we proceed with Bayesian estimate only!

########################
#### diagnostics with the DHARMa package
########################
### generate simulations
sim_uncondRE <- simulateResiduals(fittedModel = model_bayes, n = 1000, 
                                  plot = F,
                                  weights = as.integer(weights(model_simplified)))

# we ignore the warning because bglmer has the same model type as glmer

### run multiple tests:
# qq-plot and scaled residuals
plot(sim_uncondRE)
# -> qqplot looks OK 
# -> scaled residuals are OK for the untreated eggs (right 4 columns)
# -> scaled residuals are not good for the treated eggs (left 4 columns)

# over/under-dispersion?
testDispersion(sim_uncondRE)
# -> looks good -> no dispersion

# zero-inflation?
testZeroInflation(sim_uncondRE)
# -> looks good -> no zero-inflation

#### We add an interaction term between the fixed effects to the model
## Model: Y = beta0 (colony_ants, colony_eggs) + beta1 * Treatment_ants + beta2 * Treatment_eggs + beta3 * replicate + beta3 * Treatment_ants * Treatment_eggs + residual  
# Note that treatment variable = 1 for no treatment and 0 for treatment!
model_bayes_int <- bglmer(Survival_pupae_corrected ~ Treatment_ants * Treatment_eggs + Replicate + (1|Colony_ants)  + (1|Colony_eggs), 
                          data=d, family = binomial, na.action = na.omit)

summary(model_bayes_int)

## Look at the random effects
dotplot(ranef(model_bayes_int))$Colony_eggs
dotplot(ranef(model_bayes_int))$Colony_ants
# -> visually confirms that the random effects do not add much to the model

## Rerun the DHARMa diagnostics to see if the scaled residuals improved
sim_uncondRE_int <- simulateResiduals(fittedModel = model_bayes_int, n = 1000, 
                                      plot = F,
                                      weights = as.integer(weights(model_simplified)))

# qq-plot and scaled residuals
plot(sim_uncondRE_int)
# -> the interaction term indeed seems to have improved the model
# -> We can quantify this by the AIC
AIC(model_bayes,model_bayes_int)
# -> almost 50 point decrease in AIC -> strong indication for the interaction model to be
# a better fit to the data.

# for completeness we also check for over/under-dispersion
testDispersion(sim_uncondRE_int)
# -> looks good -> no dispersion

# ... and zero-inflation?
testZeroInflation(sim_uncondRE)
###### -> GLMM fit seems to be OK for the model with interactions!

# ... lastly, we also run simulations based on this model with interaction and compare the data location to the overall pattern
sim_bayes_int = simulate(model_bayes_int, nsim = 1000, weights = as.integer(weights(model_bayes)))

# Bayes model data
sim_p3 <- lapply(sim_bayes_int, aggregate,
                 by = list(AntTr = d$Treatment_ants,EggTr = d$Treatment_eggs),
                 mean)
sim_df3 <- do.call(rbind, sim_p3)

# reordering and relabeling x-axis
treat3 <- character()
for (i in 1:length(sim_df3[,1])){
  if (sim_df3$AntTr[i]=="t"){
    if (sim_df3$EggTr[i]=="t"){
      treat3 <- c(treat3,"tt")
    }
    else{
      treat3 <- c(treat3,"tu")
    }
  }
  else{
    if (sim_df3$EggTr[i]=="t"){
      treat3 <- c(treat3,"ut")
    }
    else{
      treat3 <- c(treat3,"uu")
    }
  }
}

# Adding full treatment condition to the data sets
sim_df3$Treatment <- treat3

# Relabeling the axes
sim_df3$Treatment <- factor(sim_df$Treatment, levels = c("ut","tt","uu","tu"))

# Relabeling the axes
Label_x_axis <- c("untreated workers \n eggs from treated colonies", "treated workers \n  eggs from treated colonies", 
                  "untreated workers \n eggs from untreated colonies", "treated workers \n eggs from untreated colonies") # rename x-Labels

# Plot: Bayes with interaction model
ggplot() +
  geom_jitter(aes(y = V1, x = Treatment), data = sim_df3, 
              width = 0.2, height = 0, shape = 1, alpha = 1/2) +
  geom_point(aes(y = x, x = Treatment), data = dobs, 
             size = 4, color = "red") +
  scale_x_discrete(labels = Label_x_axis) + # rename x-labels
  theme_classic() +
  labs(x = "Treatment", y = "Predicted survived pupae")

# looks very reasonable! 
# based on these diagnostics we take the model with interactions as our explanatory model

######### produce sjPlot table
tab_model(model_bayes_int)

###### compute confidence intervals of the predictions (bootstrap)
## define data treatment vector
nd <- data.frame(Treatment_ants = c("t","t","t","t","u","u","u","u"),Treatment_eggs = c("t","t","u","u","t","t","u","u"),Replicate = c("1","2","1","2","1","2","1","2"))

# need to average over the two replicates, which is done in this function
confFunc <- function(m){
  # ignore the random effects -> only take fixed effects into account.
  v <- predict(m,newdata=nd,re.form=~0,type="response")
  vN <-seq(1,4)
  for (i in seq(1,4)){
    vN[i] <- (v[2*i-1]+v[2*i])/2
  }
  return(vN)
}

# bootstrap models and apply the function above, which produces predictions for 1000 bootstrapped samples
# note that random effects are estimated anew for each sample
a2 <- bootMer(model_bayes_int,confFunc,nsim=1000,parallel="multicore",ncpus=8)

# compute confidence intervals on the response scale
predCI2 <- apply(a2$t,MARGIN=2,FUN=quantile,probs=c(0.025,0.975))
