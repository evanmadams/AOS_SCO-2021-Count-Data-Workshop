#ubms tutorial
#install.packages("ubms")
library(ubms)
library(unmarked)
library(dplyr)
library(reshape2)

#We'll start by reading in the same data from the N-mix Intro.

nobo <- read.csv('data/nobo_abund.csv')
nobo.umf <- unmarkedFramePCount(y=nobo[,2:4], 
                                siteCovs=as.data.frame(scale(nobo[,5:10])), 
                                obsCovs=list(sky=scale(nobo[,11:13]),jdate=scale(nobo[,14:16]),
                                             time=scale(nobo[,17:19])))
summary(nobo.umf)

#N-mix function from unmarked
nobo.1 <- pcount(~sky + jdate + time ~BA + Evergreen5km, 
                 data=nobo.umf,K=105) 
summary(nobo.1)

#Now we'll run the same model with Stan via ubms package. 
# This requires a few extra specifications. 'chains' is
# the number of MCMC chains we want to run while 'iter'
# is the number of iterations (MCMC samples). We would
# typically run more iterations, especially since this
# gives us a warning indicating lack of convergence, 
# but the model does run more slowly in ubms/Stan.

options(mc.cores=3)

nobo.ubms = stan_pcount(~sky + jdate + time ~BA + Evergreen5km, 
                data=nobo.umf, K=105, chains=3, iter=100)

#Compare unmarked and Stan output. Results should be
# quite similar.
cbind(unmarked=coef(nobo.1), stan=coef(nobo.ubms))

#We can now expand on the data above with the 2018 data.
# We could run this as an open N-mixture model but instead
# are only interested in running it in the traditional
# N-mixture model but for both years. Because sites are
# sampled each year, we need to account for this pseudo-
# replication with a random effect of site.

nobo.full.dat <- read.csv("data/nobo_abund_full.csv")

#data frame is similar but with site covariate
nobo.full.umf <- unmarkedFramePCount(y=nobo.full.dat[,2:4], 
                  siteCovs=cbind(as.data.frame(scale(nobo.full.dat[,6:11])),site=nobo.full.dat[,1]), 
                  obsCovs=list(sky=scale(nobo.full.dat[,12:14]),
                           time=scale(nobo.full.dat[,15:17])))
summary(nobo.full.umf)


#The formula for the random effect is implemented 
# similar to our lmer code from earlier (1|site). We
# could implement a random effect for detection probability
# if we thought it was warranted for a given variable.
nobo.ubms.rand = stan_pcount(~sky ~BA + (1|site), 
                        data=nobo.full.umf, chains=3, iter=100)

#We've been using ubms but a few models (occu or pcount) in 
# unmarked support random effects now too (v. 1.1.1)

nobo.unmark.rand <- pcount(~sky ~BA + (1|site), data=nobo.full.umf)

nobo.ubms.rand
summary(nobo.unmark.rand)
sigma(nobo.unmark.rand) #extracts estimate of random effect

#This is another example this time with mallard data
# that is included in the unmarked package.
data("mallard")
str(mallard.y)
str(mallard.obs) #ivel is measure of survey effort
str(mallard.site)

umf.mall = unmarkedFramePCount(y=mallard.y, 
                        siteCovs=mallard.site,
                        obsCovs = mallard.obs)

#Notice that we're not specifying the mixture type (Poisson, NB) here
# Currently, the only option in ubms is for Poisson. But! We can
# account for a lot of overdispersion or extra variation with 
# random effects, so this is typically not a problem.

ubms.null = stan_pcount(~1 ~1, data=umf.mall, 
                        chains=3, iter=300)

ubms.elev <- stan_pcount(~1 ~elev, data=umf.mall,
                         chains=3, iter=300)

ubms.full = stan_pcount(~ivel + date ~elev + length + forest,
                      data=umf.mall, chains=3, iter=300)

#Model comparison in a Bayesian framework doesn't use
# AIC. In this case, we'll use leave-one-out cross-validation (LOO)
# But we still use the "fitList" and "modSel" functions

mods <- fitList(ubms.null, ubms.elev, ubms.full)
round(modSel(mods),3)

#In LOO, the larger value is the better model (ubms.full)

#ubms has several models from unmarked (but not all). We
# can also run distance sampling models with random effects.

dat <- dget('data/mebba_spp_Chipping Sparrow.dat')
y <- data.frame(site = 1:dat$nsites)
ds <- data.frame(site = dat$site.idx, dc = dat$distclass)
ds <- dcast(ds, site ~ dc, fun.aggregate = length)
y <- left_join(y, ds, 'site')
y[is.na(y)] <- 0
sc <- data.frame(dat$habcovs, dat$distcovs, dat$avcovs, dat$quadID)

#Note that our distance bins are different than unmarked
# We need the bins to be < 10 to help with convergence, so
# I changed the breaks and unitsIn to km
umf.ubms <- unmarkedFrameDS(y = as.matrix(y[1:100 ,2:6]), 
                siteCovs = sc[1:100,], 
                dist.breaks = c(0, 0.020, 0.040, 0.060, 0.080, 0.100), 
                survey = 'point', unitsIn = 'km') 

head(umf.ubms)

#from DistSampling_Intro.R
fm <- distsamp(~Wind.Speed.Code + Background.Noise.Code + ydate ~eff.domhabA, data = umf.ubms, keyfun = 'halfnorm', output = 'density', unitsOut = 'ha')

ubms.dist <- stan_distsamp(~Wind.Speed.Code + Background.Noise.Code + ydate ~eff.domhabA, 
                data=umf.ubms,
                keyfun = 'halfnorm',
                output = 'density',
                unitsOut = 'ha',
                chains = 3, iter = 1000)

cbind(unmarked=coef(fm), stan=coef(ubms.dist))

#Since we're working in ubms now, we can add random 
# effects. Here, we include a random effect for our 
# spatial unit. But note this isn't a spatial random effect.

dist.rand <- stan_distsamp(~Wind.Speed.Code + Background.Noise.Code + ydate ~eff.domhabA + (1|dat.quadID),
                data=umf.ubms,
                keyfun = 'halfnorm',
                output = 'density',
                unitsOut = 'ha',
                chains = 3, iter = 1000)

dist.rand

mods <- fitList(ubms.dist, dist.rand)
round(modSel(mods),3)
