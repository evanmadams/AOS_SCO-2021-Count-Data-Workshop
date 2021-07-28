# More advanced tools for N-mixture models

## Need to run Nmix Intro.R for complete code to run! ##


#### DISTRIBUTIONS FOR ABUNDANCE ####
# Can use zero-inflated Poisson (ZIP) or negative binomial (NB) as well

# Example with previous data
nobo.nb <- pcount(~sky ~BA, data=nobo.umf, mixture="NB")

summary(nobo.nb)

nobo.zip <- pcount(~sky ~BA, data=nobo.umf, mixture="ZIP")
summary(nobo.zip)

nobo.pois <- pcount(~sky ~BA, data=nobo.umf, mixture="P")

nobo.null <- pcount(~1 ~1, data=nobo.umf, mixture="P")

#### MODEL COMPARISON ####
## We can also compare between different data models
cbind(AIC.zip=nobo.zip@AIC, AIC.nb=nobo.nb@AIC, AIC.pois=nobo.pois@AIC)

## We can compare models with different covariates using AIC
cbind(AIC.1=nobo.1@AIC, AIC.pois=nobo.pois@AIC, AIC.null=nobo.null@AIC)

#### Goodness-of-fit statistics ####
# Is the best model from AIC comparison also a good model?

library(AICcmodavg)
Nmix.gof.test(nobo.1,nsim=10) #looking for p > 0.05 - null hypothesis
#is that the data fit the model

#### REMOVAL SAMPLING MODEL ####

#using internal data to R

data(ovendata)
hist(ovendata.list$data)

head(ovendata.list$covariates) #ufc is understory foliage cover,
# TRBA is basal area of large trees

ovenFrame <- unmarkedFrameMPois(y=ovendata.list$data,
             siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
             type="removal")

oven.null <- multinomPois(~1 ~1, data=ovenFrame)
oven.ufc <- multinomPois(~1 ~ufc, data=ovenFrame)
oven.trba <- multinomPois(~1 ~trba, data=ovenFrame)
oven.interact <- multinomPois(~1 ~ufc + trba + ufc:trba, data=ovenFrame)

#Rank models by AIC
ms <- fitList(
  "lam(.) p(.)" = oven.null,
  "lam(ufc) p(.)" = oven.ufc,
  "lam(trba) p(.)" = oven.trba,
  "lam(ufc*trba) p(.)" = oven.interact)

(ms1<-modSel(ms))

#### Challenge ####

#The models for the ovenbird data didn't include any covariates for
# detection. How might ufc or trba influence detection probability?
# How do models with those covariates rank relative to the null model
# for detection probability?

#How well does the top model fit the data? What does a goodness-of-fit
# test indicate?


## Answers
oven.p.ufc <- multinomPois(~ufc ~trba, data=ovenFrame)
oven.p.trba <- multinomPois(~trba ~trba, data=ovenFrame)

ms2 <- fitList(
  "lam(.) p(.)" = oven.null,
  "lam(ufc) p(.)" = oven.ufc,
  "lam(trba) p(.)" = oven.trba,
  "lam(ufc*trba) p(.)" = oven.interact,
  "lam(trba) p(ufc)" = oven.p.ufc,
  "lam(trba) p(trba)" = oven.p.trba)

(ms3<-modSel(ms2)) #model with p(trba) ranks higher than p(.)

Nmix.gof.test(oven.p.trba,nsim=10) #data fit the model
