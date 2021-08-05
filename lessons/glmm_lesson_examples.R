#generalized linear model lesson paired R code

library(lme4)
library(pscl)
library(VGAM)
library(MASS)
library(piecewiseSEM)
library(ggplot2)
library(effects)
library(ZIM)
library(glmmTMB)
library(broom.mixed)
library(dotwhisker)
library(multcomp)
library(DHARMa)

#GLMMS

#####River Point Tick Counting Data example##########

dat <- dget('data/rp_tick.dat')

dat$Date <- as.Date(paste(dat$Month, dat$Day, dat$Year, sep = '/'), format = '%m/%d/%Y')

#we know a bit about these data already, so let's jump into how we would add a random effect to this setup

#sampling was uneven across days and we might want to account for that in our analysis

hist(table(dat$Date))

#this can influence our understanding of how important covariates are to tick relative abundance

#so let's look at what a random effect looks like using lme4

fm <- glmer(Number.of.Tick.Collected ~ (1|Date), data = dat, family = 'poisson')
summary(fm)
plot(fm)

#we know that a Poisson model doesn't fit these data very well, but we can tell that the random effect is covering a fair bit of variance

#we can isolate and explore variation in daily rates of tick exposure by pulling the random effect levels

ranef(fm)

#quick aside: what if we tried to do this via a fixed effect

fm <- glm(Number.of.Tick.Collected ~ as.factor(Date), data = dat, family = 'poisson')
summary(fm)
plot(fm)

#just a terrible result
#that is what an overparameterized model looks like, and is one of the reasons why you can't estimate complex categorical variables as a fixed effect


#let's compare the glmm to a similar glm that we looked at earlier

fm1 <- glm(Number.of.Tick.Collected ~ Band.Code + Year, data = dat, family = 'poisson')
fm2 <- glmer(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year, data = dat, family = 'poisson')

summary(fm1)
summary(fm2)

#calculate pseudo R2 for each of the models (which are not perfect estimates of model fit but they can be useful)
1 - (fm1$deviance/fm1$null.deviance)
rsquared(fm2)

#the mixed model shows a fair bit of improvement (though it's a still not great)

#note how the Year effect completely flips once we account for the individual effect (which is working across years)

#this plot can show how dramatic this effect can be
plot(predictorEffects(fm1))
plot(predictorEffects(fm2))
     
#random effects can be really important, particularly when they have the potential to interact with fixed effects


####Zero-inflation models
#we can't use pscl, so we'll be using another package called glmmTMB
#the code looks a bit different
#note that we can fit a model for both the count component and the bernoulli component of the ZIP

fm3 <- glmmTMB(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year,
              zi = ~ 1,
              data =dat,
              family = poisson)

summary(fm3)
fm3_simres <- simulateResiduals(fm3)
plot(fm3_simres)
car::Anova(fm3)
plot(allEffects(fm3))
g1 <- glht(fm3, linfct = mcp(Band.Code = "Tukey"))
summary(g1)

if (requireNamespace("broom.mixed") && requireNamespace("dotwhisker")) {
  t1 <- broom.mixed::tidy(fm3, conf.int = TRUE)
  t1 <- transform(t1,
                  term=sprintf("%s.%s", component, term))
  if (packageVersion("dotwhisker")>"0.4.1") {
    dw <- dwplot(t1)
  } else {
    owls_nb1$coefficients <- TRUE ## hack!
    dw <- dwplot(owls_nb1,by_2sd=FALSE)
  }
  print(dw+geom_vline(xintercept=0,lty=2))
}

#a more complex model where we add covariate to the ZI part of the model

fm4 <- glmmTMB(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year,
               zi = ~ Band.Code + Year,
               data =dat,
               family = poisson)

summary(fm4)
fm4_simres <- simulateResiduals(fm4)
plot(fm4_simres)
car::Anova(fm4)
plot(allEffects(fm4))

#in this case, the additional covariates don't seem that helpful

#note that this is just the glmer model above but in glmmTMB format
fm5 <- glmmTMB(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year,
               zi = ~ 0,
               data =dat,
               family = poisson)

summary(fm5)
fm5_simres <- simulateResiduals(fm5)
plot(fm5_simres)
car::Anova(fm5)
plot(allEffects(fm5))

#NB and ZINB with a random effect

fm6 <- glmmTMB(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year,
               zi = ~ 0,
               data =dat,
               family = nbinom2)

summary(fm6)
fm6_simres <- simulateResiduals(fm6)
plot(fm6_simres)
car::Anova(fm6)
plot(allEffects(fm6))

fm7 <- glmmTMB(Number.of.Tick.Collected ~ (1|Date) + Band.Code + Year,
               zi = ~ 1,
               data =dat,
               family = nbinom2)

summary(fm7)
fm7_simres <- simulateResiduals(fm7)
plot(fm7_simres)
car::Anova(fm7)
plot(allEffects(fm7))

AIC(fm3, fm4, fm5, fm6, fm7)


#in summary, count data are hard to fit and random effects can provide a powerful way to account for complex issues
#it looks like NB is best for these data
#note that we are still getting bad fit from the ZINB. This is a problem with these kinds of models, sometimes it's hard to parse all these zeros

#Challenge
#explore other GLMMs with the tick data set
#try using BandID as a random effect and see what some of your models suggest
#HINT: make sure you are getting good estimates for your parameters

