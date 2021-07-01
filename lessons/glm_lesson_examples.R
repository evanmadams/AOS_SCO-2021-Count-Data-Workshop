#generalized linear model lesson paired R code

library(lme4)
library(pscl)
library(VGAM)
library(MASS)
library(piecewiseSEM)
library(ggplot2)
library(effects)
library(ZIM)

#What is a GLM

#generate a simple regression data set

mu <- 5
sd <- 1
b1 <- .5

x <- rnorm(200, 0, 1)

y <- mu + b1*x + rnorm(length(x), 0, sd)

plot(x, y)
abline(lm(y ~ x))

#example of normally distributed data

hist(rnorm(1000, 0, 1))

#generate some poisson data to see what that looks like

d <- data.frame(y = rpois(100, 3))

ggplot(d, aes(x = y)) + geom_bar() + theme_bw()

#generate a PDF to match the random draw we took from the data

pdf <- data.frame(y = dpois(0:10, 3), x = 0:10)

ggplot(pdf, aes(x = x, y = y)) + geom_point() + theme_bw()



####Poisson simulated GLM example#########

mu <- 0.5
b1 <- 0.5

x1 <- rnorm(200, 0, 1)

y1 <- rpois(length(x), exp(mu + b1*x1))

hist(y)

plot(x, y)
abline(lm(y ~ x1))



####ZIP simulated GLM example##########

mu <- 0.5
psi <- 0.5
b1 <- 0.5

x2 <- rnorm(200, 0, 1)

y2 <- rzipois(length(x), exp(mu + b1*x2), psi)

hist(y2)

plot(x, y)
abline(lm(y ~ x2))


#compare the two simulated data sets
par(mfrow = c(1, 2))
hist(y1)
hist(y2)

dev.off()

sim <- data.frame(pois = y1, zip = y2, x1 = x1, x2 = x2)

#fit the appropriate (and inappropriate models to them)

fm.p <- glm(pois ~ x1, data = sim, family = 'poisson')
plot(fm.p)
pR2(fm.p)
summary(fm.p)

fm.zip <- zeroinfl(zip ~ x2|1, data = sim, dist = 'poisson')
pR2(fm.zip)
summary(fm.zip)

#note that you can use your poisson simulation with the zip model

fm.zip2 <- zeroinfl(pois ~ x1|1, data = sim, dist = 'poisson')
pR2(fm.zip2)
summary(fm.zip2)

#it just doesn't find any evidence of zero-inflation

#but if you do the opposite, things don't go as well

fm.p2 <- glm(zip ~ x2, data = sim, family = 'poisson')
plot(fm.p2)
pR2(fm.p2)
summary(fm.p2)

#you get biased estimates because the model can't disentangle the ZI parameter

#the model thinks your data look like this:

hist(rpois(200, exp(-0.21)))

#when they actually look like this:

hist(y2)

#that's pretty far off and these differences can lead to significant errors in inference



#####River Point Tick Counting Data example##########

dat <- dget('data/rp_tick.dat')

#let's look at our Poisson response variable, number of ticks collected

hist(dat$Number.of.Tick.Collected)

#this fits the general characteristics of a Poisson variable but wow that's a lot of zeroes
#let's test a hypothesis that tick abundance on birds is variable across months and that it depends if we've caught them before
#(we remove ticks when we sample them so the ones we catch again probably have fewer ticks to find)
#let's see what happens when we fit a Poisson GLM to these data

#fit a poisson glm

summary(dat)

fm <- glm(Number.of.Tick.Collected ~ Band.Code + Month + Year, data = dat, family = 'poisson')

summary(fm)

anova(fm, test = 'LRT')

plot(fm)

#a quick pseudo-r2 calculation for the model
pR2(fm)

#let's plot the effects so we can see how important each of the variables are

plot(predictorEffects(fm))

#so we had a small difference between years and more ticks on new birds than recaps. Cool.


#so while we saw significant results with this model we also saw fairly poor model fit, particularly underestimating high values
#this is one of the weakness of a Poisson model where the distribution mean and variance are the same parameter, it's not as flexible as other distributions
#we can do a quick comparison between the tick data and randomly generated Poisson data with the same mean to see what we mean

ggplot(dat, aes(x = Number.of.Tick.Collected)) + geom_bar() + theme_bw()

rdat <- data.frame(y = rpois(nrow(dat), mean(dat$Number.of.Tick.Collected)))

ggplot(rdat, aes(x = y)) + geom_bar() + theme_bw()

#pretty big difference right?
#the real data has a lot more zeroes and a lot more high values than the randomly generated Poisson data
#this happens pretty often in real data sets and is called a lot of things but mostly 'zero-inflation' or 'overdispersion'
#we have some kinds of distributions that can help us deal with this issue so let's try some out and see if we can improve model fit


#so let's try to fit a zero-inflated Poisson model (ZIP) and see if that does any better

fm.zip <- zeroinfl(Number.of.Tick.Collected ~ Band.Code + Month + Year|1, data = dat, dist = 'poisson')

summary(fm.zip)

pR2(fm.zip)

#note that we now have two different models, one that describes the poisson data and the other that describes the zero-inflation data
#the model seems to estimate the zero-finlation parameter well and it's significantly different than zero, but does that mean it mattered?

#we can plot the fitted/residuals plot to see if there is any evidence of correlation and thus lack of fit
plot(residuals(fm.zip) ~ fitted(fm.zip))

#let's compare the two models and see which describes the data best
#the vuong test isn't perfect for this job but it can be helpful and suggests that the zip glm is much better than the regular poisson glm
vuong(fm, fm.zip)

#we can play a similar game to before and look at what randomly generated data from zip would look like in comparison to our data

#use the glm to pull overall averages for the Poisson mean and zero-inflation probability

zeroinfl(Number.of.Tick.Collected ~ 1 | 1, data = dat, dist = 'poisson')

#plot the raw data
ggplot(dat, aes(x = Number.of.Tick.Collected)) + geom_bar() + theme_bw()

#now generate zero-inflated data from the model averages and see how similar that is to the original data set

rdat <- data.frame(y = rzipois(nrow(dat), exp(1.029), plogis(2.316)))

ggplot(rdat, aes(x = y)) + geom_bar() + theme_bw()

#it's starting to fit the data a bit better, the number of zeroes have increased and we are seeing higher maximum counts
#we have a maximum count of 82 and that is really hard to account for when your mean count is 0.25
#it's a lot of work to deal with overdispersed count data, particularly if you are interested in modeling the extremes

#so how to do it?

#last example: ZINB

fm.zinb <- zeroinfl(Number.of.Tick.Collected ~ Band.Code + Month + Year | 1, data = dat, dist = 'negbin')
summary(fm.zinb)
plot(residuals(fm.zinb) ~ fitted(fm.zinb))
pR2(fm.zinb)
vuong(fm, fm.zinb)
vuong(fm.zip, fm.zinb)
zeroinfl(Number.of.Tick.Collected ~ 1 | 1, data = dat, dist = 'negbin')
ggplot(dat, aes(x = Number.of.Tick.Collected)) + geom_bar() + theme_bw()
rdat <- data.frame(y = rzinb(nrow(dat), lambda = exp(-1.381), k = 0.0474, omega = plogis(-13.72)))
ggplot(rdat, aes(x = y)) + geom_bar() + theme_bw()

#getting close, but these data are difficult to analyze!
#note that our fitted vs. residuals plots are improving a bit, perhaps we can start thinking about believing our results