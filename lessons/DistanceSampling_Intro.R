#NAOC 2020 count data workshop distance sampling example
#By Evan Adams and Beth Ross
###########################################################

#we want to survey the different ways that count data become biased and this is an example that occurs frequently with bird survey data
#when you survey an area you are often better at surveying the area closest to you rather than further away
#distance sampling was created to try and understand how the chances of detecting a bird vary with distance of the object to the observer
#So we've added an additional source of information to the count data, for every count we have a distance to the observer
#We use the distance data to make another model that we combine with a simple Poisson regression model so that we can describe the patterns of abundance after accounting for detection bias
#once we estimate this source of bias, then we can account for the birds that we didn't see and correct our estimates

#we are using these packages today, unmarked is the package that does the distance sampling, AICcmodavg gives us some additional tools for unmarked, reshape2 and dplyr help format data, and ggplot2 is essentially useful for all R projects
library(unmarked)
library(AICcmodavg)
library(reshape2)
library(dplyr)
library(ggplot2)

####EXAMPLE USING REAL DATA#######
#let's pull in some data on Chipping Sparrow abundance from the Maine Bird Atlas
#these are in from the Github distro (https://github.com/evanmadams/NAOC2020-Count-Data-Workshop) and need to be in a folder called 'data' in your working  directory for this code to work
dat <- dget('data/mebba_spp_Chipping Sparrow.dat')

#these data are in a bit of weird format for JAGS so I'm going to reorganize them into something that unmarked will recognize
#we create a list of sites that were visited
y <- data.frame(site = 1:dat$nsites)

#then we pull all the distance classes of our CHSP observation and the site that they belong to
#here we have 5 distance classes: 20m intervals from 0-100m
ds <- data.frame(site = dat$site.idx, dc = dat$distclass)

#then we reorganize the data so I can count the number of individuals in each distance category in each site
ds <- dcast(ds, site ~ dc, fun.aggregate = length)

#finally we join back with the complete site list so that we know the sites that were surveyed
y <- left_join(y, ds, 'site')

#turn NAs in the merge into true zeroes as NA means that no birds from this species was detected
y[is.na(y)] <- 0

#pull in covariates that describe either abundance or detection probability
#these include habitat data and conditions under which the survey was conducted 
#Note that the wind speed, sky condition and background noise are categorical and the date/time data are scaled continuous

sc <- data.frame(dat$habcovs, dat$distcovs, dat$avcovs)

#then we need to convert it to a format that unmarked understands
#unmarked has a function that helps us do that specifically for distance sampling data (unmarkedFrameDS)
#I'll just be grabbing the first 100 sites so that the models converge quickly but you can explore the full data set if you wish

umf <- unmarkedFrameDS(y = as.matrix(y[1:100 ,2:6]), siteCovs = sc[1:100,], dist.breaks = c(0, 20, 40, 60, 80, 100), survey = 'point', unitsIn = 'm')   #parse the data into the unmarked data frame and subset the data to (1) have models converge faster and (2) remove birds with unknown distances

#let's have a look at these data to get a sense of what kinds of models will be useful to build with them
#first, we should have a look at the count data as that could give us some insight into what kind of error distribution to select

mean(rowSums(umf@y))

#changing the plot settings to compare two figures in the same frame
par(mfrow=c(1,2))

#plotting the histogram of the actual count data
hist(rowSums(umf@y))

#does a Poisson distribution with the same mean look similar?
hist(rpois(nrow(umf@y), mean(rowSums(umf@y))))

#next, let's look at these distance data because these are the foundation for how we estimate detection probability
barplot(colSums(y[,2:6]), ylab = 'Counts', xlab = 'Distance Band')

#well that doesn't look like it's going down with distance, what gives?
#the area of each distance band is also increasing because it's a circle
bandarea <- c( (pi * 20^2), (pi * 40^2) - (pi * 20^2), (pi * 60^2) - (pi * 40^2), (pi * 80^2) - (pi * 60^2), (pi * 100^2) - (pi * 80^2) )

barplot(bandarea, ylab = 'Area', xlab = 'Distance Band')

#so the area covered by the outer bands is significantly larger than the inner bands
#when we correct the counts for this pattern, this is what we see

barplot(colSums(y[,2:6])/bandarea, ylab = 'Count/Area', xlab = 'Distance Band')

#okay, that looks a bit more like what we expected. Now we just have to find a function that fits that decrease
#we are going to try a half-normal distribution
#briefly, this is what a half-normal curve looks like, where y is >= 0

dhalfnorm <- function(y, sigma){
  
  out <- (sqrt(2)/(sigma * sqrt(pi))) * exp(-1 * (y^2/(2 * sigma^2)))
  return(out)
  
}

barplot(dhalfnorm(y = 0:100, sigma = 40))

#so we are saying that we think that detection drop-off with distance looks like this -- the positive half of a Gaussian curve
#this seems like a good start for a detection model for these data, so now let's run a model in unmarked

#we'll use 'distsamp' in unmarked to build a model based on what we have learned from the data
#the first element of the function needs a model formula that is similar to what we were using earlier only we need multiple formulas for the distance model and the ecological model
#we also designate how we want to model the detection probability over distance (a half-normal distribution in this case)

#here we are running a null model without any covariates on either submodel

fm <- distsamp(~1 ~1, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm) #look at the parameter estimates from the model
par(mfrow=c(1,1))
plot(fm) #plot the residuals

#LOOK AT OVERALL MODEL FIT
#one of the best ways to do this is to use bootstrapping to determine how our model fit compares to all the boostrapped fits
fit <- Nmix.gof.test(fm, nsim = 100)

#in our case the observed model (the red line) is not so different from all the bootstrapped cases (hence the P > 0.05) and that overall model fit is likely good

#also can look at something called c-hat to see if there is evidence that the count data are overdispersed to the modeled error estiamtes
print(fit)

#c-hat values close to 1 indicate that overdispersion is fairly minimal that our model fits the data reasonable well


#DISTANCE MODEL FIT

#let's quantify the bias we saw in the counts
#to do this we'll look at the detection portion of the model to see how detectability varied with distance to the observer
backTransform(fm, type = 'det')   #backtransform the model estimates
hist(fm, xlab = 'Distance (m)')   #look at how well the model fits our distance data

#to determine how many birds we miss in a 100m radius point count circle, first we calculate the survey half-width
ehw <- integrate(gxhn, 0, 100, sigma = exp(coef(fm)[2]))$value  #calculate the survey half-width
ehw/100   #100 is the maximum survey distance and used to scale the half-width and calculate the proportion of birds detected


#look at our estimates of density
backTransform(fm, type = 'state')  #average across all sites
ranef(fm, K = 50)   #Empirical Bayes estimates of abundance for each site (and random effect estimated by the model), K should be a value so high that the number of individuals at that site couldn't reach it

plot(ranef(fm, K = 50))

#we can also add covariates to either the abundance or detection models to test hypotheses about the importance of those variables

#what if detection probability was influenced by wind speed at the site?
fm <- distsamp(~Wind.Speed.Code ~1, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm)

#not much evidence for that in this data set, but perhaps it's an effect that we need more data to estimate properly

#what if the presence of developed habitat affected abundance at the site?
fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'halfnorm', output = 'density', unitsOut = 'ha')

summary(fm)

#Well that's pretty strong evidence that CHSP densities are increasing in developed habitat

#let's try a more complicated model that suggests that detections are a function of wind speed, background noise, and time of year; and abundance is forest and developed habitat

fm <- distsamp(~Wind.Speed.Code + Background.Noise.Code + ydate ~eff.domhabA + eff.domhabF, data = umf, keyfun = 'halfnorm', output = 'density', unitsOut = 'ha')

summary(fm)

#the most important variable for abundance appears to be agricultural habitat while date is very important to detection
#note that the parameter estimates for habitat and wind speed changed when compared to the earlier models

#let's check the fit of this model
fit <- Nmix.gof.test(fm, nsim = 100)

print(fit)

#even better than before, so that suggests that the covariates we added were perhaps useful (thought it's hard to tell those things with a test like this, essentially we know that the new model isn't worse)

#now what do we need to visualize these effects?
#to do this, we want to use the predict function in unmarked to show these effects
#so first let's create some new data that we want to predict to

#for abundance we want to predict for different combinations of the dominant habitat types

newdat <- data.frame(eff.domhabF = c(0, 0, 1, 1), eff.domhabA = c(0, 1, 0, 1))

ps <- predict(fm, newdat, type = 'state')

#now we have predicted densities and 95% CIs for each combination of mixed forest and developed habitat
#let's add a couple things to the data frame and then plot the differences

ps <- data.frame(ps, newdat)
ps$Habitat <- c('None', 'Developed', 'Mixed Forest', 'Both')


p <- ggplot(ps, aes(x = Habitat, y = Predicted, ymin = lower, ymax = upper))
p <- p + geom_point() + geom_errorbar(aes(width = 0.2))
p <- p + theme_bw() + ylab('Predicted Density (birds/ha)')

plot(p)

#so now we can see how both developed and mixed forest habitat contributes to the density of CHSP

#we can also visualize how covariates influence detection probability

#effect plots for detection
#This time we will make predictions for each factor while holding the others constant. So we vary one and use the mean for the others

#wind speed

newdat <- data.frame(Wind.Speed.Code = c(0, 1, 2, 3, 4), Background.Noise.Code = mean(sc$Background.Noise.Code), ydate = 0)

ps <- predict(fm, newdat, type = 'det')
ps <- data.frame(ps, newdat)

#We estimated sigma but we still need to convert this to detection probability, so we do that by using a similar process to above
#only I'll use lapply to speed up the process across multiple predictions

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = Wind.Speed.Code, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Wind Speed (Beaufort scale)')

plot(p)

#Background noise

newdat <- data.frame(Wind.Speed.Code = mean(sc$Wind.Speed.Code), Background.Noise.Code = c(0, 1, 2), ydate = 0)

ps <- predict(fm, newdat, type = 'det')

ps <- data.frame(ps, newdat)

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = Background.Noise.Code, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Background Noise (categorical)')

plot(p)

#Date

newdat <- data.frame(Wind.Speed.Code = mean(sc$Wind.Speed.Code), Background.Noise.Code = mean(sc$Background.Noise.Code), ydate = seq(from = min(sc$ydate), to = max(sc$ydate), by = .1))

ps <- predict(fm, newdat, type = 'det')

ps <- data.frame(ps, newdat)

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = ydate, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Date (scaled)')

plot(p)

#Okay, now you've seen the basics of what you can do with distance sampling in unmarked, let's test some of your new skills

####CHALLENGE####

#Often we are interested in finding the detection function that fits the data that we have the best
#Try your hand at doing this. Make a few models to compare different detection functions

#note please use 'aictab' in the AICcmodavg library to compare multiple models

#check out ?distsamp for the other types of detection functions in can fit

#[see below for one way to do this]








##answer

fm1 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'halfnorm', output = 'density')
fm2 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'exp', output = 'density')
fm3 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'hazard', output = 'density')
fm4 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'uniform', output = 'density')


aictab(list(fm1, fm2, fm3, fm4), second.ord = FALSE)  #create an AIC table to compare the results

#compare the fits of the various models

par(mfrow=c(2,2))
hist(fm1, xlab = 'Distance (m)')
hist(fm2, xlab = 'Distance (m)')
hist(fm3, xlab = 'Distance (m)')
hist(fm4, xlab = 'Distance (m)')

#quick note: notice that the exponential detection function didn't work that well (i.e., produced NANs)
#exp models are often poor performers and you should probably try to avoid them
#hn and hazard models are much more useful

####SIMULATION EXAMPLE#######

#let's simulate some data to show you how to construct a scenario where we can use distance sampling

#simulate the total N at 100 sites using a Poisson

N <- rpois(100, 3)

par(mfrow = c(1,1))
hist(N)

#let's assume that we have a 30% detection probability within the survey area assuming a half-normal detection function
#so the observations (ys) would look something like this

Y <- rbinom(100, N, 0.5)

hist(Y)

Y <- data.frame(site = 1:length(Y), Y)

nobs <- sum(Y$Y)

obs.idx <- rep(1:nrow(Y), times = Y$Y)

#now that we have our observations we can assign them a distance of first detection
#first we need to determine what 30% detection probability would be in our survey

#assuming a 500m transect, a 30% detection probability would be achieved with a sigma 

ehw <- integrate(gxhn, 0, 500, sigma = 200)$value  #calculate the survey half-width
ehw/500   

rhalfnorm <- function(n = 1, mean = 0, sigma = 1){
  
  abs(rnorm(n, mean, sigma))
  
}

dhalfnorm <- function(y, sigma){
  
  out <- (sqrt(2)/(sigma * sqrt(pi))) * exp(-1 * (y^2/(2 * sigma^2)))
  return(out)
  
}


d <- rhalfnorm(nobs, sigma = 200)

hist(d)

breaks <- c(0, 100, 200, 300, 400, 500) #create distance breaks to divide the observations

d[d > max(breaks)] <- max(breaks) # truncate the data so that observations beyond the final break are included in the final category

dclass <- cut(d, breaks = breaks)
dclass <- data.frame(site = obs.idx, dc = dclass)
dclass <- dcast(dclass, site ~ dc, fun.aggregate = length)
Y <- left_join(Y, dclass, 'site')
Y[is.na(Y)] <- 0

#create our unmarked data frames as we did before
umf <- unmarkedFrameDS(y = as.matrix(Y[1:100 ,3:ncol(Y)]), dist.breaks = breaks, survey = 'line', tlength = rep(1000, 100), unitsIn = 'm')
summary(umf)
fm <- distsamp(~1 ~1, data = umf, keyfun = 'halfnorm', output = 'abund')
summary(fm)
plot(fm)

#let's check the outputs to see how similar they are to our known parameter values
backTransform(fm, type='state')
#that's pretty close, but let's look at the confidence intervals to see if our initial values are inside it
exp(confint(fm, type = 'state'))
#they are, so this model seems fairly reasonable

#let's do the same thing for detection probability
backTransform(fm, type="det") 

# Effective strip half-width
eshw <- integrate(gxhn, 0, 500, sigma=exp(coef(fm)[2]))$value

# Detection probability
eshw / 500 # 500 is strip-width

#yep, that is solid also

#and overall model fit also looks fine
fit <- Nmix.gof.test(fm, nsim = 100)
print(fit)

#note that you if want to simulate more complicated distance sampling data, check out library(DSsim), it's very cool

###END LESSON###
