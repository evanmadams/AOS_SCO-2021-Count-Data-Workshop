
# We'll start by simulating the data for a basic N-mixture model. In our simulation, abundance
# will be affected by a variable called "vegHt" (vegetation height) and detection will be affected
# by a variabled called "wind" (wind speed).

#### Mallard Data Intro ####

library(unmarked)
data(mallard)

#look at the data
str(mallard.y)
head(mallard.y)

str(mallard.obs) #ivel is a measure of effort
str(mallard.site) #length is transect length

#create a special data structure required for unmarked
umf.mall = unmarkedFramePCount(y=mallard.y, 
                             siteCovs=mallard.site,
                             obsCovs = mallard.obs)

head(umf.mall)

#can now run some basic models to look at detection and 
# abundance in relation to covariates

# When writing formula, detection covariates follow 
# first tilde, then come abundance covariates

#Note that we're using the pcount function for an N-mixture model
mall.nmix1 <- pcount(~ivel ~elev, data=umf.mall) 

#K is the maximum number of individuals you would expect
# to see at a site. You typically don't need to set this
# parameter, but it doesn't hurt. Larger values of K 
# will mean longer computational time.
#mall.nmix1 <- pcount(~ivel ~elev, data=umf.mall, K=104) 

summary(mall.nmix1)
# things to look for: big or small SE values (+10 or ~0),
# sites removed, values of betas

mall.full <- pcount(~ivel + date ~elev+length+forest,
                    data=umf.mall, K=104)

summary(mall.full)

#The null model is the model with only intercepts, no
# covariates
mall.null <- pcount(~1 ~1, data=umf.mall, K=104)
summary(mall.null)

#### Simulating N-Mixture Data ####  

# We start by simulating the number of sites we need
nSites <- 100
set.seed(443)  # we use this function so that we all get the same values for vegHt when simulating random variables

# Create a covariate called vegHt
vegHt <- rnorm(nSites, 10, 3) # Normal distribution with mean 10 and variation 3

# We now standardize vegHt for use in our model (set to mean 0 with sd 1)
vegHt <- scale(vegHt)

# We now need to create our linear model for abundance with vegHt
# The relationship is described by an intercept of -1 and
#    a slope parameter of 2 on the log scale
lambda <- exp(-1 + 2*vegHt)

# Now we simulate abundace at each site. The 
N <- rpois(nSites, lambda)

plot(N)
plot(vegHt,N)

# We can fit a model without detection probability that relates abundance to vegHt 
# using the glm() function with "family=Poisson":

N.glm <- summary(fm.glm1 <- glm(N ~ vegHt, family=poisson))
N.glm

# Do some analysis of the results
plot(vegHt, N, xlab="Vegetation height", ylab="Abundance (N)")

glm1.est <- coef(fm.glm1) #returns intercept and vegHt beta coefficients. should be around -1 and 2.

lines(sort(vegHt), exp(-1 +2*sort(vegHt)), add=TRUE, lwd=3)
lines(sort(vegHt), exp(glm1.est[1] + glm1.est[2]*sort(vegHt)), add=TRUE,
     lwd=3, col="blue")
legend(1, 20, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=3)


#### ADDING DETECTION PROBABILITY ####

#Now we incorporate the detection/observation portion
nVisits <- 3    # number of repeat samples at each site

wind <- rnorm(300,15,5)
wind <- scale(wind)

#now use inverse logit to calculate linear model for detection probability
p <- 1/(1+exp(-(0.5 + 1*wind)))
plot(p,wind)

#arrange p and wind into dataframes
p <- data.frame(J1=p[1:100],J2=p[101:200],J3=p[201:300])
wind.df <- data.frame(wind1=wind[1:100],wind2=wind[101:200],wind3=wind[201:300])

#create matrix to hold each observation at each site for each visit
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
  for(j in 1:nVisits){
  y[i,j] <- rbinom(1, N[i], p[i,j])
}}

# Look at the data
cbind(N=N, y1=y[,1], y2=y[,2], y3=y[,3] ,p1=p[,1], p2=p[,2],p3=p[,3])

#Create data structure unique to unmarked package. Inputs data (y) and site covariates (siteCovs)
umf <- unmarkedFramePCount(y=y, siteCovs=as.data.frame(vegHt), 
                           obsCovs=list(wind=wind.df))
summary(umf)

# Fit a model and extract estimates
# When writing formula, detection covariates follow first tilde, then come abundance covariates

#Note that we're using the pcount function for an N-mixture model
fm.nmix1 <- pcount(~wind ~vegHt, data=umf) 
summary(fm.nmix1)

# Note, estimates of detection coefficients are on the logit-scale
# When covariates are in the model we can use the following to
# backtransform them
beta1 <- coef(fm.nmix1) #coef() extracts coefficients from the nmix model
beta1

#We can now plot the predicted response between the covariates and 
#abundance
veg.plot=data.frame(vegHt=seq(min(vegHt),max(vegHt),length=100))
veg.pred <- predict(fm.nmix1,newdata=veg.plot,type="state")
veg.plot<- data.frame(veg.plot,veg.pred)

library(ggplot2)

ggplot(veg.plot, aes(x=vegHt,y=Predicted,ymin=lower,ymax=upper)) + 
  geom_line() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=15),axis.title=element_text(size=20)) +
  xlab("Vegetation Height") + ylab("Expected Abundance")

#or covariates and detection probability
wind.plot=data.frame(wind=seq(min(wind),max(wind),length=300))
wind.pred<-predict(fm.nmix1,newdata=wind.plot,type="det")
wind.plot<-data.frame(wind.plot,wind.pred)

ggplot(wind.plot, aes(x=wind,y=Predicted,ymin=lower,ymax=upper)) + 
  geom_line(aes(y=Predicted)) + geom_line(aes(y=lower)) + 
  geom_line(aes(y=upper)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=15),axis.title=element_text(size=20)) +
  xlab("Wind Speed (km/hr)") + ylab("Detection Probability")


# Or suppose you want predictions for new values of vegHt, say 1.2 and 3.1
newdat <- data.frame(vegHt=c(3.1,5))
predict(fm.nmix1, type="state", newdata=newdat)

# We can calculate posterior distributions for latent abundance at
# each site as well.
site.N <- ranef(fm.nmix1)
# and get confidence intervals for each of those estimates
confint(site.N,level=0.95)

# Or compare the estimate the overall total abundance to truth
sum(N) #original simulated abundance
sum(bup(site.N)) #total estimated abundance

# Analyze with real data
nobo <- read.csv("nobo_abund.csv")
head(nobo)

# Create unmarked dataframe
#Create data structure unique to unmarked package. Inputs data (y) and site covariates (siteCovs)
nobo.umf <- unmarkedFramePCount(y=nobo[,2:4], 
            siteCovs=as.data.frame(scale(nobo[,5:10])), 
            obsCovs=list(sky=scale(nobo[,11:13]),jdate=scale(nobo[,14:16]),
                         time=scale(nobo[,17:19])))
summary(nobo.umf)

# Fit a model to evaluate covariate affects on NOBO abundance
nobo.1 <- pcount(~sky + jdate + time ~BA + Evergreen5km, 
                 data=nobo.umf,K=105) 
summary(nobo.1)

#We can now plot the predicted response between the covariates and 
#abundance
BA.plot=data.frame(BA=seq(min(scale(nobo[,6])),max(scale(nobo[,6]))),
                   Evergreen5km=0)

nobo.pred=predict(nobo.1,newdata=BA.plot,type="state")

BA.plot <- data.frame(nobo.pred,BA.plot)

nobo.plot <- ggplot(data=BA.plot,aes(x=BA,y=Predicted,ymin=lower,ymax=upper))
nobo.plot <- nobo.plot + geom_line(aes(x=BA,y=Predicted))
nobo.plot <- nobo.plot + geom_line(aes(x=BA,y=lower),color=2) +
  geom_line(aes(x=BA,y=upper),color=2)
nobo.plot <- nobo.plot + theme_bw() + ylab("Predicted Abundance")
plot(nobo.plot)

#### Challenge ####

#How many bobwhites were there overall?

#How many bobwhites would you expect to see at a site with a BA of -3
# and evergreen5km of 0? How does the error of this estimate change?

#How would we create a plot to look at the effect of wind 
# on detection probability?


## Answers

#Total bobwhite
sum(bup(ranef(nobo.1)))

#Prediction at particular point
newdat <- data.frame(BA=-3,Evergreen5km=0)
predict(nobo.1, type="state", newdata=newdat) #Note the high SE!

#plot for detection probability and wind
sky.plot=data.frame(sky=seq(min(scale(nobo[,11:13])),max(scale(nobo[,11:13])),length=100),
                    jdate=0,time=0)

p.pred <- predict(nobo.1,newdata=sky.plot,type="det")

sky.plot<-data.frame(sky.plot,p.pred)

p.plot <- ggplot(data=sky.plot,aes(x=sky,y=Predicted,ymin=lower,ymax=upper)) 
p.plot <- p.plot + geom_line(aes(x=sky,y=Predicted)) + 
  geom_line(aes(y=lower),color=2) + geom_line(aes(y=upper),color=2)
p.plot <- p.plot + theme_bw() + ylab("Detection Probability")
plot(p.plot)
