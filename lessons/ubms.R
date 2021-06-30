#ubms tutorial
install.packages("ubms")
library(ubms)
library(unmarked)

nobo <- read.csv("nobo_abund.csv")
nobo.umf <- unmarkedFramePCount(y=nobo[,2:4], 
                                siteCovs=cbind(as.data.frame(scale(nobo[,5:10])),site=c(seq(1,25),seq(1,26))), 
                                obsCovs=list(sky=scale(nobo[,11:13]),jdate=scale(nobo[,14:16]),
                                             time=scale(nobo[,17:19])))
summary(nobo.umf)

nobo.1 <- pcount(~sky + jdate + time ~BA + Evergreen5km, 
                 data=nobo.umf,K=105) 
summary(nobo.1)

#run with Stan via ubms package
options(mc.cores=3)

nobo.ubms = stan_pcount(~sky + jdate + time ~BA + Evergreen5km, 
                data=nobo.umf, K=105, chains=3, iter=300)

#compare unmarked and Stan output
cbind(unmarked=coef(nobo.1), stan=coef(nobo.ubms))

#Random effect of site
nobo.ubms.rand = stan_pcount(~sky ~BA + (1|site), 
                        data=nobo.umf, chains=3, iter=300)

#another example this time with mallard data
data("mallard")
str(mallard.y)
str(mallard.obs) #ivel is measure of survey effort
str(mallard.site)

umf.mall = unmarkedFramePCount(y=mallard.y, 
                        siteCovs=mallard.site,
                        obsCovs = mallard.obs)

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