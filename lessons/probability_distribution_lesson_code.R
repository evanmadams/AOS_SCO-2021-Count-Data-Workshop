#NAOC 2020 count data workshop distance probability distribution lesson
#By Evan Adams and Beth Ross
###########################################################

#Not a lot of code here, just showing how we created the distributions shown in the presentation

#coin flipping

hist(rbinom(10, 1, 0.5))


#probability density of the binomial distribution

plot(dbinom(0:10, 10, 0.5))


#distribution types

#gaussian
hist(rnorm(100, 0, 1))

#lognormal
hist(rlnorm(100, 0, 1))

#gamma
hist(rgamma(100, 1, 1))

#binomial
hist(rbinom(100, 10, 0.5))

#poisson
hist(rpois(100, 3))


#probability mixtures

#showing two unmixed gaussians
hist(c(rnorm(100, -25, 3), rnorm(100, 25, 7)), freq = FALSE, breaks = -50:50)

#then combining two gaussians
hist(c(rnorm(100, -25, 3) + rnorm(100, 25, 7)), freq = FALSE, breaks = -50:50)

#now let's create a zero-inflated poisson

hist(rpois(100, 4))

hist(rbinom(100, 1, 0.1))

rzipois <- function(n,p,lambda) {
  z <- rbinom(n,size=1,prob=p)
  y <- (1-z)*rpois(n,lambda)
  return(y)
}


hist(rzipois(100, 0.9, 4))

#negative binomial

hist(rnbinom(100, mu = 4, size = 1.4))


