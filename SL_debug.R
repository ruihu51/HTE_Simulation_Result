rm(list=ls())

library(devtools)
devtools::install_github("ruihu51/CausalSim")
require(CausalSim)

# 900712
# load("Data/Estimator/ests.sim.2.25.75.w.sl3.c.RData")
# ret <- subset(ests.sim.2.25.75.w.sl3.c,
#               (type %in% c('psi.est', 'theta.est')))[,c("type", "est", "se", "ll", "ul", "n", "j" , "seed")]
# optional.ret <- subset(ests.sim.2.25.75.w.sl3.c,
#                        (type %in% c('psi.est', 'theta.est', 'psi.plug.in.est', 'psi.one.step.est',
#                                     'theta.plug.in.est', 'theta.one.step.est')))[,c("type", "est", "n", "j" , "seed")]
# pars.ret <- subset(ests.sim.2.25.75.w.sl3.c,
#                    (type %in% 'pars'))[,c("n", "j" , "seed", "psi0", "theta0")]
#
# ests.sim.2.25.75.w.sl.ret <- join(ret, pars.ret, by=c("n", "j", "seed"))
# ests.sim.2.25.75.w.sl.optional.ret <- join(optional.ret, pars.ret, by=c("n", "j", "seed"))
#
# ests.sim.2.25.75.w.sl.ret %>%
#   filter((type %in% c("theta.est")) & (n==2000)) %>%
#   filter(se <= 1.84) %>%
#   head()

######################
# run one single case
######################
# 1) generate data
set.seed(128873)
n <-2000
beta <- c(0.25, 0.75)
pi0 <- function(w) expit(0.5+(w[,1]/3))
mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
A <- rbinom(n, size=1, prob=pi0(W))
Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

psi0 <- (beta[1]+beta[2]*(5/6))^2+(beta[2]^2*(1/5-1/9))+beta[2]^2*1/4
theta0 <- (beta[2]^2*(1/5-1/9))+beta[2]^2*1/4

mu1.0 <- mu0(a = 1, w = W, beta = beta)
mu0.0 <- mu0(a = 0, w = W, beta = beta)
tau0 <- mu1.0 - mu0.0

# 2) update control for outcome regression
# est mu.hat together
# include AW interaction term
SL.gam.interaction <- function(Y, X, newX, family, obsWeights, id, deg.gam = 2, cts.num = 4, ...) {
  # load required packages
  # require('pkg')
  # if(family$family == 'gaussian') {
  #
  # }
  # if(family$family == 'binomial') {
  #
  # }
  if(!require('gam')) {stop("SL.gam requires the gam package, but it isn't available")}
  if("mgcv" %in% loadedNamespaces()) warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  # create the formula for gam with a spline for each continuous variable
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  gam.model <- as.formula("Y ~ s(W1, 2) + s(W2, 2) + s(W1):A + s(W2):A + A*W3")
  fit.gam <- gam::gam(gam.model, data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50), weights = obsWeights)
  if(packageVersion('gam') >= 1.15) {
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  # pred is the predicted responses for newX (on the scale of the outcome)
  # pred <- numeric()
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = fit.gam)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.gam.interaction'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}
learners = create.Learner("SL.earth", params = list(penalty=-1))

control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald',
               mu.SL.library = c("SL.gam.interaction", "SL.glm.interaction", learners$names))
control <- hte.measure.NullTest.control(control)
control

# 3) estimated propensity
prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                         newX = data.frame(W),
                         SL.library = control$pi.SL.library,
                         family = binomial(),
                         obsWeights=rep(1,n),
                         id=1:n)
control$pi.hat <- prop.reg$SL.predict
# plot(pi0(W), control$pi.hat)
# abline(0, 1, col='red')
# prop.reg$coef
# prop.reg$cvRisk

# 4) estimated outcome regression
AW <- cbind(A, data.frame(W))
if(length(setdiff(Y, c(0,1))) == 0) {
  family = 'binomial'
} else {
  family = 'gaussian'
}

# a) estimate mu1.hat and mu0.hat in one model
mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                       newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                       SL.library = control$mu.SL.library,
                       family = family,
                       obsWeights=rep(1,n),
                       id=1:n)
control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0
tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0

# plot(mu1.0, control$mu.hats$mu1)
# abline(0, 1, col='red')
# plot(mu0.0, control$mu.hats$mu0)
# abline(0, 1, col='red')

# b) estimate mu1.hat and mu0.hat separately
mu0.reg <- SuperLearner(Y=Y[which(A==0)], X = data.frame(W[which(A==0),]),
                        newX = data.frame(W),
                        SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                        family = family,
                        obsWeights=rep(1,length(which(A==0))),
                        id=1:length(which(A==0)))
mu0.hat2 <- mu0.reg$SL.predict
mu1.reg <- SuperLearner(Y=Y[which(A==1)], X = data.frame(W[which(A==1),]),
                        newX = data.frame(W),
                        SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                        # SL.library = c("SL.gam"),
                        family = family,
                        obsWeights=rep(1,length(which(A==1))),
                        id=1:length(which(A==1)))
mu1.hat2 <- mu1.reg$SL.predict
mu.hat2 <- A * mu1.hat2 + (1-A) * mu0.hat2
tau.hat2 <- mu1.hat2 - mu0.hat2

plot(tau0, tau.hat, ylim=c(-0.5, 2.1))
points(tau0, tau.hat2, col='blue')
abline(0,1, col='red')

# estimated Z.hat
Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))
Z0 <- (2*A - 1) / (A * pi0(W) + (1-A) * (1-pi0(W)))

###################
# estimated theta
##################
gamma.hat <- mean(tau.hat)
gamma.hat2 <- mean(tau.hat2)
gamma0 <- mean(tau0)

##############
# se(eif.hat)
#############
# a) estimate mu1.hat and mu0.hat in one model
theta.plug.in <- mean((tau.hat-gamma.hat)^2) # why not var(tau.hat)?

theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2
sd(theta.eif.hat)
var(theta.eif.hat)
var(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
var(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat) + var((tau.hat-gamma.hat)^2) + 2*var(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat, (tau.hat-gamma.hat)^2)

var(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat)
var((tau.hat-gamma.hat) * (Y - control$mu.hat))
var(tau.hat*(Y - control$mu.hat) - gamma.hat*(Y - control$mu.hat))
var(tau.hat*(Y - control$mu.hat)) + var(gamma.hat*(Y - control$mu.hat)) -
  2*cov(tau.hat*(Y - control$mu.hat), gamma.hat*(Y - control$mu.hat))
var(tau.hat*Y)
var(gamma.hat*Y)
var(tau.hat*control$mu.hat)
var(gamma.hat*control$mu.hat)

var(Y*control$mu.hats$mu1)
var(Y*control$mu.hats$mu0)
var(A*control$mu.hats$mu1^2)
var((1-A)*control$mu.hats$mu0^2)
var(control$mu.hats$mu1*control$mu.hats$mu0*(2*A-1))

# b) estimate mu1.hat and mu0.hat separately
theta.eif.hat2 <- 2 * (tau.hat2-mean(tau.hat2)) * (Y - mu.hat2) * Z.hat + (tau.hat2-mean(tau.hat2))^2
sd(theta.eif.hat2)
var(theta.eif.hat2)
var(2 * (tau.hat2-mean(tau.hat2)) * (Y - mu.hat2) * Z.hat + (tau.hat2-mean(tau.hat2))^2)
var(2 * (tau.hat2-mean(tau.hat2)) * (Y - mu.hat2) * Z.hat) + var((tau.hat2-mean(tau.hat2))^2) +
  2*cov(2 * (tau.hat2-mean(tau.hat2)) * (Y - mu.hat2) * Z.hat, (tau.hat2-mean(tau.hat2))^2)

var(tau.hat2* (Y - mu.hat2)-mean(tau.hat2) * (Y - mu.hat2))
var(tau.hat2* (Y - mu.hat2)) + var(mean(tau.hat2) * (Y - mu.hat2)) -
  2*cov(tau.hat2* (Y - mu.hat2), mean(tau.hat2) * (Y - mu.hat2))

var(Y*mu1.hat2)

# c)
theta.eif <- 2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0 + (tau0-gamma0)^2
sd(theta.eif)
var(theta.eif)

var((tau0-gamma0) * (Y - mu0(A, W, beta)))

var(tau0*Y)
var(gamma0*Y)
var(tau0*mu0(A, W, beta))
var(gamma0*mu0(A, W, beta))

var(Y*mu1.0)
var(Y*mu0.0)
var(A*mu1.0^2)
var((1-A)*mu0.0^2)
var(mu1.0*mu0.0*(2*A-1))


var(2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0 + (tau0-gamma0)^2)
var(2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0) + var((tau0-gamma0)^2) + 2*cov(2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0, (tau0-gamma0)^2)
var(2 * (tau0) * (Y - mu0(A, W, beta)) * Z0)

###############
# var(tau.hat)
##############
# a)
var(tau.hat)
var(control$mu.hats$mu1) + var(control$mu.hats$mu0) - 2*cov(control$mu.hats$mu1, control$mu.hats$mu0)

# b)
var(tau.hat2)
var(mu1.hat2) + var(mu0.hat2) - 2*cov(mu1.hat2, mu0.hat2)

# c)
var(tau0)
var(mu1.0) + var(mu0.0) - 2* cov(mu1.0, mu0.0)

##############
# os.est bias
#############
# a)
theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - control$mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
# b)
mean(2 * (tau.hat2-gamma.hat2) * (Y - mu.hat2) * Z.hat + (tau.hat2-gamma.hat2)^2)
# c)
mean(2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0 + (tau0-gamma0)^2)
sd(2 * (tau0-gamma0) * (Y - mu0(A, W, beta)) * Z0 + (tau0-gamma0)^2)
