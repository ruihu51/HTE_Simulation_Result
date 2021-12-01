rm(list=ls())

library(devtools)
devtools::install_github("ruihu51/CausalSim")
require(CausalSim)

# 1) generate data function
generate.data <- function(n){
  beta <- c(0.25, 0.25)
  pi0 <- function(w) 0.5+(w[,1]/3)
  mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

  psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
  theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))

  simulated.data <- list(Y=Y, A=A, W=W, psi0=psi0, theta0=theta0)
  return(simulated.data)
}
set.seed(2298)
n <-1000
beta <- c(0.25, 0.25)
pi0 <- function(w) 0.5+(w[,1]/3)
mu0 <- function(a, w, beta) 0.1 + beta[1]*a + beta[2]*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
A <- rbinom(n, size=1, prob=pi0(W))
Y <- rnorm(n, mean=mu0(A, W, beta), sd=1)

# psi0 <- mean((mu0(1, W, beta) - mu0(0, W, beta))^2)
psi0 <- (beta[1]+beta[2]*(5/6))^2+(beta[2]^2*(1/5-1/9))+beta[2]^2*1/4
# theta0 <- var((mu0(1, W, beta) - mu0(0, W, beta)))
theta0 <- (beta[2]^2*(1/5-1/9))+beta[2]^2*1/4

# update control
control = list(est.type = list('psi.est'='all', 'theta.est'='all'),
               conf.int = TRUE,
               conf.int.type = 'Wald')
control <- hte.measure.NullTest.control(control)
control

# estimated propensity
prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                         newX = data.frame(W),
                         SL.library = control$pi.SL.library,
                         family = binomial(),
                         obsWeights=rep(1,n),
                         id=1:n)
control$pi.hat <- prop.reg$SL.predict
plot(pi0(W), control$pi.hat)
abline(0, 1, col='red')
prop.reg$coef
prop.reg$cvRisk

# estimated outcome regression
AW <- cbind(A, data.frame(W))
if(length(setdiff(Y, c(0,1))) == 0) {
  family = 'binomial'
} else {
  family = 'gaussian'
}
# control$mu.SL.library <- c("SL.mean", "SL.glm.interaction", "SL.gam", "SL.earth")
mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                       newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                       SL.library = control$mu.SL.library,
                       family = family,
                       obsWeights=rep(1,n),
                       id=1:n)
control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0

mu1.0 <- mu0(a = 1, w = W, beta = beta)
mu0.0 <- mu0(a = 0, w = W, beta = beta)
tau0 <- mu1.0 - mu0.0

plot(mu1.0, control$mu.hats$mu1)
abline(0, 1, col='red')
plot(mu0.0, control$mu.hats$mu0)
abline(0, 1, col='red')

# estimated tau
tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))

plot(tau0, tau.hat)

# 1) Why is the tau.hat a constant?
#
mu1.library <- mu.reg$library.predict[1:n,]
mu0.library <- mu.reg$library.predict[-(1:n),]
mu.SL <- data.frame(mu1.hat=mu1.library[,"SL.gam_All"], mu0.hat=mu0.library[,"SL.gam_All"])
# mu.SL <- data.frame(mu1.hat=mu1.library[,"SL.glm.interaction_All"], mu0.hat=mu0.library[,"SL.glm.interaction_All"])
mu.SL %>%
  mutate(tau.hat = mu1.hat-mu0.hat) %>%
  head()
mu.reg$cvRisk
mu.reg$coef

# where did tau.hat come from?
X = data.frame(cbind(A, W))
newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W)))
cts.num <- 4
deg.gam = 2
cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
if (sum(!cts.x) > 0) {
  gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
} else {
  gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
}
fit.gam <- gam::gam(gam.model, data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50), weights = rep(1,n))
fit.gam$coefficients[4]

j <- 3
fit.gam$coefficients[1] + fit.gam$coefficients[4]*newX[j,1] + fit.gam$coefficients[5]*newX[j,4] +
  fit.gam$coefficients[2]*newX[j,][2] + fit.gam$coefficients[3]*newX[j,][3] +
  fit.gam$smooth[j,1] + fit.gam$smooth[j,2]

# let's try some different SL models
# c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
learners = create.Learner("SL.earth", params = list(penalty=-1)) # penalty=-1 keep all variables

# control$mu.SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
control$mu.SL.library <- c("SL.earth")
control$mu.SL.library <- c(learners$names)
# control$mu.SL.library <- c("SL.mean", "SL.glm.interaction", "SL.gam", learners$names)
mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                       newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                       SL.library = control$mu.SL.library,
                       family = family,
                       obsWeights=rep(1,n),
                       id=1:n)
control$mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
control$mu.hat <- A * control$mu.hats$mu1 + (1-A) * control$mu.hats$mu0

mu1.0 <- mu0(a = 1, w = W, beta = beta)
mu0.0 <- mu0(a = 0, w = W, beta = beta)
tau0 <- mu1.0 - mu0.0

plot(mu1.0, control$mu.hats$mu1)
abline(0, 1, col='red')
plot(mu0.0, control$mu.hats$mu0)
abline(0, 1, col='red')

# estimated tau
tau.hat <- control$mu.hats$mu1 - control$mu.hats$mu0
Z.hat <- (2*A - 1) / (A * control$pi.hat + (1-A) * (1-control$pi.hat))

plot(tau0, tau.hat, col=W[,3]+1)

mu.reg$fitLibrary
mu.reg$fitLibrary$SL.earth_All$object$bx %>%
  head()

mu.reg$fitLibrary$SL.earth_1_All$object$bx %>%
  head()
