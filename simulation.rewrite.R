rm(list=ls())

# packages
library(SuperLearner)

# utils
expit <- function(x) 1 / (1 + exp(-x))

n <- 2000
#################
# generate data
#################
# W: continuous
# A: unchanged
# Y: no interaction term AW
pi0 <- function(w) expit(0.5+(w[,1]/3))
# mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1] - 0.3*w[,2] + 0.1*w[,3]
# mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1]^2 - 0.3*w[,2] + 0.1*w[,3]
# mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*a*w[,1] - 0.3*w[,2] + 0.1*w[,3]
mu0 <- function(a, w) 0.1 + 0.25*a + 0.75*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

set.seed(95359812)
# W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=runif(n, -1, 1))
W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
A <- rbinom(n, size=1, prob=pi0(W))
Y <- rnorm(n, mean=mu0(A, W), sd=1)

psi0 <- 0.04
theta0 <- 0

psi0 <- 0.2^2 + (0.5^2)/3
theta0 <- (0.5^2)/3

psi0 <- (0.25+0.75*(5/6))^2+(0.75^2*(1/5-1/9))+0.75^2*1/4
theta0 <- (0.75^2*(1/5-1/9))+0.75^2*1/4

# quick check true psi and theta
# psi0 <- c()
# theta0 <- c()
# n <- 2000
# for (i in 1:1000){
#   pi0 <- function(w) expit(0.5+(w[,1]/3))
#   # mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1] - 0.3*w[,2] + 0.1*w[,3]
#   mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*a*w[,1] - 0.3*w[,2] + 0.1*w[,3]
#
#   W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=runif(n, -1, 1))
#   A <- rbinom(n, size=1, prob=pi0(W))
#   Y <- rnorm(n, mean=mu0(A, W), sd=1)
#
#   psi0[i]<- mean((mu0(1, W) - mu0(0, W))^2)
#   theta0[i] <- var((mu0(1, W) - mu0(0, W)))
# }
# mean(psi0)
# mean(theta0)

############
# estimate
###########
# 1) estimated propensity
# a) use glm
AW <- cbind(data.frame(A=A), data.frame(W))
prop.reg <- glm(A ~ ., data = AW, family = "binomial")
pi.hat <- predict(prop.reg, newdata = data.frame(W), type = "response")
# b) use superlearner
n <- length(A)
prop.reg <- SuperLearner(Y=A, X = data.frame(W),
                         newX = data.frame(W),
                         SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                         family = binomial(),
                         obsWeights=rep(1,n),
                         id=1:n)
pi.hat <- prop.reg$SL.predict
# expit(prop.reg$coefficients[1] + prop.reg$coefficients[2]*W[1,1] +
#   prop.reg$coefficients[3]*W[1,2] + prop.reg$coefficients[4]*W[1,3])

# 2) estimated outcome regression
# a) estimate mu1 and mu0 in one model
## a.1 glm
AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

mu.reg <- glm(Y ~ ., data=AW, family='gaussian')
# mu.reg <- glm(Y ~ A + A*W1 + W2 + W3, data=AW, family='gaussian')
mu.reg <- glm(Y ~ .^2, data = AW, family = 'gaussian')
#mu.reg <-  lm(Y ~ A + W)
# mu.hat <- mu.reg$fitted.values
# gam.model <- as.formula("Y ~ s(W1, 2) + s(W2, 2) + s(W1):A + s(W2):A + A*W3")
# fit.gam <- gam::gam(gam.model, data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50), weights = obsWeights)
gam.model <- as.formula("Y ~ s(W1, 2) + s(W2, 2) + s(W1):A + s(W2):A + A*W3")
gam.model <- as.formula("Y ~ s(W1, 1) + s(W2, 5) + s(W1,2):A + s(W2,1):A + A*W3")
gam.model <- as.formula("Y ~ s(W1, 4) + s(W2, 4) + s(W1,4):A + s(W2,4):A + A*W3")
gam.model <- as.formula("Y ~ W1 + I(W2^2) + I(W1^2):A + A*W3")
gam.model <- as.formula("Y ~ s(W1, 2) + lo(W2, span=0.5) + lo(W1, span=0.5):A + lo(W2, span=0.5):A + A*W3") # gam3
gam.model <- as.formula("Y ~ s(W1, 2) + I(W2^2) + I(W1^2):A + s(W1, 2):A + A*W3") # gam4
gam.model <- as.formula("Y ~ s(W1, df=4) + s(W2, df=4) + I(W1^2):A + W2:A + A*W3") # GAM5

mu.reg <- gam::gam(gam.model, data = AW, family = 'gaussian', control = gam::gam.control(maxit = 30, bf.maxit = 30))
plot(c(1:2000), mu.reg$smooth[,2])
points(c(1:2000), W$W2^2)

mu.reg$coefficients
plot(mu.reg)
plot(W$W2, Y)

x <- W$W2
y <- Y
ss1<-smooth.spline(x,y,df=3)
ss2<-smooth.spline(x,y,df=15)
ss<-smooth.spline(x,y)
plot(x,y)
lines(ss1$x,ss1$y,col="blue",lwd=2)
lines(ss2$x,ss2$y,col="yellow",lwd=2,lty=2)
lines(ss$x,ss$y,col="red",lwd=2)

plot(W$W2, mu.reg$smooth[,2])
plot(W$W1, mu.reg$smooth[,1])

gam.s.ret <- gam.s(W$W2, Y, df=4)

mu.reg <- earth::earth(x = data.frame(cbind(A, W)), y = Y, degree = 2, penalty = -1)
mu.reg <- earth::earth(x = data.frame(cbind(A, W)), y = Y, degree = 2, penalty = 8)
mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
names(mu.hats) <- c("mu1", "mu0")

# mean square error
mean((mu.reg$fitted.values - mu0(A, W))^2)
mean((mu1.hat - mu0(1, W))^2)
mean((mu0.hat - mu0(0, W))^2)


## a.2 superlearner
learners = create.Learner("SL.earth", params = list(penalty=-1))
mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                       newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                       SL.library = c("SL.earth"),
                       # SL.library = c("SL.glm.interaction", "SL.gam.interaction", learners$names),
                       family = 'gaussian',
                       obsWeights=rep(1,n),
                       id=1:n)
mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

# mean square error

# mean((mu.reg$fitted.values - mu0(A, W))^2)
mean((mu.hats$mu1 - mu0(1, W))^2)
mean((mu.hats$mu0 - mu0(0, W))^2)

tau.hat <- mu.hats$mu1 - mu.hats$mu0
Z.hat <- (2*A - 1) / (A * pi.hat + (1-A) * (1-pi.hat))
var(tau.hat)

psi.plug.in <- mean(tau.hat^2)
psi.eif.hat <- 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2 - psi.plug.in
psi.se <- sd(psi.eif.hat)
psi.one.step.est <- mean(2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2)
psi.est <- ifelse(psi.one.step.est>0, psi.one.step.est, psi.plug.in)

gamma.hat <- mean(tau.hat)
theta.plug.in <- mean((tau.hat-gamma.hat)^2)
theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - theta.plug.in
theta.se <- sd(theta.eif.hat)
theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
theta.est <- ifelse(theta.one.step.est>0, theta.one.step.est, theta.plug.in)

# b) estimate mu1 and mu0 separately
## b.1 glm
AW1 <- cbind(data.frame(Y=Y[which(A==1)]), data.frame(W[which(A==1),]))
mu1.reg  <- glm(Y ~ ., data=AW1, family='gaussian')
mu1.hat2 <- predict(mu1.reg, newdata = data.frame(W), type = 'response')
AW0 <- cbind(data.frame(Y=Y[which(A==0)]), data.frame(W[which(A==0),]))
mu0.reg  <- glm(Y ~ ., data=AW0, family='gaussian')
mu0.hat2 <- predict(mu0.reg, newdata = data.frame(W), type = 'response')

## b.2 superlearner
mu0.reg <- SuperLearner(Y=Y[which(A==0)], X = data.frame(W[which(A==0),]),
                        newX = data.frame(W),
                        # SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                        SL.library = c("SL.glm", "SL.gam"),
                        family = 'gaussian',
                        obsWeights=rep(1,length(which(A==0))),
                        id=1:length(which(A==0)))
mu0.hat2 <- mu0.reg$SL.predict
mu1.reg <- SuperLearner(Y=Y[which(A==1)], X = data.frame(W[which(A==1),]),
                        newX = data.frame(W),
                        # SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                        SL.library = c("SL.gam", "SL.glm"),
                        family = 'gaussian',
                        obsWeights=rep(1,length(which(A==1))),
                        id=1:length(which(A==1)))
mu1.hat2 <- mu1.reg$SL.predict

mu.hat2 <- A * mu1.hat2 + (1-A) * mu0.hat2
tau.hat2 <- mu1.hat2 - mu0.hat2

psi.plug.in2 <- mean(tau.hat2^2)
psi.eif.hat2 <- 2 * tau.hat2 * (Y - mu.hat2) * Z.hat + tau.hat2^2 - psi.plug.in2
psi.se2 <- sd(psi.eif.hat2)
psi.one.step.est2 <- mean(2 * tau.hat2 * (Y - mu.hat2) * Z.hat + tau.hat2^2)
psi.est2 <- ifelse(psi.one.step.est2>0, psi.one.step.est2, psi.plug.in2)

gamma.hat2 <- mean(tau.hat2)
theta.plug.in2 <- mean((tau.hat2-gamma.hat2)^2)
theta.eif.hat2 <- 2 * (tau.hat2-gamma.hat2) * (Y - mu.hat2) * Z.hat + (tau.hat2-gamma.hat2)^2 - theta.plug.in2
theta.se2 <- sd(theta.eif.hat2)
theta.one.step.est2 <- mean(2 * (tau.hat2-gamma.hat2) * (Y - mu.hat2) * Z.hat + (tau.hat2-gamma.hat2)^2)
theta.est2 <- ifelse(theta.one.step.est2>0, theta.one.step.est2, theta.plug.in2)



# confidence interval

##############
# simulation
#############
load("est.sim.rst.redo1.RData")
seed.fixed <- est.sim.rst.redo1$seed
seed.list <- c()
psi.one.step.est.list <- c()
psi.plug.in.list <- c()
psi.se.list <- c()
psi.est.list <- c()
theta.one.step.est.list <- c()
theta.plug.in.list <- c()
theta.se.list <- c()
theta.est.list <- c()
psi.one.step.est2.list <- c()
psi.plug.in2.list <- c()
psi.se2.list <- c()
psi.est2.list <- c()
theta.one.step.est2.list <- c()
theta.plug.in2.list <- c()
theta.se2.list <- c()
theta.est2.list <- c()
for (j in 1:1000){
  n <- 2000
  # seed <- sample(1e3:1e8, 1)
  seed <- seed.fixed[j]
  cat(n, j, seed, '\n')
  set.seed(seed)

  #################
  # generate data
  #################
  # W: continuous
  # A: unchanged
  # Y: no interaction term AW
  pi0 <- function(w) expit(0.5+(w[,1]/3))
  # mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1] - 0.3*w[,2] + 0.1*w[,3]
  # mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1]^2 - 0.3*w[,2] + 0.1*w[,3]
  # mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*a*w[,1] - 0.3*w[,2] + 0.1*w[,3]
  mu0 <- function(a, w) 0.1 + 0.25*a + 0.75*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2

  # W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=runif(n, -1, 1))
  W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
  A <- rbinom(n, size=1, prob=pi0(W))
  Y <- rnorm(n, mean=mu0(A, W), sd=1)

  # psi0 <- 0.04
  # theta0 <- 0

  # psi0 <- 0.2^2 + (0.5^2)/3
  # theta0 <- (0.5^2)/3

  psi0 <- (0.25+0.75*(5/6))^2+(0.75^2*(1/5-1/9))+0.75^2*1/4
  theta0 <- (0.75^2*(1/5-1/9))+0.75^2*1/4

  ############
  # estimate
  ###########
  # 1) estimated propensity
  # a) use glm
  AW <- cbind(data.frame(A=A), data.frame(W))
  prop.reg <- glm(A ~ ., data = AW, family = "binomial")
  pi.hat <- predict(prop.reg, newdata = data.frame(W), type = "response")
  # b) use superlearner
  n <- length(A)
  # prop.reg <- SuperLearner(Y=A, X = data.frame(W),
  #                          newX = data.frame(W),
  #                          SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
  #                          family = binomial(),
  #                          obsWeights=rep(1,n),
  #                          id=1:n)
  # pi.hat <- prop.reg$SL.predict

  # 2) estimated outcome regression
  # a) estimate mu1 and mu0 in one model
  ## a.1 glm
  AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

  # mu.reg  <- glm(Y ~ ., data=AW, family='gaussian')
  # mu.reg  <- glm(Y ~ A + A*W1 + W2 + W3, data=AW, family='gaussian')
  # mu.reg <- glm(Y ~ .^2, data = AW, family = 'gaussian')
  gam.model <- as.formula("Y ~ s(W1, df=4) + s(W2, df=4) + I(W1^2):A + W2:A + A*W3")
  mu.reg <- gam::gam(gam.model, data = AW, family = 'gaussian', control = gam::gam.control(maxit = 50, bf.maxit = 50))
  #mu.reg <-  lm(Y ~ A + W)
  # mu.hat <- mu.reg$fitted.values
  mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
  mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
  mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)

  ## a.2 superlearner
  # learners = create.Learner("SL.earth", params = list(penalty=-1))
  # mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
  #                       newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
  #                       # SL.library = c("SL.earth"),
  #                       SL.library = c("SL.earth", learners$names),
  #                       family = 'gaussian',
  #                       obsWeights=rep(1,n),
  #                       id=1:n)
  # mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
  mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

  tau.hat <- mu.hats$mu1 - mu.hats$mu0
  Z.hat <- (2*A - 1) / (A * pi.hat + (1-A) * (1-pi.hat))

  psi.plug.in <- mean(tau.hat^2)
  psi.eif.hat <- 2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2 - psi.plug.in
  psi.se <- sd(psi.eif.hat)
  psi.one.step.est <- mean(2 * tau.hat * (Y - mu.hat) * Z.hat + tau.hat^2)
  psi.est <- ifelse(psi.one.step.est>0, psi.one.step.est, psi.plug.in)

  gamma.hat <- mean(tau.hat)
  theta.plug.in <- mean((tau.hat-gamma.hat)^2)
  theta.eif.hat <- 2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2 - theta.plug.in
  theta.se <- sd(theta.eif.hat)
  theta.one.step.est <- mean(2 * (tau.hat-gamma.hat) * (Y - mu.hat) * Z.hat + (tau.hat-gamma.hat)^2)
  theta.est <- ifelse(theta.one.step.est>0, theta.one.step.est, theta.plug.in)

  # b) estimate mu1 and mu0 separately
  ## b.1 glm
  AW1 <- cbind(data.frame(Y=Y[which(A==1)]), data.frame(W[which(A==1),]))
  mu1.reg  <- glm(Y ~ ., data=AW1, family='gaussian')
  mu1.hat2 <- predict(mu1.reg, newdata = data.frame(W), type = 'response')
  AW0 <- cbind(data.frame(Y=Y[which(A==0)]), data.frame(W[which(A==0),]))
  mu0.reg  <- glm(Y ~ ., data=AW0, family='gaussian')
  mu0.hat2 <- predict(mu0.reg, newdata = data.frame(W), type = 'response')

  ## b.2 superlearner
  # mu0.reg <- SuperLearner(Y=Y[which(A==0)], X = data.frame(W[which(A==0),]),
  #                         newX = data.frame(W),
  #                         # SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
  #                         SL.library = c("SL.glm"),
  #                         family = 'gaussian',
  #                         obsWeights=rep(1,length(which(A==0))),
  #                         id=1:length(which(A==0)))
  # mu0.hat2 <- mu0.reg$SL.predict
  # mu1.reg <- SuperLearner(Y=Y[which(A==1)], X = data.frame(W[which(A==1),]),
  #                         newX = data.frame(W),
  #                         # SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
  #                         SL.library = c("SL.glm"),
  #                         family = 'gaussian',
  #                         obsWeights=rep(1,length(which(A==1))),
  #                         id=1:length(which(A==1)))
  # mu1.hat2 <- mu1.reg$SL.predict

  mu.hat2 <- A * mu1.hat2 + (1-A) * mu0.hat2
  tau.hat2 <- mu1.hat2 - mu0.hat2

  psi.plug.in2 <- mean(tau.hat2^2)
  psi.eif.hat2 <- 2 * tau.hat2 * (Y - mu.hat2) * Z.hat + tau.hat2^2 - psi.plug.in2
  psi.se2 <- sd(psi.eif.hat2)
  psi.one.step.est2 <- mean(2 * tau.hat2 * (Y - mu.hat2) * Z.hat + tau.hat2^2)
  psi.est2 <- ifelse(psi.one.step.est2>0, psi.one.step.est2, psi.plug.in2)

  gamma.hat2 <- mean(tau.hat2)
  theta.plug.in2 <- mean((tau.hat2-gamma.hat2)^2)
  theta.eif.hat2 <- 2 * (tau.hat2-gamma.hat2) * (Y - mu.hat2) * Z.hat + (tau.hat2-gamma.hat2)^2 - theta.plug.in2
  theta.se2 <- sd(theta.eif.hat2)
  theta.one.step.est2 <- mean(2 * (tau.hat2-gamma.hat2) * (Y - mu.hat2) * Z.hat + (tau.hat2-gamma.hat2)^2)
  theta.est2 <- ifelse(theta.one.step.est2>0, theta.one.step.est2, theta.plug.in2)

  seed.list[j] <- seed
  psi.one.step.est.list[j] <- psi.one.step.est
  psi.plug.in.list[j] <- psi.plug.in
  psi.se.list[j] <- psi.se
  psi.est.list[j] <- psi.est
  theta.one.step.est.list[j] <- theta.one.step.est
  theta.plug.in.list[j] <- theta.plug.in
  theta.se.list[j] <- theta.se
  theta.est.list[j] <- theta.est
  psi.one.step.est2.list[j] <- psi.one.step.est2
  psi.plug.in2.list[j] <- psi.plug.in2
  psi.se2.list[j] <- psi.se2
  psi.est2.list[j] <- psi.est2
  theta.one.step.est2.list[j] <- theta.one.step.est2
  theta.plug.in2.list[j] <- theta.plug.in2
  theta.se2.list[j] <- theta.se2
  theta.est2.list[j] <- theta.est2
}

# save simulation data

est.sim.rst.redo4.gam5 <- data.frame(
  seed = seed.list,
  n = rep(2000, 1000),
  psi.one.step.est = psi.one.step.est.list,
  psi.plug.in = psi.plug.in.list,
  psi.se = psi.se.list,
  psi.est = psi.est.list,
  theta.one.step.est = theta.one.step.est.list,
  theta.plug.in = theta.plug.in.list,
  theta.se = theta.se.list,
  theta.est = theta.est.list,
  psi.one.step.est2 = psi.one.step.est2.list,
  psi.plug.in2 = psi.plug.in2.list,
  psi.se2 = psi.se2.list,
  psi.est2 = psi.est2.list,
  theta.one.step.est2 = theta.one.step.est2.list,
  theta.plug.in2 = theta.plug.in2.list,
  theta.se2 = theta.se2.list,
  theta.est2 = theta.est2.list
)
save(est.sim.rst.redo4.8, file="est.sim.rst.redo4.8.RData")

# analysis
rst <- est.sim.rst.redo4.gam5
psi0
# 1) bias
# psi
sqrt(n)*mean(rst$psi.plug.in - psi0)
sqrt(n)*mean(rst$psi.one.step.est - psi0)
sqrt(n)*mean(rst$psi.est - psi0)
sqrt(n)*mean(rst$psi.plug.in2 - psi0)
sqrt(n)*mean(rst$psi.one.step.est2 - psi0)
sqrt(n)*mean(rst$psi.est2 - psi0)
# theta
sqrt(n)*mean(rst$theta.plug.in - theta0)
sqrt(n)*mean(rst$theta.one.step.est - theta0)
sqrt(n)*mean(rst$theta.est - theta0)
sqrt(n)*mean(rst$theta.plug.in2 - theta0)
sqrt(n)*mean(rst$theta.one.step.est2 - theta0)
sqrt(n)*mean(rst$theta.est2 - theta0)
# 2) se
# empirical.sd
sd(sqrt(n)*(rst$psi.plug.in - psi0))
sd(sqrt(n)*(rst$psi.est - psi0))
mean(rst$psi.se)

sd(sqrt(n)*(rst$psi.plug.in2 - psi0))
sd(sqrt(n)*(rst$psi.est2 - psi0))
mean(rst$psi.se2)

sd(sqrt(n)*(rst$theta.plug.in - theta0))
sd(sqrt(n)*(rst$theta.est - theta0))
mean(rst$theta.se)

sd(sqrt(n)*(rst$theta.plug.in2 - theta0))
sd(sqrt(n)*(rst$theta.est2 - theta0))
mean(rst$theta.se2)

# mean(se)

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
  gam.model <- as.formula("Y ~ s(W1, 4) + s(W2, 4) + s(W1, 4):A + s(W2, 4):A + A*W3")
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
