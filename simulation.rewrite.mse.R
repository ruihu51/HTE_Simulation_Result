rm(list=ls())

library(SuperLearner)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mgcv)

# utils
expit <- function(x) 1 / (1 + exp(-x))


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
mse.mu.list <- c()
mse.mu1.list <- c()
mse.mu0.list <- c()
mse.tau.list <- c()
n.list <- c()
j.list <- c()

# s1 glm simulation
rst.s1.glm <- data.frame()
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)){
  for (j in 1:1000){
    # n <- 2000
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
    mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1] - 0.3*w[,2] + 0.1*w[,3] # scenario1

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=runif(n, -1, 1))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W), sd=1)

    psi0 <- 0.04
    theta0 <- 0

    ############
    # estimate
    ###########
    # 1) estimated propensity
    AW <- cbind(data.frame(A=A), data.frame(W))
    prop.reg <- glm(A ~ ., data = AW, family = "binomial")
    pi.hat <- predict(prop.reg, newdata = data.frame(W), type = "response")

    # 2) estimated outcome regression
    AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

    mu.reg  <- glm(Y ~ ., data=AW, family='gaussian')
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
    mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)

    mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

    # mean square error
    mse.mu <- mean((mu.hat - mu0(A, W))^2)
    mse.mu1 <- mean((mu.hats$mu1 - mu0(1, W))^2)
    mse.mu0 <- mean((mu.hats$mu0 - mu0(0, W))^2)

    tau.hat <- mu.hats$mu1 - mu.hats$mu0
    tau0 <- mu0(1, W) - mu0(0, W)
    mse.tau <- mean((tau.hat - tau0)^2)

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

    seed.list[j] <- seed
    psi.one.step.est.list[j] <- psi.one.step.est
    psi.plug.in.list[j] <- psi.plug.in
    psi.se.list[j] <- psi.se
    psi.est.list[j] <- psi.est
    theta.one.step.est.list[j] <- theta.one.step.est
    theta.plug.in.list[j] <- theta.plug.in
    theta.se.list[j] <- theta.se
    theta.est.list[j] <- theta.est
    mse.mu.list[j] <- mse.mu
    mse.mu1.list[j] <- mse.mu1
    mse.mu0.list[j] <- mse.mu0
    mse.tau.list[j] <- mse.tau
    n.list[j] <- n
    j.list[j] <- j
  }
  d1 <- data.frame(
    seed = seed.list,
    n = n.list,
    j = j.list,
    psi.one.step.est = psi.one.step.est.list,
    psi.plug.in = psi.plug.in.list,
    psi.se = psi.se.list,
    psi.est = psi.est.list,
    theta.one.step.est = theta.one.step.est.list,
    theta.plug.in = theta.plug.in.list,
    theta.se = theta.se.list,
    theta.est = theta.est.list,
    mse.mu = mse.mu.list,
    mse.mu1 = mse.mu1.list,
    mse.mu0 = mse.mu0.list,
    mse.tau = mse.tau.list
  )
  rst.s1.glm <- rbind(rst.s1.glm, d1)
}

# s1 earth simulation
rst.s1.earth <- data.frame()
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)){
  for (j in 1:1000){
    # n <- 2000
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
    mu0 <- function(a, w) 0.1 + 0.2*a + 0.5*w[,1] - 0.3*w[,2] + 0.1*w[,3]

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=runif(n, -1, 1))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W), sd=1)

    psi0 <- 0.04
    theta0 <- 0

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

    # 2) estimated outcome regression
    # a) estimate mu1 and mu0 in one model
    ## a.1 glm
    AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

    ## a.2 superlearner
    mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                           newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                           SL.library = c("SL.earth"),
                           family = 'gaussian',
                           obsWeights=rep(1,n),
                           id=1:n)
    mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
    mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

    # mean square error
    mse.mu <- mean((mu.hat - mu0(A, W))^2)
    mse.mu1 <- mean((mu.hats$mu1 - mu0(1, W))^2)
    mse.mu0 <- mean((mu.hats$mu0 - mu0(0, W))^2)

    tau.hat <- mu.hats$mu1 - mu.hats$mu0
    tau0 <- mu0(1, W) - mu0(0, W)
    mse.tau <- mean((tau.hat - tau0)^2)

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

    seed.list[j] <- seed
    psi.one.step.est.list[j] <- psi.one.step.est
    psi.plug.in.list[j] <- psi.plug.in
    psi.se.list[j] <- psi.se
    psi.est.list[j] <- psi.est
    theta.one.step.est.list[j] <- theta.one.step.est
    theta.plug.in.list[j] <- theta.plug.in
    theta.se.list[j] <- theta.se
    theta.est.list[j] <- theta.est
    mse.mu.list[j] <- mse.mu
    mse.mu1.list[j] <- mse.mu1
    mse.mu0.list[j] <- mse.mu0
    mse.tau.list[j] <- mse.tau
    n.list[j] <- n
    j.list[j] <- j
  }
  d1 <- data.frame(
    seed = seed.list,
    n = n.list,
    j = j.list,
    psi.one.step.est = psi.one.step.est.list,
    psi.plug.in = psi.plug.in.list,
    psi.se = psi.se.list,
    psi.est = psi.est.list,
    theta.one.step.est = theta.one.step.est.list,
    theta.plug.in = theta.plug.in.list,
    theta.se = theta.se.list,
    theta.est = theta.est.list,
    mse.mu = mse.mu.list,
    mse.mu1 = mse.mu1.list,
    mse.mu0 = mse.mu0.list,
    mse.tau = mse.tau.list
  )
  rst.s1.earth <- rbind(rst.s1.earth, d1)
}

save(rst.s1.glm, file = "rst.s1.glm.RData")
save(rst.s1.earth, file = "rst.s1.earth.RData")

# s3 gam simulation
load("Data/Testing/testing.sim.2.25.75.sl3.Rdata")
rst.s3.gam.correct <- data.frame()
create.seed.dict <- function(x){
  x %>%
    subset(type %in% "Gamma.stat") %>%
    select(c("n", "j", "seed")) %>%
    rename(sample_n=n, j=j, seed=seed)
}
seed.data <- create.seed.dict(testing.sim.2.25.75.sl3)
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)){
  for (j in 1:1000){
    # n <- 2000
    # seed <- sample(1e3:1e8, 1)
    # seed <- seed.fixed[j]
    seed <- subset(seed.data, (sample_n==n))[j,"seed"]
    cat(n, j, seed, '\n')
    set.seed(seed)

    #################
    # generate data
    #################
    # W: continuous
    # A: unchanged
    # Y: no interaction term AW
    pi0 <- function(w) expit(0.5+(w[,1]/3))
    mu0 <- function(a, w) 0.1 + 0.25*a + 0.75*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2 #3

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W), sd=1)

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

    # 2) estimated outcome regression
    # a) estimate mu1 and mu0 in one model
    ## a.1 glm
    AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

    gam.model <- as.formula("Y ~ W1 + I(W2^2) + I(W1^2):A + A*W3")
    mu.reg <- gam::gam(gam.model, data = AW, family = 'gaussian', control = gam::gam.control(maxit = 50, bf.maxit = 50))

    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
    mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
    names(mu.hats) <- c("mu1", "mu0")

    ## a.2 superlearner
    # mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
    #                        newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
    #                        SL.library = c("SL.earth"),
    #                        # SL.library = c("SL.earth", learners$names),
    #                        family = 'gaussian',
    #                        obsWeights=rep(1,n),
    #                        id=1:n)
    # mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
    mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

    # mean square error
    mse.mu <- mean((mu.hat - mu0(A, W))^2)
    mse.mu1 <- mean((mu.hats$mu1 - mu0(1, W))^2)
    mse.mu0 <- mean((mu.hats$mu0 - mu0(0, W))^2)

    tau.hat <- mu.hats$mu1 - mu.hats$mu0
    tau0 <- mu0(1, W) - mu0(0, W)
    mse.tau <- mean((tau.hat - tau0)^2)

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

    seed.list[j] <- seed
    psi.one.step.est.list[j] <- psi.one.step.est
    psi.plug.in.list[j] <- psi.plug.in
    psi.se.list[j] <- psi.se
    psi.est.list[j] <- psi.est
    theta.one.step.est.list[j] <- theta.one.step.est
    theta.plug.in.list[j] <- theta.plug.in
    theta.se.list[j] <- theta.se
    theta.est.list[j] <- theta.est
    mse.mu.list[j] <- mse.mu
    mse.mu1.list[j] <- mse.mu1
    mse.mu0.list[j] <- mse.mu0
    mse.tau.list[j] <- mse.tau
    n.list[j] <- n
    j.list[j] <- j
  }
  d1 <- data.frame(
    seed = seed.list,
    n = n.list,
    j = j.list,
    psi.one.step.est = psi.one.step.est.list,
    psi.plug.in = psi.plug.in.list,
    psi.se = psi.se.list,
    psi.est = psi.est.list,
    theta.one.step.est = theta.one.step.est.list,
    theta.plug.in = theta.plug.in.list,
    theta.se = theta.se.list,
    theta.est = theta.est.list,
    mse.mu = mse.mu.list,
    mse.mu1 = mse.mu1.list,
    mse.mu0 = mse.mu0.list,
    mse.tau = mse.tau.list
  )
  rst.s3.gam.correct <- rbind(rst.s3.gam.correct, d1)
}

rst.s3.earth <- data.frame()
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)){
  for (j in 1:1000){
    # n <- 2000
    # seed <- sample(1e3:1e8, 1)
    # seed <- seed.fixed[j]
    seed <- subset(seed.data, (sample_n==n))[j,"seed"]
    cat(n, j, seed, '\n')
    set.seed(seed)

    #################
    # generate data
    #################
    # W: continuous
    # A: unchanged
    # Y: no interaction term AW
    pi0 <- function(w) expit(0.5+(w[,1]/3))
    mu0 <- function(a, w) 0.1 + 0.25*a + 0.75*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2 #3

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W), sd=1)

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

    # 2) estimated outcome regression
    # a) estimate mu1 and mu0 in one model
    ## a.1 glm
    AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))

    ## a.2 superlearner
    mu.reg <- SuperLearner(Y=Y, X = data.frame(cbind(A, W)),
                           newX = rbind(data.frame(cbind(A=1, W)), data.frame(cbind(A=0, W))),
                           SL.library = c("SL.earth"),
                           # SL.library = c("SL.earth", learners$names),
                           family = 'gaussian',
                           obsWeights=rep(1,n),
                           id=1:n)
    mu.hats <- data.frame(mu1=mu.reg$SL.predict[1:n], mu0=mu.reg$SL.predict[-(1:n)])
    names(mu.hats) <- c("mu1", "mu0")
    mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

    # mean square error
    mse.mu <- mean((mu.hat - mu0(A, W))^2)
    mse.mu1 <- mean((mu.hats$mu1 - mu0(1, W))^2)
    mse.mu0 <- mean((mu.hats$mu0 - mu0(0, W))^2)

    tau.hat <- mu.hats$mu1 - mu.hats$mu0
    tau0 <- mu0(1, W) - mu0(0, W)
    mse.tau <- mean((tau.hat - tau0)^2)

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

    seed.list[j] <- seed
    psi.one.step.est.list[j] <- psi.one.step.est
    psi.plug.in.list[j] <- psi.plug.in
    psi.se.list[j] <- psi.se
    psi.est.list[j] <- psi.est
    theta.one.step.est.list[j] <- theta.one.step.est
    theta.plug.in.list[j] <- theta.plug.in
    theta.se.list[j] <- theta.se
    theta.est.list[j] <- theta.est
    mse.mu.list[j] <- mse.mu
    mse.mu1.list[j] <- mse.mu1
    mse.mu0.list[j] <- mse.mu0
    mse.tau.list[j] <- mse.tau
    n.list[j] <- n
    j.list[j] <- j
  }
  d1 <- data.frame(
    seed = seed.list,
    n = n.list,
    j = j.list,
    psi.one.step.est = psi.one.step.est.list,
    psi.plug.in = psi.plug.in.list,
    psi.se = psi.se.list,
    psi.est = psi.est.list,
    theta.one.step.est = theta.one.step.est.list,
    theta.plug.in = theta.plug.in.list,
    theta.se = theta.se.list,
    theta.est = theta.est.list,
    mse.mu = mse.mu.list,
    mse.mu1 = mse.mu1.list,
    mse.mu0 = mse.mu0.list,
    mse.tau = mse.tau.list
  )
  rst.s3.earth <- rbind(rst.s3.earth, d1)
}

rst.s3.gam.mgcv <- data.frame()
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)){
  for (j in 1:1000){
    # n <- 2000
    # seed <- sample(1e3:1e8, 1)
    # seed <- seed.fixed[j]
    seed <- subset(seed.data, (sample_n==n))[j,"seed"]
    cat(n, j, seed, '\n')
    set.seed(seed)

    #################
    # generate data
    #################
    # W: continuous
    # A: unchanged
    # Y: no interaction term AW
    pi0 <- function(w) expit(0.5+(w[,1]/3))
    mu0 <- function(a, w) 0.1 + 0.25*a + 0.75*a*(w[,1]^2+w[,3]) + w[,1] + w[,2]^2 #3

    W <- data.frame(W1=runif(n, -1, 1), W2=runif(n, -1, 1), W3=rbinom(n, size=1, prob=0.5))
    A <- rbinom(n, size=1, prob=pi0(W))
    Y <- rnorm(n, mean=mu0(A, W), sd=1)

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

    # 2) estimated outcome regression
    # a) estimate mu1 and mu0 in one model
    ## a.1 glm
    AW <- cbind(data.frame(Y=Y, A=A), data.frame(W))
    gam.model <- as.formula("Y ~ s(W1) + s(W2) + s(W1, by=A) + s(W2, by=A) + A*W3")
    mu.reg <- mgcv::gam(gam.model, data = AW, method = "REML")
    mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1), data.frame(W)), type = 'response')
    mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0), data.frame(W)), type = 'response')
    mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
    names(mu.hats) <- c("mu1", "mu0")

    ## a.2 superlearner
    mu.hat <- A * mu.hats$mu1 + (1-A) * mu.hats$mu0

    # mean square error
    mse.mu <- mean((mu.hat - mu0(A, W))^2)
    mse.mu1 <- mean((mu.hats$mu1 - mu0(1, W))^2)
    mse.mu0 <- mean((mu.hats$mu0 - mu0(0, W))^2)

    tau.hat <- mu.hats$mu1 - mu.hats$mu0
    tau0 <- mu0(1, W) - mu0(0, W)
    mse.tau <- mean((tau.hat - tau0)^2)

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

    seed.list[j] <- seed
    psi.one.step.est.list[j] <- psi.one.step.est
    psi.plug.in.list[j] <- psi.plug.in
    psi.se.list[j] <- psi.se
    psi.est.list[j] <- psi.est
    theta.one.step.est.list[j] <- theta.one.step.est
    theta.plug.in.list[j] <- theta.plug.in
    theta.se.list[j] <- theta.se
    theta.est.list[j] <- theta.est
    mse.mu.list[j] <- mse.mu
    mse.mu1.list[j] <- mse.mu1
    mse.mu0.list[j] <- mse.mu0
    mse.tau.list[j] <- mse.tau
    n.list[j] <- n
    j.list[j] <- j
  }
  d1 <- data.frame(
    seed = seed.list,
    n = n.list,
    j = j.list,
    psi.one.step.est = psi.one.step.est.list,
    psi.plug.in = psi.plug.in.list,
    psi.se = psi.se.list,
    psi.est = psi.est.list,
    theta.one.step.est = theta.one.step.est.list,
    theta.plug.in = theta.plug.in.list,
    theta.se = theta.se.list,
    theta.est = theta.est.list,
    mse.mu = mse.mu.list,
    mse.mu1 = mse.mu1.list,
    mse.mu0 = mse.mu0.list,
    mse.tau = mse.tau.list
  )
  rst.s3.gam.mgcv <- rbind(rst.s3.gam.mgcv, d1)
}

save(rst.s3.gam.correct, file = "rst.s3.gam.correct.RData")
save(rst.s3.earth, file = "rst.s3.earth.RData")
save(rst.s3.gam.mgcv, file = "rst.s3.gam.mgcv.RData")

# plot
# s1 mse
load("rst.s1.glm.RData")
load("rst.s1.earth.RData")
rst.s1 <- rbind(rst.s1.glm %>%
                  gather(mu.type, mse, mse.mu:mse.tau) %>%
                  mutate(est.method="glm"),
                rst.s1.earth %>%
                  gather(mu.type, mse, mse.mu:mse.tau) %>%
                  mutate(est.method="earth"))

rst.s1 <- rst.s1 %>%
  mutate(normL2=sqrt(mse*n)) %>%
  mutate(rmse=normL2/sqrt(n))

rst.s1 %>%
  ggplot() +
  geom_boxplot(aes(as.factor(n), y = rmse, color=est.method)) +
  facet_wrap(~mu.type, scales='free_x')

rst.s1 %>%
  filter((est.method=="glm") & (mu.type=="mse.mu")) %>%
  group_by(n) %>%
  summarise(n())

for (ml in c("glm", "earth")){
  reg.data <- rst.s1 %>%
    filter((est.method==ml) & (mu.type=="mse.mu"))

  reg<-lm(formula = log(rmse) ~ log(n),
          data=reg.data)
  coeff <- coefficients(reg)
  intercept <- coeff[1]
  slope <- coeff[2]

  rst.s1 %>%
    filter((est.method==ml) & (mu.type=="mse.mu")) %>%
    ggplot() +
    geom_point(aes(log(n), log(rmse))) +
    geom_abline(intercept = intercept, slope = slope, color="red",
                linetype="dashed", size=1.5) +
    ggtitle(paste(ml,as.character(slope)))
}

for (ml in c("glm", "earth")){
  reg.data <- rst.s1 %>%
    filter((est.method==ml) & (mu.type=="mse.tau"))

  reg<-lm(formula = log(rmse) ~ log(n),
          data=reg.data)
  coeff <- coefficients(reg)
  intercept <- coeff[1]
  slope <- coeff[2]

  p <- rst.s1 %>%
    filter((est.method==ml) & (mu.type=="mse.tau")) %>%
    ggplot() +
    geom_point(aes(log(n), log(rmse))) +
    geom_abline(intercept = intercept, slope = slope, color="red",
                linetype="dashed", size=1.5) +
    ggtitle(paste(ml,as.character(slope)))

  assign(paste0("plot_", ml), p)
}

plot_glm
plot_earth

# s1 bias
psi0 <- 0.04
theta0 <- 0

rst.s1.bias.psi <- rbind(rst.s1.glm %>%
                           gather(est.type, est, psi.one.step.est:psi.plug.in) %>%
                           mutate(est.method="glm"),
                         rst.s1.earth %>%
                           gather(est.type, est, psi.one.step.est:psi.plug.in) %>%
                           mutate(est.method="earth"))

rst.s1.bias.psi %>%
  group_by(n, est.type, est.method) %>%
  summarise(bias=mean(est - psi0, na.rm=TRUE)) %>%
  ggplot() +
  geom_line(aes(n, sqrt(n) * bias, color=est.type)) +
  facet_wrap(~est.method, scales='free_x') +
  theme(legend.key.size = unit(0.1,"line"))

# s3 mse
load("rst.s3.earth.RData")
load("rst.s3.gam.correct.RData")
load("rst.s3.gam.mgcv.RData")

rst.s3.gam.correct %>%
  head()
rst.s3 <- bind_rows(rst.s3.gam.correct %>%
                  gather(mu.type, mse, mse.mu:mse.tau) %>%
                  mutate(est.method="gam.correct"),
                rst.s3.earth %>%
                  gather(mu.type, mse, mse.mu:mse.tau) %>%
                  mutate(est.method="earth"),
                rst.s3.gam.mgcv %>%
                  gather(mu.type, mse, mse.mu:mse.tau) %>%
                  mutate(est.method="gam.mgcv"))

rst.s3 <- rst.s3 %>%
  mutate(normL2=sqrt(mse*n)) %>%
  mutate(rmse=normL2/sqrt(n))

rst.s3 %>%
  ggplot() +
  geom_boxplot(aes(as.factor(n), y = rmse, color=est.method)) +
  facet_wrap(~mu.type, scales='free_x')

plot.list <- c()
for (ml in c("gam.correct", "earth", "gam.mgcv")){
  reg.data <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.mu"))

  reg<-lm(formula = log(rmse) ~ log(n),
          data=reg.data)
  coeff <- coefficients(reg)
  intercept <- coeff[1]
  slope <- coeff[2]


  p <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.mu")) %>%
    ggplot() +
    geom_point(aes(log(n), log(rmse))) +
    geom_abline(intercept = intercept, slope = slope, color="red",
                linetype="dashed", size=1.5) +
    ggtitle(paste(ml,as.character(slope)))

  assign(paste0("plot_", ml), p)
}
plot_earth
plot_gam.correct
plot_gam.mgcv

plot.list <- c()
for (ml in c("gam.correct", "earth", "gam.mgcv")){
  reg.data <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.tau"))

  reg<-lm(formula = log(rmse) ~ log(n),
          data=reg.data)
  coeff <- coefficients(reg)
  intercept <- coeff[1]
  slope <- coeff[2]


  p <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.tau")) %>%
    ggplot() +
    geom_point(aes(log(n), log(rmse))) +
    geom_abline(intercept = intercept, slope = slope, color="red",
                linetype="dashed", size=1.5) +
    ggtitle(paste(ml,as.character(slope)))

  assign(paste0("plot_", ml), p)
}
plot_earth
plot_gam.correct
plot_gam.mgcv

rst.s3 %>%
  ggplot() +
  geom_boxplot(aes(as.factor(n), y = mse, color=est.method)) +
  facet_wrap(~mu.type, scales='free_x')

for (ml in c("gam.correct", "earth", "gam.mgcv")){
  reg.data <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.mu0"))

  reg<-lm(formula = log(rmse) ~ log(n),
          data=reg.data)
  coeff <- coefficients(reg)
  intercept <- coeff[1]
  slope <- coeff[2]

  p <- rst.s3 %>%
    filter((est.method==ml) & (mu.type=="mse.mu0")) %>%
    ggplot() +
    geom_point(aes(log(n), log(rmse))) +
    geom_abline(intercept = intercept, slope = slope, color="red",
                linetype="dashed", size=1.5) +
    ggtitle(paste(ml,as.character(round(slope, 4))))

  assign(paste0("plot_", ml), p)
}

library(gridExtra)
grid.arrange(plot_earth, plot_gam.correct, plot_gam.mgcv, ncol=1)



# s3 bias
psi0 <- (0.25+0.75*(5/6))^2+(0.75^2*(1/5-1/9))+0.75^2*1/4
theta0 <- (0.75^2*(1/5-1/9))+0.75^2*1/4

rst.s3.bias.psi <- bind_rows(rst.s3.gam.correct %>%
                           gather(est.type, est, psi.one.step.est:psi.plug.in) %>%
                           mutate(est.method="gam.correct"),
                         rst.s3.earth %>%
                           gather(est.type, est, psi.one.step.est:psi.plug.in) %>%
                           mutate(est.method="earth"),
                         rst.s3.gam.mgcv %>%
                           gather(est.type, est, psi.one.step.est:psi.plug.in) %>%
                           mutate(est.method="gam.mgcv"))

rst.s3.bias.psi %>%
  group_by(n, est.type, est.method) %>%
  summarise(bias=mean(est - psi0, na.rm=TRUE),
            abs.bias=mean(abs(est - psi0), na.rm=TRUE)) %>%
  ggplot() +
  geom_line(aes(n, sqrt(n) * bias, color=est.type)) +
  geom_line(aes(n, sqrt(n) * abs.bias, color=est.type)) +
  facet_wrap(~est.method, scales='free_x') +
  theme(legend.key.size = unit(0.1,"line"))

rst.s3.bias.theta <- bind_rows(rst.s3.gam.correct %>%
                           gather(est.type, est, theta.one.step.est:theta.plug.in) %>%
                           mutate(est.method="gam.correct"),
                         rst.s3.earth %>%
                           gather(est.type, est, theta.one.step.est:theta.plug.in) %>%
                           mutate(est.method="earth"),
                         rst.s3.gam.mgcv %>%
                           gather(est.type, est, theta.one.step.est:theta.plug.in) %>%
                           mutate(est.method="gam.mgcv"))

rst.s3.bias.theta %>%
  group_by(n, est.type, est.method) %>%
  summarise(bias=mean(est - theta0, na.rm=TRUE)) %>%
  ggplot() +
  geom_line(aes(n, sqrt(n) * bias, color=est.type)) +
  facet_wrap(~est.method, scales='free_x') +
  theme(legend.key.size = unit(0.1,"line"))

