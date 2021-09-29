rm(list=ls())

library(devtools)
# load_all()
devtools::install_github("ruihu51/CausalSim")


require(CausalSim)

# unit test
null.sims <- FALSE

n <- 1000
set.seed(39338064)
W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
if(null.sims) {
  Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
  psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
  theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
} else {
  Y <- rbinom(n, size = 1, prob = mu0(A, W))
  psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
  theta0 <- var((mu0(1,W) - mu0(0,W)))
}


ret1 <- htem.estimator(A, W, Y, control = list())
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE))
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type = 'boot'))

# use `glm` for both pi.hat and mu.hat in simple unit test
prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
pi.hat <- prop.reg$fitted.values
AW <- cbind(A, data.frame(W))
mu.reg <- glm(Y ~ ., data=AW, family='binomial')
mu.hat <- mu.reg$fitted.values
mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                           data.frame(W)), type = 'response')
mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                           data.frame(W)), type = 'response')
mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
ret1 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE,
                                       pi.hat = pi.hat, mu.hats = mu.hats))
ret2 <- htem.estimator(A, W, Y, control = list(conf.int = TRUE, conf.int.type='boot',
                                               pi.hat = pi.hat, mu.hats = mu.hats))
theta0
# why bootstrap CI didn't work?
theta.test.vector <- rep(0, 500)
for (i in 1:500){
  boot.inds <- sample(1:n, n, replace=TRUE)

  boot.pi.hat <- pi.hat[boot.inds]
  boot.mu.hats <- mu.hats[boot.inds,]
  boot.ret <- htem.estimator(A[boot.inds], W[boot.inds], Y[boot.inds],
                             control = list(pi.hat = boot.pi.hat,
                                            mu.hats = boot.mu.hats,
                                            conf.int = FALSE))
  theta.test.vector[i] <-boot.ret$est[2]
}
density(theta.test.vector)
quantile(theta.test.vector, c(0.025,0.975))
sort(theta.test.vector)
summary(theta.test.vector)
?density

# numerical analysis
## simulation function
ests.sim <-  function(n_range, j_range, control, null.sims=FALSE, out.glm=TRUE){
  ests <- ldply(n_range, function(n) {
    ldply(j_range, function(j) {
      # print(j)
      if(j %% 100 == 0) cat(n, j, '\n')
      seed <- sample(1e3:1e8, 1)
      set.seed(seed)

      W <- matrix(runif(n*3, 0, 1), ncol=3)
      A <- rbinom(n, size = 1, prob = pi0(W))
      if(null.sims) {
        Y <- rbinom(n, size = 1, prob = mu0.null(A, W))
        psi0 <- mean((mu0.null(1,W) - mu0.null(0,W))^2)
        theta0 <- var((mu0.null(1,W) - mu0.null(0,W)))
      } else {
        Y <- rbinom(n, size = 1, prob = mu0(A, W))
        psi0 <- mean((mu0(1,W) - mu0(0,W))^2)
        theta0 <- var((mu0(1,W) - mu0(0,W)))
      }

      prop.reg <- gam(A ~ s(W[,1]) + s(W[,2]) + s(W[,3]), family = 'binomial')
      pi.hat <- prop.reg$fitted.values
      AW <- cbind(A, data.frame(W))
      mu.reg <- glm(Y ~ ., data=AW, family='binomial')
      mu.hat <- mu.reg$fitted.values
      mu1.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 1),
                                                 data.frame(W)), type = 'response')
      mu0.hat <- predict(mu.reg, newdata = cbind(data.frame(A = 0),
                                                 data.frame(W)), type = 'response')
      mu.hats <- data.frame(mu1=mu1.hat, mu0=mu0.hat)
      if (out.glm){
        control$pi.hat = pi.hat
        control$mu.hats = mu.hats
      }
      #print(control)
      ret <- htem.estimator(A, W, Y, control = control)

      ret$pi0 <- rep(mean(pi0(W)), 2)
      ret$pi.hat <- rep(mean(pi.hat), 2)
      ret$mu0 <- rep(mean(mu0(A, W)), 2)
      ret$mu.hat <- rep(mean(A*mu.hats$mu1 + (1-A)*mu.hats$mu0), 2)
      ret$n <- rep(n, 2)
      ret$j <- rep(j, 2)
      ret$seed <- rep(seed, 2)

      ret$psi0 <- rep(psi0, 2)
      ret$theta0 <- rep(theta0, 2)

      return(ret)
    })
  })

  return(ests)
}

## case 1
## Wald-type
ests.sim.1 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE), null.sims=FALSE, out.glm=TRUE)
## Bootstrap
ests.sim.2 <- ests.sim(200*c(1:10), 1:500, control = list(conf.int = TRUE, conf.int.type = 'boot'), null.sims=FALSE, out.glm=TRUE)

ests.sim.1[ests.sim.1$type == 'psi.est', c('est', 'se',
                                           'pi0', 'pi.hat',
                                           'mu0', 'mu.hat', 'n',
                                           'psi0')]

ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((pi0 - pi.hat), color=type)) +
  facet_wrap(~n, scales='free')

ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((mu0 - mu.hat), color=type)) +
  facet_wrap(~n, scales='free')

sum(ests.sim.1$est<0)
p <- ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes(sqrt(n) * (est - psi0), color=type)) +
  facet_wrap(~n, scales='free')

ggplot(psi.summaries.1) +
  geom_line(aes(n,  bias, color=type))

ggplot(subset(ests.sim.1, (type %in% 'psi.est'))) +
  geom_density(aes((pi0), color=type)) +
  facet_wrap(~n, scales='free')

ests.sim.1[ests.sim.1$type == 'psi.est', ]$ll

ggplot(subset(ests.sim.1, (type %in% 'psi.est') & (n %in% c(200, 600, 1000, 2000)) &
                (psi.coverage ==FALSE))) +
  geom_linerange(mapping = aes(x=j, ymin=ll, ymax=ul)) +
  geom_point(mapping=aes(x=j, y=psi0), size=2, col="red") +
  facet_wrap(~n, scales='free')


ests.sim.1$psi.coverage
ests.sim.1$psi.coverage <- (ests.sim.1$ll <= ests.sim.1$psi0) & (ests.sim.1$psi0 <= ests.sim.1$ul)
ests.sim.1$theta.coverage <- (ests.sim.1$ll <= ests.sim.1$theta0) & (ests.sim.1$theta0 <= ests.sim.1$ul)


ggplot() +
  geom_line(data = psi.summaries.1, aes(n, coverage, colour='wald-type')) +
  geom_line(data = psi.summaries.2, aes(n, coverage, colour='bootstrap')) +
  geom_hline(yintercept=.95)

load("ests.sim.1.Rdata")
# analysis
psi.summaries.1 <- ddply(subset(ests.sim.1, (type %in% 'psi.est')), .(n, type), summarize,
                   na = sum(is.na(est)),
                   coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                   bias = mean(est - psi0, na.rm=TRUE),
                   var = var(est, na.rm=TRUE),
                   mse = mean((est - psi0)^2, na.rm=TRUE))
library(ggplot2)
p1 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

p2 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

p3 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * mse, color=type))

p4 <- ggplot(psi.summaries.1) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)

psi.summaries.2 <- ddply(subset(ests.sim.2, (type %in% 'psi.est')), .(n, type), summarize,
                       na = sum(is.na(est)),
                       coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                       bias = mean(est - psi0, na.rm=TRUE),
                       var = var(est, na.rm=TRUE),
                       mse = mean((est - psi0)^2, na.rm=TRUE))

theta.summaries <- ddply(subset(ests, (type %in% 'theta.est')), .(n, type), summarize,
                       na = sum(is.na(est)),
                       coverage = mean(ll <= theta0 & theta0 <= ul, na.rm=TRUE),
                       bias = mean(est - theta0, na.rm=TRUE),
                       var = var(est, na.rm=TRUE),
                       mse = mean((est - theta0)^2, na.rm=TRUE))

#

# simulation 1: SL.gam, SL.gam
ests.sim.1 <-  est.psi.sim(200*c(1:10), 1:500, func_1 = "SL.gam", func_2 = "SL.gam", null.sims=FALSE)
load("ests.sim.1.RData")
est.psi.plot(ests.sim.1, plot.type='density')
est.theta.plot(ests.sim.1, plot.type='density')
summaries.psi <- est.psi.summary(ests.sim.1)
summaries.theta <- est.theta.summary(ests.sim.1)

save(ests.sim.1, file = "ests.sim.1.RData")

save(ests.sim.2, file = "ests.sim.2.RData")
save(ests.sim.1.null, file = "ests.sim.1.null.RData")

library(ggplot2)
ggplot(summaries.psi) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

ggplot(summaries.psi) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * mse, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

# simulation 1.1: SL.gam, SL.gam, null
ests.sim.1.null <-  est.psi.sim(200*c(1:10), 1:500, func_1 = "SL.gam", func_2 = "SL.gam", null.sims=TRUE)
est.psi.plot(ests.sim.1.null, plot.type='density')
est.theta.plot(ests.sim.1.null, plot.type='density')
summaries.psi <- est.psi.summary(ests.sim.1.null)
summaries.theta <- est.theta.summary(ests.sim.1.null)

library(ggplot2)
ggplot(summaries.psi) +
  geom_line(aes(n, coverage, color=type)) +
  geom_hline(yintercept=.95)

ggplot(summaries.psi) +
  geom_line(aes(n, sqrt(n) * bias, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * mse, color=type))

ggplot(summaries.psi) +
  geom_line(aes(n, n * var, color=type)) +
  coord_cartesian(ylim=c(0,1))

# simulation 2: SL.gam, SL.gam
ests.sim.2 <-  est.psi.sim(200*c(1:10), 1:500, func_1 = "SL.glm", func_2 = "SL.glm", null.sims=FALSE)
est.psi.plot(ests.sim.2, plot.type='density')
est.theta.plot(ests.sim.2, plot.type='density')
summaries.psi <- est.psi.summary(ests.sim.2)
summaries.theta <- est.theta.summary(ests.sim.2)

# simulation 3: SL.earth, SL.earth
ests.sim.3 <-  est.psi.sim(200*c(1:10), 1:500, func_1 = "SL.earth", func_2 = "SL.earth", null.sims=FALSE)
est.psi.plot(ests.sim.3, plot.type='density')
est.theta.plot(ests.sim.3, plot.type='density')
summaries.psi <- est.psi.summary(ests.sim.3)
summaries.theta <- est.theta.summary(ests.sim.3)

save(ests.sim.3, file = "ests.sim.3.RData")

# true data for simulation
n <- 200
W <- matrix(runif(n*3, 0, 1), ncol=3)
plot(pi0(W), ylim=c(0,1), col=1, pch=20)
for (i in 2:6){
  W <- matrix(runif(n*3, 0, 1), ncol=3)
  points(pi0(W), ylim=c(0,1), col=i, pch=20)
}

W <- matrix(runif(n*3, 0, 1), ncol=3)
A <- rbinom(n, size = 1, prob = pi0(W))
plot(mu0(A, W), ylim=c(0,1), col=1, pch=20)
for (i in 2:6){
  W <- matrix(runif(n*3, 0, 1), ncol=3)
  A <- rbinom(n, size = 1, prob = pi0(W))
  points(pi0(W), ylim=c(0,1), col=i, pch=20)
}
